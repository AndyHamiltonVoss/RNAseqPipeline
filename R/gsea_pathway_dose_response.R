################################################################################
# RNAseqPipeline
#
# File: gsea_pathway_dose_response.R
#
# Description:
#   For GSEA GO/KEGG results, some pathways show compounds split in
#   OPPOSITE directions -- some up-regulating, others down-regulating the
#   same pathway (visible as mixed red/blue rows in the top-variable-pathway
#   NES heatmap in R/gsea_compound_correlation.R). This asks two further
#   questions about that split:
#     - Which pathways show the strongest such split, and is there a
#       clearer way to see it than one row of a big heatmap?
#     - Is a compound's up/down direction for a pathway DOSE-DEPENDENT --
#       does NES trend consistently with dose, or is it noisy/flat?
#
#   Answered with two views per pathway family -- a line plot (NES vs.
#   dose, one line per compound, faceted by pathway) and a heatmap
#   (pathway x dose, one small block per compound) -- for a selected set
#   of pathways. selection_method picks which pathways:
#     - "bidirectional" (default): the strongest cross-compound NES SIGN
#       SPLIT at top dose (max(NES) * min(NES) most negative) -- pathways
#       where some compounds go up and others go down.
#     - "variance": the highest cross-compound NES VARIANCE at top dose,
#       same criterion as the pathway x compound NES heatmap in
#       R/gsea_compound_correlation.R, with no sign-crossing requirement
#       -- a pathway every compound pushes hard in the SAME direction can
#       rank highly here even though it wouldn't under "bidirectional".
#   Mirrors R/dose_response.R's dose-axis/color/line styling, but at the
#   pathway-NES level instead of gene-count level, and reuses
#   parse_dose_group() from there.
#
#   gsea_dir selects the data source (same parameterization as
#   R/gsea_compound_correlation.R):
#     - "GSEA_LFC" (default argument used in DGE_pipeline.R): already
#       complete at every dose for every compound, but significant-only --
#       a dose where a pathway isn't significant simply has no row for it,
#       which shows as a gap in the line (not interpolated across), which
#       is itself informative (the dose where a pathway "turns on").
#     - "GSEA_LFC_Dense": every tested pathway at every dose (once a full
#       all-dose dense batch run has been done -- see project plan), no
#       gaps. Same code, just point gsea_dir there once available.
#
################################################################################

#===============================================================================
# Build an Extreme-Dose (Top or Lowest) Compound x Pathway NES Matrix
#===============================================================================
#
# Same dose discovery + compound x pathway NES matrix build as
# run_gsea_pathway_family_correlation() in R/gsea_compound_correlation.R
# (duplicated here rather than imported, consistent with this codebase's
# one-file-per-view convention), generalized to pick either each
# compound's highest or lowest tested dose. Shared by every
# pathway-selection criterion below.
#===============================================================================

build_extreme_dose_nes_matrix <- function(config, pathway_family, gsea_dir, dose_selector = c("max", "min")) {

    dose_selector <- match.arg(dose_selector)
    pick_dose_row <- if (dose_selector == "max") which.max else which.min

    table_dir <- file.path(config$output_directory, gsea_dir, pathway_family, "Tables")

    result_files <- list.files(
        table_dir,
        pattern = paste0("^GSEA_", pathway_family, "_.*\\.csv$")
    )

    group_names <- sub(
        "\\.csv$", "",
        sub(paste0("^GSEA_", pathway_family, "_"), "", result_files)
    )

    parsed <- lapply(group_names, parse_dose_group)

    keep <- !vapply(parsed, is.null, logical(1))

    dose_table <- data.frame(
        group = group_names[keep],
        compound = vapply(parsed[keep], `[[`, character(1), "compound"),
        dose_nM = vapply(parsed[keep], `[[`, numeric(1), "dose_nM")
    )

    extreme_dose_table <- do.call(
        rbind,
        lapply(split(dose_table, dose_table$compound), function(d) {
            d[pick_dose_row(d$dose_nM), ]
        })
    )

    compounds <- sort(extreme_dose_table$compound)

    per_compound_nes <- lapply(compounds, function(cmpd) {
        group <- extreme_dose_table$group[extreme_dose_table$compound == cmpd]
        results <- utils::read.csv(file.path(table_dir, paste0("GSEA_", pathway_family, "_", group, ".csv")))
        stats::setNames(results$NES, results$ID)
    })

    names(per_compound_nes) <- compounds

    nes_matrix <- Reduce(
        function(a, b) merge(a, b, by = "ID", all = TRUE),
        lapply(compounds, function(cmpd) {
            df <- data.frame(ID = names(per_compound_nes[[cmpd]]), stringsAsFactors = FALSE)
            df[[cmpd]] <- as.numeric(per_compound_nes[[cmpd]])
            df
        })
    )

    pathway_ids <- nes_matrix$ID
    nes_matrix$ID <- NULL
    nes_matrix <- as.matrix(nes_matrix[, compounds, drop = FALSE])
    rownames(nes_matrix) <- pathway_ids

    list(nes_matrix = nes_matrix, compounds = compounds, extreme_dose_table = extreme_dose_table, table_dir = table_dir)

}

#===============================================================================
# Build the Top-Dose / Lowest-Dose Compound x Pathway NES Matrix
#===============================================================================
#
# Thin wrappers over build_extreme_dose_nes_matrix() kept for readability
# at call sites; top_dose_table is just extreme_dose_table under its
# original name, for backward compatibility with existing callers.
#===============================================================================

build_top_dose_nes_matrix <- function(config, pathway_family, gsea_dir) {

    built <- build_extreme_dose_nes_matrix(config, pathway_family, gsea_dir, dose_selector = "max")
    built$top_dose_table <- built$extreme_dose_table
    built

}

build_lowest_dose_nes_matrix <- function(config, pathway_family, gsea_dir) {

    build_extreme_dose_nes_matrix(config, pathway_family, gsea_dir, dose_selector = "min")

}

#===============================================================================
# Look Up Pathway Descriptions for a Set of Selected IDs
#===============================================================================

lookup_pathway_descriptions <- function(table_dir, pathway_family, top_dose_table, compounds, selected_ids) {

    description_lookup <- character(0)

    for (cmpd in compounds) {
        group <- top_dose_table$group[top_dose_table$compound == cmpd]
        results <- utils::read.csv(file.path(table_dir, paste0("GSEA_", pathway_family, "_", group, ".csv")))
        new_ids <- setdiff(intersect(selected_ids, results$ID), names(description_lookup))
        if (length(new_ids) > 0) {
            description_lookup[new_ids] <- stats::setNames(results$Description, results$ID)[new_ids]
        }
    }

    unname(description_lookup[selected_ids])

}

#===============================================================================
# Select the Pathways With the Strongest Cross-Compound Sign Split
#===============================================================================
#
# Ranks by max(NES) * min(NES) across compounds at top dose -- the more
# negative, the stronger the split (one compound strongly positive,
# another strongly negative). This is a different criterion from
# select_top_variance_pathways() below (or R/gsea_compound_correlation.R's
# "top variance" heatmap selection), which doesn't guarantee a sign
# crossing -- a pathway that's high-magnitude but same-direction in every
# compound has a high variance but a POSITIVE (not negative) product, so
# it would NOT be selected here.
#===============================================================================

select_bidirectional_pathways <- function(config, pathway_family, gsea_dir, top_n) {

    built <- build_top_dose_nes_matrix(config, pathway_family, gsea_dir)
    nes_matrix <- built$nes_matrix

    max_nes <- apply(nes_matrix, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    min_nes <- apply(nes_matrix, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))

    bidirectional_score <- max_nes * min_nes

    ranked_ids <- names(sort(bidirectional_score, na.last = NA))

    n_select <- min(top_n, length(ranked_ids))
    selected_ids <- ranked_ids[seq_len(n_select)]

    data.frame(
        pathway_ID = selected_ids,
        Description = lookup_pathway_descriptions(built$table_dir, pathway_family, built$top_dose_table, built$compounds, selected_ids),
        rank_score = bidirectional_score[selected_ids],
        max_NES = max_nes[selected_ids],
        min_NES = min_nes[selected_ids],
        stringsAsFactors = FALSE
    )

}

#===============================================================================
# Select the Pathways With the Highest Cross-Compound NES Variance
#===============================================================================
#
# Ranks by var(NES) across compounds at top dose, descending -- the same
# criterion as the pathway x compound NES heatmap in
# R/gsea_compound_correlation.R, but applied here to feed the dose-level
# line/heatmap views instead of a single-dose snapshot. Unlike
# select_bidirectional_pathways(), does not require a sign crossing --
# a pathway that every compound pushes strongly in the SAME direction
# (e.g. all compounds +2 to +3) can still rank highly here.
#===============================================================================

select_top_variance_pathways <- function(config, pathway_family, gsea_dir, top_n) {

    built <- build_top_dose_nes_matrix(config, pathway_family, gsea_dir)
    nes_matrix <- built$nes_matrix

    pathway_vars <- apply(nes_matrix, 1, stats::var, na.rm = TRUE)
    pathway_vars <- pathway_vars[!is.na(pathway_vars) & is.finite(pathway_vars)]

    max_nes <- apply(nes_matrix, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    min_nes <- apply(nes_matrix, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))

    ranked_ids <- names(sort(pathway_vars, decreasing = TRUE))

    n_select <- min(top_n, length(ranked_ids))
    selected_ids <- ranked_ids[seq_len(n_select)]

    data.frame(
        pathway_ID = selected_ids,
        Description = lookup_pathway_descriptions(built$table_dir, pathway_family, built$top_dose_table, built$compounds, selected_ids),
        rank_score = pathway_vars[selected_ids],
        max_NES = max_nes[selected_ids],
        min_NES = min_nes[selected_ids],
        stringsAsFactors = FALSE
    )

}

#===============================================================================
# Select Pathways With a "Clean" Dose-Response Divergence
#===============================================================================
#
# Specifically targets pathways where compounds start out CLUSTERED at a
# similar NES at their lowest tested dose, then split into clearly
# up-regulated vs. down-regulated groups by their highest tested dose --
# a "fan out" pattern, as opposed to select_bidirectional_pathways()
# (which only looks at top dose, so would equally rank a pathway that was
# ALREADY split at the lowest dose and stayed split).
#
# Two requirements, both at each compound's own extreme doses:
#   - A genuine sign crossing at TOP dose (max(NES) * min(NES) < 0 there)
#     -- same bidirectional test as select_bidirectional_pathways(), used
#     here as a hard filter, not just part of the score.
#   - Ranked by range_at_lowest_dose - range_at_top_dose (more negative
#     = tighter clustering at the lowest dose AND a wider split at the
#     top dose -- the "nicest" fan-out shape).
#===============================================================================

select_dose_response_divergence_pathways <- function(config, pathway_family, gsea_dir, top_n) {

    built_top <- build_top_dose_nes_matrix(config, pathway_family, gsea_dir)
    built_low <- build_lowest_dose_nes_matrix(config, pathway_family, gsea_dir)

    common_ids <- intersect(rownames(built_top$nes_matrix), rownames(built_low$nes_matrix))
    common_compounds <- intersect(built_top$compounds, built_low$compounds)

    nes_top <- built_top$nes_matrix[common_ids, common_compounds, drop = FALSE]
    nes_low <- built_low$nes_matrix[common_ids, common_compounds, drop = FALSE]

    range_or_na <- function(x) if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
    extreme_or_na <- function(x, fn) if (all(is.na(x))) NA_real_ else fn(x, na.rm = TRUE)

    max_top <- apply(nes_top, 1, extreme_or_na, fn = max)
    min_top <- apply(nes_top, 1, extreme_or_na, fn = min)
    range_top <- max_top - min_top
    range_low <- apply(nes_low, 1, range_or_na)

    crosses_zero_at_top <- !is.na(max_top) & !is.na(min_top) & (max_top * min_top) < 0

    divergence_score <- range_low - range_top
    divergence_score[!crosses_zero_at_top | is.na(range_low) | is.na(range_top)] <- NA_real_

    ranked_ids <- names(sort(divergence_score, na.last = NA))

    n_select <- min(top_n, length(ranked_ids))
    selected_ids <- ranked_ids[seq_len(n_select)]

    data.frame(
        pathway_ID = selected_ids,
        Description = lookup_pathway_descriptions(built_top$table_dir, pathway_family, built_top$top_dose_table, built_top$compounds, selected_ids),
        rank_score = divergence_score[selected_ids],
        range_at_lowest_dose = range_low[selected_ids],
        range_at_top_dose = range_top[selected_ids],
        max_NES = max_top[selected_ids],
        min_NES = min_top[selected_ids],
        stringsAsFactors = FALSE
    )

}

#===============================================================================
# Build a Long Table of NES vs. Dose for Selected Pathways
#===============================================================================
#
# Unlike run_gsea_pathway_family_correlation(), keeps EVERY dose level per
# compound (not collapsed to top dose).
#===============================================================================

build_pathway_dose_nes_table <- function(config, pathway_family, gsea_dir, pathway_ids) {

    table_dir <- file.path(config$output_directory, gsea_dir, pathway_family, "Tables")

    result_files <- list.files(
        table_dir,
        pattern = paste0("^GSEA_", pathway_family, "_.*\\.csv$")
    )

    group_names <- sub(
        "\\.csv$", "",
        sub(paste0("^GSEA_", pathway_family, "_"), "", result_files)
    )

    parsed <- lapply(group_names, parse_dose_group)

    keep <- !vapply(parsed, is.null, logical(1))

    dose_table <- data.frame(
        group = group_names[keep],
        compound = vapply(parsed[keep], `[[`, character(1), "compound"),
        dose_nM = vapply(parsed[keep], `[[`, numeric(1), "dose_nM")
    )

    rows <- lapply(seq_len(nrow(dose_table)), function(i) {

        results <- utils::read.csv(file.path(table_dir, paste0("GSEA_", pathway_family, "_", dose_table$group[i], ".csv")))

        matched <- results[match(pathway_ids, results$ID), c("ID", "Description", "NES")]

        data.frame(
            compound = dose_table$compound[i],
            dose_nM = dose_table$dose_nM[i],
            group = dose_table$group[i],
            pathway_ID = pathway_ids,
            Description = matched$Description,
            NES = matched$NES,
            stringsAsFactors = FALSE
        )

    })

    nes_table <- do.call(rbind, rows)

    # Description can be NA for a dose where the pathway wasn't
    # significant (significant-only directory) -- backfill from any dose
    # that does have it, so facet labels are never blank.
    description_lookup <- stats::setNames(
        rep(NA_character_, length(pathway_ids)),
        pathway_ids
    )

    known <- !is.na(nes_table$Description)

    description_lookup[nes_table$pathway_ID[known]] <- nes_table$Description[known]

    nes_table$Description <- description_lookup[nes_table$pathway_ID]

    nes_table

}

#===============================================================================
# Plot NES vs. Dose, Faceted by Pathway
#===============================================================================

plot_pathway_dose_response <- function(nes_table, pathway_family, gsea_dir, config, selection_label = "Most Bidirectional Pathways") {

    compounds <- sort(unique(nes_table$compound))

    base_colors <- c(Rosi = "#1b9e77", MRL24 = "#d95f02")

    unmatched <- setdiff(compounds, names(base_colors))

    if (length(unmatched) > 0) {

        extra_palette <- grDevices::colorRampPalette(
            c("#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#1f78b4", "#33a02c", "#e31a1c")
        )(length(unmatched))

        compound_colors <- c(base_colors, stats::setNames(extra_palette, unmatched))

    } else {

        compound_colors <- base_colors

    }

    compound_colors <- compound_colors[compounds]

    gsea_label <- switch(
        gsea_dir,
        "GSEA" = "GSEA",
        "GSEA_LFC" = "GSEA (LFC-Ranked)",
        "GSEA_Dense" = "GSEA (Dense, All Tested Pathways)",
        "GSEA_LFC_Dense" = "GSEA (LFC-Ranked, Dense)"
    )

    n_pathways <- length(unique(nes_table$pathway_ID))
    n_cols <- min(5, n_pathways)
    n_rows <- ceiling(n_pathways / n_cols)

    plot <- ggplot2::ggplot(
        nes_table,
        ggplot2::aes(x = dose_nM, y = NES, color = compound, group = compound)
    ) +
        ggplot2::geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
        ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
        ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_color_manual(values = compound_colors) +
        ggplot2::facet_wrap(~Description, scales = "free_y", ncol = n_cols) +
        ggplot2::labs(
            title = paste0(gsea_label, " ", pathway_family, " Pathway NES vs. Dose (", selection_label, ")"),
            x = "Dose (nM, log10 scale)",
            y = "NES",
            color = "Compound"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(strip.text = ggplot2::element_text(size = 7))

    attr(plot, "n_rows") <- n_rows
    attr(plot, "n_cols") <- n_cols

    plot

}

#===============================================================================
# Plot NES Heatmap: Pathways x Dose, One Facet Block Per Compound
#===============================================================================
#
# Same underlying nes_table as plot_pathway_dose_response(), shown as a
# grid instead of lines: y = pathway, x = dose, one small heatmap "block"
# per compound (facet_wrap(~compound)). Dose position is ordinal (1st
# through 4th lowest-to-highest dose tested FOR THAT COMPOUND, since
# compounds differ in both absolute dose range and number of dose levels
# tested -- 3 for SR1664/SR16832, 4 for everyone else) rather than a
# shared nM axis, so every compound's block is directly comparable by
# dose rank; the actual nM value is still in the saved CSV (dose_nM
# column) for reference. Pathways are ordered top-to-bottom by rank
# (strongest-ranked first, per whichever selection criterion produced
# `selected`), shared across all blocks.
#===============================================================================

plot_pathway_dose_heatmap <- function(nes_table, selected, pathway_family, gsea_dir, config, selection_label = "Most Bidirectional Pathways") {

    compounds <- sort(unique(nes_table$compound))

    nes_table$dose_rank <- stats::ave(
        nes_table$dose_nM, nes_table$compound,
        FUN = function(x) match(x, sort(unique(x)))
    )

    max_dose_rank <- max(nes_table$dose_rank)

    nes_table$dose_label <- factor(
        paste0("Dose ", nes_table$dose_rank),
        levels = paste0("Dose ", seq_len(max_dose_rank))
    )

    # Top-to-bottom reading order = strongest-ranked pathway first;
    # `selected` is already sorted that way by whichever select_*()
    # function produced it. ggplot's discrete y-axis plots the first
    # factor level at the BOTTOM by default, so the level order is
    # reversed here to get that.
    pathway_order <- selected$Description
    nes_table$Description <- factor(nes_table$Description, levels = rev(pathway_order))

    gsea_label <- switch(
        gsea_dir,
        "GSEA" = "GSEA",
        "GSEA_LFC" = "GSEA (LFC-Ranked)",
        "GSEA_Dense" = "GSEA (Dense, All Tested Pathways)",
        "GSEA_LFC_Dense" = "GSEA (LFC-Ranked, Dense)"
    )

    max_abs_nes <- max(abs(nes_table$NES), na.rm = TRUE)

    n_cols <- min(5, length(compounds))
    n_facet_rows <- ceiling(length(compounds) / n_cols)
    n_pathways <- length(pathway_order)

    plot <- ggplot2::ggplot(
        nes_table,
        ggplot2::aes(x = dose_label, y = Description, fill = NES)
    ) +
        ggplot2::geom_tile(color = "grey90") +
        ggplot2::scale_fill_gradient2(
            low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
            limits = c(-max_abs_nes, max_abs_nes), na.value = "grey85"
        ) +
        ggplot2::facet_wrap(~compound, ncol = n_cols) +
        ggplot2::labs(
            title = paste0(gsea_label, " ", pathway_family, " Pathway NES by Dose (", selection_label, ")"),
            x = "Dose (lowest -> highest tested)",
            y = NULL,
            fill = "NES"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            strip.text = ggplot2::element_text(size = 8, face = "bold"),
            axis.text.x = ggplot2::element_text(size = 6, angle = 45, hjust = 1),
            axis.text.y = ggplot2::element_text(size = 6.5),
            panel.grid = ggplot2::element_blank()
        )

    attr(plot, "n_facet_rows") <- n_facet_rows
    attr(plot, "n_cols") <- n_cols
    attr(plot, "n_pathways") <- n_pathways

    plot

}

#===============================================================================
# Run Pathway Dose-Response for One Pathway Family (GO/KEGG)
#===============================================================================

run_gsea_pathway_dose_response <- function(config, pathway_family, gsea_dir = "GSEA_LFC", selection_method = c("bidirectional", "variance", "divergence")) {

    selection_method <- match.arg(selection_method)

    method_info <- switch(
        selection_method,
        "bidirectional" = list(
            select_fn = select_bidirectional_pathways,
            top_n = config$gsea_bidirectional_top_pathways,
            label = "Most Bidirectional Pathways",
            file_tag = "bidirectional_pathways"
        ),
        "variance" = list(
            select_fn = select_top_variance_pathways,
            top_n = config$gsea_nes_heatmap_top_pathways,
            label = "Top Variance Pathways",
            file_tag = "top_variance_pathways"
        ),
        "divergence" = list(
            select_fn = select_dose_response_divergence_pathways,
            top_n = config$gsea_divergence_top_pathways,
            label = "Clean Dose-Response Divergence: Similar at Lowest Dose, Split at Top Dose",
            file_tag = "dose_divergence_pathways"
        )
    )

    message(sprintf("Selecting %s %s %s pathways...", selection_method, gsea_dir, pathway_family))

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    selected <- method_info$select_fn(config, pathway_family, gsea_dir, top_n = method_info$top_n)

    message(sprintf(
        "Building dose-response NES table for %d %s pathways...",
        nrow(selected), pathway_family
    ))

    nes_table <- build_pathway_dose_nes_table(config, pathway_family, gsea_dir, selected$pathway_ID)

    plot <- plot_pathway_dose_response(nes_table, pathway_family, gsea_dir, config, selection_label = method_info$label)

    table_dir <- file.path(config$output_directory, gsea_dir, pathway_family, "Tables")
    figure_dir <- file.path(config$output_directory, gsea_dir, pathway_family, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        selected,
        file.path(table_dir, paste0("compound_gsea_", pathway_family, "_", method_info$file_tag, ".csv")),
        row.names = FALSE
    )

    save_csv(
        nes_table,
        file.path(table_dir, paste0("compound_gsea_", pathway_family, "_", method_info$file_tag, "_dose_nes.csv")),
        row.names = FALSE
    )

    n_rows <- attr(plot, "n_rows")
    n_cols <- attr(plot, "n_cols")

    plot_width <- max(config$figure_width, n_cols * 2.8)
    plot_height <- max(config$figure_height, n_rows * 2.3)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0("GSEA_", pathway_family, "_", method_info$file_tag, "_dose_response.png")),
            plot = plot, width = plot_width, height = plot_height, dpi = config$figure_dpi, limitsize = FALSE
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0("GSEA_", pathway_family, "_", method_info$file_tag, "_dose_response.pdf")),
            plot = plot, width = plot_width, height = plot_height, limitsize = FALSE
        )

    }

    heatmap_plot <- plot_pathway_dose_heatmap(nes_table, selected, pathway_family, gsea_dir, config, selection_label = method_info$label)

    n_facet_rows <- attr(heatmap_plot, "n_facet_rows")
    n_heatmap_cols <- attr(heatmap_plot, "n_cols")
    n_pathways <- attr(heatmap_plot, "n_pathways")

    heatmap_width <- max(config$figure_width, n_heatmap_cols * 2.2)
    heatmap_height <- max(config$figure_height, n_facet_rows * (2 + n_pathways * 0.18))

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0("GSEA_", pathway_family, "_", method_info$file_tag, "_dose_heatmap.png")),
            plot = heatmap_plot, width = heatmap_width, height = heatmap_height, dpi = config$figure_dpi, limitsize = FALSE
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0("GSEA_", pathway_family, "_", method_info$file_tag, "_dose_heatmap.pdf")),
            plot = heatmap_plot, width = heatmap_width, height = heatmap_height, limitsize = FALSE
        )

    }

    message(sprintf("%s pathway dose-response complete.", pathway_family))

    list(selected = selected, nes_table = nes_table, plot = plot, heatmap_plot = heatmap_plot)

}

#===============================================================================
# Run Pathway Dose-Response (GO + KEGG)
#===============================================================================

run_gsea_pathway_dose_response_all <- function(config, gsea_dir = "GSEA_LFC", selection_method = c("bidirectional", "variance", "divergence")) {

    selection_method <- match.arg(selection_method)

    message(sprintf("Running GSEA pathway dose-response (%s selection)...", selection_method))

    go_result <- run_gsea_pathway_dose_response(config, "GO", gsea_dir = gsea_dir, selection_method = selection_method)
    kegg_result <- run_gsea_pathway_dose_response(config, "KEGG", gsea_dir = gsea_dir, selection_method = selection_method)

    message("GSEA pathway dose-response complete.")

    list(GO = go_result, KEGG = kegg_result)

}
