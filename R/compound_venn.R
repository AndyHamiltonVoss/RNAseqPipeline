################################################################################
# RNAseqPipeline
#
# File: compound_venn.R
#
# Description:
#   Venn diagram comparing the consensus regulated gene set between two
#   compounds, each taken at its own peak-regulation dose. Reuses
#   classify_genes() / load_method_results() from dose_response.R and
#   save_venn_plot() from venn_diagrams.R.
#
################################################################################

#===============================================================================
# Get Consensus Regulated Genes for One Contrast
#===============================================================================

get_consensus_genes <- function(group, direction, config, log2fc_cutoff = config$log2fc_cutoff) {

    if (identical(direction, "combined")) {
        return(
            union(
                get_consensus_genes(group, "up", config, log2fc_cutoff),
                get_consensus_genes(group, "down", config, log2fc_cutoff)
            )
        )
    }

    results <- load_method_results(group, config$output_directory, config$reference_group)

    empty <- list(up = character(0), down = character(0))

    classified <- list(
        DESeq2 = if (!is.null(results$DESeq2))
            classify_genes(results$DESeq2, "gene", "log2FoldChange", "padj", config$alpha, log2fc_cutoff)
        else empty,

        edgeR = if (!is.null(results$edgeR))
            classify_genes(results$edgeR, "gene", "logFC", "FDR", config$alpha, log2fc_cutoff)
        else empty,

        limma = if (!is.null(results$limma))
            classify_genes(results$limma, "gene", "logFC", "adj.P.Val", config$alpha, log2fc_cutoff)
        else empty
    )

    Reduce(
        intersect,
        list(classified$DESeq2[[direction]], classified$edgeR[[direction]], classified$limma[[direction]])
    )

}

#===============================================================================
# Find the Dose With the Most Consensus Regulated Genes for a Compound
#===============================================================================

find_peak_dose <- function(compound, groups, direction, config) {

    series <- get_dose_series(groups, compound)

    if (nrow(series) == 0) {
        stop("No dose series found for compound: ", compound)
    }

    series$n_genes <- vapply(
        series$group,
        function(group) length(get_consensus_genes(group, direction, config)),
        integer(1)
    )

    peak <- series[which.max(series$n_genes), ]

    return(peak)

}

#===============================================================================
# Format a Direction Value Into a Plot-Title Label
#===============================================================================

format_direction_label <- function(direction) {

    if (identical(direction, "combined")) "up+down regulated" else paste0(direction, "regulated")

}

#===============================================================================
# Build a Two-Way Venn Diagram Grob
#===============================================================================

build_two_way_venn <- function(gene_sets, title, margin = 0.1) {

    # VennDiagram::draw.pairwise.venn always redraws with the larger set at
    # its internal "position 1", reassigning category/fill/count labels to
    # match -- so which set lands left/right (and which count follows which
    # label) depends on set size, not input order, and fighting that via its
    # `inverted` flag desyncs the count text from the category labels. A
    # small custom two-circle plot avoids the issue entirely: position,
    # color, and count placement are all fixed by input order here.

    category_names <- names(gene_sets)

    n_a_only <- length(setdiff(gene_sets[[1]], gene_sets[[2]]))
    n_b_only <- length(setdiff(gene_sets[[2]], gene_sets[[1]]))
    n_shared <- length(intersect(gene_sets[[1]], gene_sets[[2]]))

    n_a_total <- n_a_only + n_shared
    n_b_total <- n_b_only + n_shared

    # Circle radius scales with sqrt(gene count) so that circle *area* --
    # what the eye actually compares -- is proportional to set size, floored
    # so an empty/tiny set is still visible and ceilinged for a stable layout.
    r_min <- 0.7
    r_max <- 1.6

    scale_radius <- function(n_total, n_max) {
        if (n_max <= 0) return(r_min)
        r_min + (r_max - r_min) * sqrt(n_total / n_max)
    }

    n_max <- max(n_a_total, n_b_total, 1)

    r_a <- scale_radius(n_a_total, n_max)
    r_b <- scale_radius(n_b_total, n_max)

    make_circle <- function(x0, y0, r, group, n = 200) {
        theta <- seq(0, 2 * pi, length.out = n)
        data.frame(x = x0 + r * cos(theta), y = y0 + r * sin(theta), group = group)
    }

    offset <- 0.58 * (r_a + r_b) / 2

    circles <- rbind(
        make_circle(-offset, 0, r_a, category_names[1]),
        make_circle(offset, 0, r_b, category_names[2])
    )
    circles$group <- factor(circles$group, levels = category_names)

    fill_colors <- stats::setNames(c("#1b9e77", "#d95f02"), category_names)

    count_labels <- data.frame(
        x = c(-offset - r_a * 0.4, 0, offset + r_b * 0.4),
        y = 0,
        label = c(n_a_only, n_shared, n_b_only)
    )

    category_labels <- data.frame(
        x = c(-offset, offset),
        y = c(r_a, r_b) + 0.35,
        label = category_names
    )

    plot <- ggplot2::ggplot() +
        ggplot2::geom_polygon(
            data = circles,
            ggplot2::aes(x = x, y = y, group = group, fill = group),
            alpha = 0.5, color = "black", linewidth = 0.6
        ) +
        ggplot2::geom_text(
            data = count_labels,
            ggplot2::aes(x = x, y = y, label = label),
            size = 5
        ) +
        ggplot2::geom_text(
            data = category_labels,
            ggplot2::aes(x = x, y = y, label = label),
            size = 4.2
        ) +
        ggplot2::scale_fill_manual(values = fill_colors, guide = "none") +
        ggplot2::coord_fixed(clip = "off") +
        ggplot2::labs(title = title) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
            plot.margin = ggplot2::margin(20, 20, 20, 20)
        )

    return(plot)

}

#===============================================================================
# Save a Two-Way Venn Diagram (ggplot2)
#===============================================================================

save_two_way_venn_plot <- function(
    venn_plot,
    filename_base,
    figure_dir,
    config,
    width = config$figure_width,
    height = config$figure_height
) {

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_base, ".png")),
            plot = venn_plot,
            width = width,
            height = height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_base, ".pdf")),
            plot = venn_plot,
            width = width,
            height = height
        )

    }

    invisible(TRUE)

}

#===============================================================================
# Compare Two Compounds at Their Peak-Regulation Dose
#===============================================================================

run_compound_comparison <- function(
    compound_a,
    compound_b,
    config,
    metadata,
    direction = "up"
) {

    message("Running compound comparison: ", compound_a, " vs ", compound_b)

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    groups <- setdiff(
        levels(factor(metadata$group)),
        config$reference_group
    )

    peak_a <- find_peak_dose(compound_a, groups, direction, config)
    peak_b <- find_peak_dose(compound_b, groups, direction, config)

    message(
        compound_a, " peak dose: ", peak_a$group,
        " (", peak_a$n_genes, " ", direction, "regulated consensus genes)"
    )
    message(
        compound_b, " peak dose: ", peak_b$group,
        " (", peak_b$n_genes, " ", direction, "regulated consensus genes)"
    )

    genes_a <- get_consensus_genes(peak_a$group, direction, config)
    genes_b <- get_consensus_genes(peak_b$group, direction, config)

    gene_sets <- stats::setNames(
        list(genes_a, genes_b),
        c(peak_a$group, peak_b$group)
    )

    title <- paste0(
        peak_a$group, " vs ", peak_b$group,
        " (", direction, "regulated, consensus)"
    )

    venn_plot <- build_two_way_venn(gene_sets, title)

    figure_dir <- file.path(config$output_directory, "Figures", "Venn")

    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    filename_base <- paste0(
        "Venn_", compound_a, "_vs_", compound_b, "_", direction, "_peak_dose"
    )

    save_two_way_venn_plot(venn_plot, filename_base, figure_dir, config)

    only_a <- setdiff(genes_a, genes_b)
    only_b <- setdiff(genes_b, genes_a)
    shared <- intersect(genes_a, genes_b)

    overlap_summary <- data.frame(
        region = c(
            paste(peak_a$group, "only"),
            paste(peak_b$group, "only"),
            "Shared"
        ),
        n_genes = c(length(only_a), length(only_b), length(shared)),
        stringsAsFactors = FALSE
    )

    all_genes <- union(genes_a, genes_b)

    gene_membership <- data.frame(
        gene = all_genes,
        stringsAsFactors = FALSE
    )

    gene_membership[[peak_a$group]] <- gene_membership$gene %in% genes_a
    gene_membership[[peak_b$group]] <- gene_membership$gene %in% genes_b

    save_csv(
        overlap_summary,
        file.path(figure_dir, paste0(filename_base, "_overlap_counts.csv")),
        row.names = FALSE
    )

    save_csv(
        gene_membership,
        file.path(figure_dir, paste0(filename_base, "_gene_membership.csv")),
        row.names = FALSE
    )

    message("Compound comparison complete.")

    return(
        list(
            peak_a = peak_a,
            peak_b = peak_b,
            gene_sets = gene_sets,
            overlap_summary = overlap_summary,
            gene_membership = gene_membership
        )
    )

}

#===============================================================================
# Compare Two Groups Across a Range of log2FC Cutoffs for One Side
#===============================================================================

run_threshold_sweep_comparison <- function(
    group_a,
    group_b,
    config,
    log2fc_cutoff_b_values,
    log2fc_cutoff_a = config$log2fc_cutoff,
    direction = "up"
) {

    message(
        "Running threshold sweep: ", group_a, " (log2FC>", log2fc_cutoff_a, ") vs ",
        group_b, " across log2FC cutoffs: ", paste(log2fc_cutoff_b_values, collapse = ", ")
    )

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    genes_a <- get_consensus_genes(group_a, direction, config, log2fc_cutoff = log2fc_cutoff_a)

    figure_dir <- file.path(config$output_directory, "Figures", "Venn")

    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    summary_rows <- list()

    for (cutoff_b in log2fc_cutoff_b_values) {

        message("Processing threshold: ", group_b, " log2FC>", cutoff_b)

        genes_b <- get_consensus_genes(group_b, direction, config, log2fc_cutoff = cutoff_b)

        gene_sets <- stats::setNames(
            list(genes_a, genes_b),
            c(
                paste0(group_a, " (log2FC>", log2fc_cutoff_a, ")"),
                paste0(group_b, " (log2FC>", cutoff_b, ")")
            )
        )

        direction_label <- format_direction_label(direction)

        title <- paste0(
            group_a, " vs ", group_b,
            " (", direction_label, ", consensus)\n",
            group_a, ": log2FC>", log2fc_cutoff_a,
            "  |  ", group_b, ": log2FC>", cutoff_b
        )

        venn_plot <- build_two_way_venn(gene_sets, title, margin = 0.15)

        filename_base <- paste0(
            "Venn_", group_a, "_lfc", gsub("\\.", "", sprintf("%.2f", log2fc_cutoff_a)),
            "_vs_", group_b, "_lfc", gsub("\\.", "", sprintf("%.2f", cutoff_b)),
            "_", direction
        )

        save_two_way_venn_plot(
            venn_plot, filename_base, figure_dir, config,
            width = config$figure_width * 1.5,
            height = config$figure_height
        )

        only_a <- setdiff(genes_a, genes_b)
        only_b <- setdiff(genes_b, genes_a)
        shared <- intersect(genes_a, genes_b)

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
            group_a = group_a,
            log2fc_cutoff_a = log2fc_cutoff_a,
            group_b = group_b,
            log2fc_cutoff_b = cutoff_b,
            n_a_only = length(only_a),
            n_b_only = length(only_b),
            n_shared = length(shared),
            n_a_total = length(genes_a),
            n_b_total = length(genes_b),
            stringsAsFactors = FALSE
        )

        gene_membership <- data.frame(
            gene = union(genes_a, genes_b),
            stringsAsFactors = FALSE
        )

        gene_membership[[group_a]] <- gene_membership$gene %in% genes_a
        gene_membership[[group_b]] <- gene_membership$gene %in% genes_b

        save_csv(
            gene_membership,
            file.path(figure_dir, paste0(filename_base, "_gene_membership.csv")),
            row.names = FALSE
        )

    }

    summary_table <- do.call(rbind, summary_rows)

    save_csv(
        summary_table,
        file.path(
            figure_dir,
            paste0(
                "Venn_", group_a, "_lfc", gsub("\\.", "", sprintf("%.2f", log2fc_cutoff_a)),
                "_vs_", group_b, "_", direction,
                "_threshold_sweep_summary.csv"
            )
        ),
        row.names = FALSE
    )

    message("Threshold sweep comparison complete.")

    return(summary_table)

}

#===============================================================================
# Build a Four-Way (Dose Series) Venn Diagram Grob
#===============================================================================

build_dose_series_venn <- function(gene_sets, title) {

    # draw.quad.venn (unlike draw.pairwise.venn) keeps its four ellipses in
    # fixed, input-order positions -- there's no size-based left/right
    # reordering to fight here, so the VennDiagram package is safe to use
    # directly for this one.

    futile.logger::flog.threshold(
        futile.logger::ERROR,
        name = "VennDiagramLogger"
    )

    dose_labels <- names(gene_sets)

    fill_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

    grob <- VennDiagram::draw.quad.venn(
        area1 = length(gene_sets[[1]]),
        area2 = length(gene_sets[[2]]),
        area3 = length(gene_sets[[3]]),
        area4 = length(gene_sets[[4]]),
        n12 = length(intersect(gene_sets[[1]], gene_sets[[2]])),
        n13 = length(intersect(gene_sets[[1]], gene_sets[[3]])),
        n14 = length(intersect(gene_sets[[1]], gene_sets[[4]])),
        n23 = length(intersect(gene_sets[[2]], gene_sets[[3]])),
        n24 = length(intersect(gene_sets[[2]], gene_sets[[4]])),
        n34 = length(intersect(gene_sets[[3]], gene_sets[[4]])),
        n123 = length(Reduce(intersect, gene_sets[c(1, 2, 3)])),
        n124 = length(Reduce(intersect, gene_sets[c(1, 2, 4)])),
        n134 = length(Reduce(intersect, gene_sets[c(1, 3, 4)])),
        n234 = length(Reduce(intersect, gene_sets[c(2, 3, 4)])),
        n1234 = length(Reduce(intersect, gene_sets)),
        category = dose_labels,
        fill = fill_colors,
        alpha = 0.5,
        cex = 1.1,
        cat.cex = 1.1,
        cat.col = fill_colors,
        main = title,
        main.cex = 1.3,
        ind = FALSE
    )

    return(grob)

}

#===============================================================================
# Break Four Gene Sets Into All Non-Overlapping (Exclusive) Venn Regions
#===============================================================================

venn_regions_n <- function(gene_sets) {

    n <- length(gene_sets)
    set_names <- names(gene_sets)

    rows <- list()

    for (k in seq_len(n)) {

        combos <- utils::combn(seq_len(n), k, simplify = FALSE)

        for (idx in combos) {

            included <- Reduce(intersect, gene_sets[idx])

            outside_idx <- setdiff(seq_len(n), idx)

            exclusive <- if (length(outside_idx) > 0) {
                setdiff(included, Reduce(union, gene_sets[outside_idx]))
            } else {
                included
            }

            rows[[length(rows) + 1]] <- data.frame(
                region = paste(set_names[idx], collapse = " & "),
                n_genes = length(exclusive),
                stringsAsFactors = FALSE
            )

        }

    }

    do.call(rbind, rows)

}

#===============================================================================
# Compare a Compound's Own Dose Series (Four-Way Venn)
#===============================================================================

run_dose_series_venn <- function(compound, config, metadata, direction = "combined") {

    message("Running dose series Venn comparison: ", compound)

    if (!requireNamespace("VennDiagram", quietly = TRUE)) {
        stop("Package 'VennDiagram' is required.")
    }

    groups <- setdiff(
        levels(factor(metadata$group)),
        config$reference_group
    )

    series <- get_dose_series(groups, compound)

    if (nrow(series) != 4) {
        stop(
            "run_dose_series_venn() expects exactly 4 doses for compound '",
            compound, "', found ", nrow(series), ": ",
            paste(series$group, collapse = ", ")
        )
    }

    dose_labels <- sub(paste0("^", compound, "_"), "", series$group)

    gene_sets <- stats::setNames(
        lapply(series$group, function(g) get_consensus_genes(g, direction, config)),
        dose_labels
    )

    title <- paste0(
        compound, " Dose Series (", format_direction_label(direction), ", consensus)"
    )

    venn_grob <- build_dose_series_venn(gene_sets, title)

    figure_dir <- file.path(config$output_directory, "Figures", "Venn")

    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    filename_base <- paste0("Venn_", compound, "_dose_series_", direction)

    save_venn_plot(
        venn_grob, filename_base, figure_dir, config,
        width = config$figure_width, height = config$figure_height
    )

    overlap_summary <- venn_regions_n(gene_sets)

    all_genes <- Reduce(union, gene_sets)

    gene_membership <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

    for (i in seq_along(gene_sets)) {
        gene_membership[[dose_labels[i]]] <- gene_membership$gene %in% gene_sets[[i]]
    }

    save_csv(
        overlap_summary,
        file.path(figure_dir, paste0(filename_base, "_overlap_counts.csv")),
        row.names = FALSE
    )

    save_csv(
        gene_membership,
        file.path(figure_dir, paste0(filename_base, "_gene_membership.csv")),
        row.names = FALSE
    )

    message("Dose series Venn comparison complete.")

    return(
        list(
            gene_sets = gene_sets,
            overlap_summary = overlap_summary,
            gene_membership = gene_membership
        )
    )

}
