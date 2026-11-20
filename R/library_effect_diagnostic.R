################################################################################
# RNAseqPipeline
#
# File: library_effect_diagnostic.R
#
# Description:
#   Quantifies how much of each compound's apparent vs-Control signal is
#   actually just its plate/library identity. In AMP0961, each compound's
#   entire dose series sits on a fixed pair of the 10 libraries (plates)
#   -- e.g. Dar only ever appears on libraries AV1/AV2 -- so any
#   plate-to-plate technical variation is inseparable from that compound's
#   "real" biology using a single compound-vs-Control contrast alone.
#   R/control_plate_correlation.R already shows Control-well profiles
#   differ modestly but measurably by plate (0.95-1.00 correlation, with
#   real dendrogram structure); this module goes further and asks how
#   much that plate structure specifically overlaps with each compound's
#   OWN apparent signal.
#
#   `Control` is the only group spread across all 10 libraries, so a
#   Control-only, library-vs-library DESeq2 contrast isolates a "pure
#   plate effect" genome-wide LFC profile with zero drug involved. That
#   profile is then compared (Pearson correlation of LFC vectors, plus a
#   secondary Jaccard overlap of significant up/down genes) against each
#   compound's own real vs-Control LFC profile already on disk from
#   run_deseq2(). A high correlation means a large share of that
#   compound's "signal" is actually just its plate pair.
#
#   Bypasses run_deseq2() deliberately: that function is a single
#   hardcoded full-dataset pass (one fitted model, contrasts selected via
#   a "^group_.*_vs_" name regex) not meant to be re-entered for an ad
#   hoc Control-only contrast -- reusing it here would mean mutating
#   config fields other code in the same run depends on, and colliding
#   with its own output files. This module calls
#   DESeq2::DESeqDataSetFromMatrix()/DESeq()/results() directly instead,
#   on a small Control-only subset.
#
################################################################################

#===============================================================================
# Discover Each Compound's Primary Library Pair
#===============================================================================
#
# For each non-Control, non-T007-arm group (parse_dose_group() reused
# from R/dose_response.R), tally which libraries its samples fall on.
# The 2 most-sample-heavy libraries become that compound's "primary
# pair"; is_multi_library flags compounds where a third library also
# contributes samples (e.g. one dose shifted onto an adjacent plate
# pair), since the primary-pair simplification is an approximation for
# those. Discovered dynamically from the metadata, not hardcoded.
#===============================================================================

discover_compound_library_pairs <- function(metadata, exclude_pattern = "T007") {

    groups <- unique(as.character(metadata$group))
    groups <- groups[groups != "Control"]

    if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {
        groups <- groups[!grepl(exclude_pattern, groups)]
    }

    parsed <- lapply(groups, parse_dose_group)
    keep <- !vapply(parsed, is.null, logical(1))

    dose_table <- data.frame(
        group = groups[keep],
        compound = vapply(parsed[keep], `[[`, character(1), "compound"),
        stringsAsFactors = FALSE
    )

    compounds <- sort(unique(dose_table$compound))

    rows <- lapply(compounds, function(cmpd) {

        cmpd_groups <- dose_table$group[dose_table$compound == cmpd]
        cmpd_samples <- metadata[metadata$group %in% cmpd_groups, , drop = FALSE]

        # table() on a factor column retains unused levels (e.g. all 10
        # libraries, most with a zero count, after subsetting to one
        # compound's rows) -- coerce to character first so n_libs below
        # only counts libraries actually present for this compound.
        lib_counts <- sort(table(as.character(cmpd_samples$library)), decreasing = TRUE)

        n_libs <- length(lib_counts)
        primary <- names(lib_counts)[seq_len(min(2, n_libs))]

        data.frame(
            compound = cmpd,
            library_a = sort(primary)[1],
            library_b = sort(primary)[2],
            is_multi_library = n_libs > 2,
            n_libraries_total = n_libs,
            stringsAsFactors = FALSE
        )

    })

    do.call(rbind, rows)

}

#===============================================================================
# Compute the Pure Plate-Effect LFC Profile for One Library Pair
# (Control Samples Only, No Drug Involved)
#===============================================================================

compute_library_pair_effect <- function(counts, metadata, library_a, library_b, config) {

    control_metadata <- metadata[
        metadata$group == config$reference_group & metadata$library %in% c(library_a, library_b),
        ,
        drop = FALSE
    ]

    control_metadata$library <- factor(control_metadata$library, levels = c(library_a, library_b))

    control_counts <- counts[, rownames(control_metadata), drop = FALSE]

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(as.matrix(control_counts)),
        colData = control_metadata,
        design = ~library
    )

    dds <- DESeq2::DESeq(dds)

    res <- DESeq2::results(dds, contrast = c("library", library_b, library_a), alpha = config$alpha)

    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

    res_df

}

#===============================================================================
# Run the Library-Effect Diagnostic (All Compounds)
#===============================================================================

run_library_effect_diagnostic <- function(counts, metadata, config) {

    message("Running library/plate-effect diagnostic...")

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    pairs_table <- discover_compound_library_pairs(metadata, exclude_pattern = config$exclude_groups_pattern)

    table_dir <- file.path(config$output_directory, "Library_Effect", "Tables")
    figure_dir <- file.path(config$output_directory, "Library_Effect", "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    deseq2_table_dir <- file.path(config$output_directory, "DESeq2", "Tables")

    # Discover each compound's top tested dose, same convention used
    # throughout (parse_dose_group() + which.max(dose_nM)).
    all_groups <- unique(as.character(metadata$group))
    all_groups <- all_groups[all_groups != "Control"]
    if (!is.null(config$exclude_groups_pattern) && nzchar(config$exclude_groups_pattern)) {
        all_groups <- all_groups[!grepl(config$exclude_groups_pattern, all_groups)]
    }
    parsed_all <- lapply(all_groups, parse_dose_group)
    keep_all <- !vapply(parsed_all, is.null, logical(1))
    dose_table <- data.frame(
        group = all_groups[keep_all],
        compound = vapply(parsed_all[keep_all], `[[`, character(1), "compound"),
        dose_nM = vapply(parsed_all[keep_all], `[[`, numeric(1), "dose_nM"),
        stringsAsFactors = FALSE
    )
    top_dose_table <- do.call(
        rbind,
        lapply(split(dose_table, dose_table$compound), function(d) d[which.max(d$dose_nM), ])
    )

    # Compute each unique library pair's pure plate effect once (several
    # compounds share a pair, e.g. T007 arms, SR1664/SR16832/T007-alone).
    pair_cache <- list()

    pair_key <- function(a, b) paste0(a, "_vs_", b)

    summary_rows <- list()

    for (i in seq_len(nrow(pairs_table))) {

        cmpd <- pairs_table$compound[i]
        library_a <- pairs_table$library_a[i]
        library_b <- pairs_table$library_b[i]
        key <- pair_key(library_a, library_b)

        if (is.null(pair_cache[[key]])) {

            message(sprintf("Computing pure plate effect for %s...", key))

            pair_effect <- compute_library_pair_effect(counts, metadata, library_a, library_b, config)
            pair_cache[[key]] <- pair_effect

            save_csv(
                pair_effect,
                file.path(table_dir, paste0("library_pair_effect_", key, ".csv")),
                row.names = FALSE
            )

        }

        library_effect <- pair_cache[[key]]

        compound_group <- top_dose_table$group[top_dose_table$compound == cmpd]
        compound_file <- file.path(
            deseq2_table_dir,
            paste0("DESeq2_group_", compound_group, "_vs_", config$reference_group, "_raw.csv")
        )

        if (!file.exists(compound_file)) {
            warning("No DESeq2 results found for compound: ", cmpd, "; skipping library-effect comparison.")
            next
        }

        compound_results <- utils::read.csv(compound_file)

        merged <- merge(
            compound_results[, c("gene", "log2FoldChange", "padj")],
            library_effect[, c("gene", "log2FoldChange", "padj")],
            by = "gene", suffixes = c("_compound", "_library")
        )

        pearson_correlation <- stats::cor(
            merged$log2FoldChange_compound, merged$log2FoldChange_library,
            method = "pearson", use = "pairwise.complete.obs"
        )

        compound_sig <- classify_genes(merged, "gene", "log2FoldChange_compound", "padj_compound", config$alpha, config$log2fc_cutoff)
        library_sig <- classify_genes(merged, "gene", "log2FoldChange_library", "padj_library", config$alpha, config$log2fc_cutoff)

        jaccard <- function(a, b) {
            union_n <- length(union(a, b))
            if (union_n == 0) NA_real_ else length(intersect(a, b)) / union_n
        }

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
            compound = cmpd,
            top_dose_group = compound_group,
            library_pair = key,
            is_multi_library = pairs_table$is_multi_library[i],
            n_genes_compared = nrow(merged),
            n_significant_library_genes = length(library_sig$up) + length(library_sig$down),
            pearson_correlation = pearson_correlation,
            jaccard_up = jaccard(compound_sig$up, library_sig$up),
            jaccard_down = jaccard(compound_sig$down, library_sig$down),
            stringsAsFactors = FALSE
        )

    }

    summary_table <- do.call(rbind, summary_rows)
    summary_table <- summary_table[order(-abs(summary_table$pearson_correlation)), ]

    save_csv(
        summary_table,
        file.path(table_dir, "compound_library_effect_summary.csv"),
        row.names = FALSE
    )

    summary_table$compound <- factor(summary_table$compound, levels = rev(summary_table$compound))

    plot_height <- max(config$figure_height, nrow(summary_table) * 0.3)

    bar_plot <- ggplot2::ggplot(
        summary_table,
        ggplot2::aes(x = pearson_correlation, y = compound, fill = is_multi_library)
    ) +
        ggplot2::geom_col() +
        ggplot2::geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
        ggplot2::scale_fill_manual(
            values = c(`FALSE` = "steelblue", `TRUE` = "firebrick"),
            labels = c(`FALSE` = "Clean 2-library pair", `TRUE` = "Spans >2 libraries (approximate)")
        ) +
        ggplot2::labs(
            title = "Compound Signal vs. Pure Plate Effect (Top Dose)",
            x = "Pearson correlation: compound LFC vs. library-pair LFC (Control-only)",
            y = NULL,
            fill = NULL
        ) +
        ggplot2::theme_bw()

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "Compound_Library_Effect_Correlation.png"),
            plot = bar_plot, width = config$figure_width, height = plot_height, dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "Compound_Library_Effect_Correlation.pdf"),
            plot = bar_plot, width = config$figure_width, height = plot_height
        )

    }

    message("Library/plate-effect diagnostic complete.")

    list(summary_table = summary_table, pair_effects = pair_cache, pairs_table = pairs_table)

}

#===============================================================================
# Control-Only PCA, Colored by Library
#===============================================================================
#
# Same subset -> VST -> embed pattern as run_umap_top_dose() in R/umap.R,
# but PCA (matching R/qc.R's existing PCA style) instead of UMAP, and
# restricted to Control samples colored by library instead of group --
# the most direct visual check for whether plate identity creates
# structure in otherwise-untreated samples.
#===============================================================================

run_control_library_pca <- function(counts, metadata, config) {

    message("Running Control-only PCA by library...")

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    control_metadata <- metadata[metadata$group == config$reference_group, , drop = FALSE]
    control_counts <- counts[, rownames(control_metadata), drop = FALSE]

    vst_data <- tryCatch(
        DESeq2::vst(round(as.matrix(control_counts)), blind = TRUE),
        error = function(e) {
            warning("DESeq2::vst() failed (", conditionMessage(e), "); falling back to log2(counts + 1).")
            log2(as.matrix(control_counts) + 1)
        }
    )

    gene_vars <- apply(vst_data, 1, stats::var)
    n_top <- min(config$pca_top_variable_genes, nrow(vst_data))
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_top)]

    pca <- stats::prcomp(
        t(vst_data[top_genes, , drop = FALSE]),
        center = isTRUE(config$pca_center),
        scale. = isTRUE(config$pca_scale)
    )

    variance_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100

    scores <- as.data.frame(pca$x)
    scores$sample <- rownames(scores)
    scores <- merge(scores, control_metadata, by = "sample")

    output_dir <- file.path(config$output_directory, "Library_Effect")
    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(scores, file.path(table_dir, "control_pca_scores.csv"), row.names = FALSE)

    pca_plot <- ggplot2::ggplot(
        scores,
        ggplot2::aes(x = PC1, y = PC2, color = library)
    ) +
        ggplot2::geom_point(size = 3, alpha = 0.85) +
        ggplot2::labs(
            title = "Control Samples Only: PCA by Library",
            x = sprintf("PC1 (%.1f%% variance)", variance_explained[1]),
            y = sprintf("PC2 (%.1f%% variance)", variance_explained[2]),
            color = "Library"
        ) +
        ggplot2::theme_bw()

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "Control_PCA_by_Library.png"),
            plot = pca_plot, width = config$figure_width, height = config$figure_height, dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "Control_PCA_by_Library.pdf"),
            plot = pca_plot, width = config$figure_width, height = config$figure_height
        )

    }

    message("Control-only PCA by library complete.")

    list(scores = scores, plot = pca_plot)

}
