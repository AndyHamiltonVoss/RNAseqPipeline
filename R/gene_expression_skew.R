################################################################################
# RNAseqPipeline
#
# File: gene_expression_skew.R
#
# Description:
#   Diagnostic: does gene-level log2FoldChange (apeglm-shrunk DESeq2 results
#   already written by run_deseq2(), vs. the reference group) skew toward
#   down-regulation at each compound's top tested dose, the way GSVA KEGG
#   pathway activity does (R/pathway_downregulation_skew.R)? Same two views,
#   at the gene level instead of the pathway level:
#     - Per-compound % of all tested genes with negative vs. positive
#       log2FoldChange (not just the ones crossing the significance
#       cutoff), grouped bars.
#     - Pooled distribution of log2FoldChange across every compound x gene
#       combination.
#   This is the control for R/pathway_downregulation_skew.R: if genes show
#   little or no analogous skew, the shared negative shift found at the
#   pathway level is introduced by GSVA's aggregation/scoring, not inherited
#   from the underlying differential expression calls.
#
################################################################################

#===============================================================================
# Run Gene Expression Skew Diagnostic
#===============================================================================

run_gene_expression_skew <- function(config) {

    message("Running gene-level expression skew diagnostic...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    deseq2_table_dir <- file.path(config$output_directory, "DESeq2", "Tables")
    diagnostic_dir <- file.path(config$output_directory, "gene_diagnostic")

    #===========================================================================
    # Discover Each Compound's Top Dose
    #===========================================================================

    shrunk_files <- list.files(
        deseq2_table_dir,
        pattern = "^DESeq2_group_.*_shrunk\\.csv$"
    )

    group_names <- sub(
        paste0("_vs_", config$reference_group, "_shrunk\\.csv$"), "",
        sub("^DESeq2_group_", "", shrunk_files)
    )

    parsed <- lapply(group_names, parse_dose_group)

    keep <- !vapply(parsed, is.null, logical(1))

    dose_table <- data.frame(
        group = group_names[keep],
        compound = vapply(parsed[keep], `[[`, character(1), "compound"),
        dose_nM = vapply(parsed[keep], `[[`, numeric(1), "dose_nM")
    )

    top_dose_table <- do.call(
        rbind,
        lapply(split(dose_table, dose_table$compound), function(d) {
            d[which.max(d$dose_nM), ]
        })
    )

    shrunk_file_for <- function(group) {
        file.path(deseq2_table_dir, paste0("DESeq2_group_", group, "_vs_", config$reference_group, "_shrunk.csv"))
    }

    top_dose_table$has_deseq2 <- vapply(top_dose_table$group, function(g) file.exists(shrunk_file_for(g)), logical(1))
    eligible <- top_dose_table[top_dose_table$has_deseq2, ]
    compounds <- sort(eligible$compound)

    message(sprintf("Checking gene log2FoldChange skew across %d compounds at their highest tested dose.", length(compounds)))

    #===========================================================================
    # Read log2FoldChange for All Genes (Not Just Significant Ones)
    #===========================================================================

    long_data <- do.call(rbind, lapply(compounds, function(cmpd) {
        group <- eligible$group[eligible$compound == cmpd]
        res <- read.csv(shrunk_file_for(group))
        data.frame(compound = cmpd, gene = res$gene, log2FoldChange = res$log2FoldChange)
    }))

    long_data <- long_data[!is.na(long_data$log2FoldChange), ]

    #===========================================================================
    # Per-Compound % of Genes With Negative vs. Positive log2FoldChange
    #===========================================================================

    pct_negative_df <- do.call(rbind, lapply(compounds, function(cmpd) {
        sub <- long_data[long_data$compound == cmpd, ]
        data.frame(
            compound = cmpd,
            pct_negative = mean(sub$log2FoldChange < 0) * 100,
            pct_positive = mean(sub$log2FoldChange > 0) * 100
        )
    }))

    overall_mean_pct_negative <- mean(pct_negative_df$pct_negative)
    overall_mean_pct_positive <- mean(pct_negative_df$pct_positive)

    pct_negative_df <- pct_negative_df[order(-pct_negative_df$pct_negative), ]

    pct_long_df <- rbind(
        data.frame(compound = pct_negative_df$compound, direction = "Down", pct = pct_negative_df$pct_negative),
        data.frame(compound = pct_negative_df$compound, direction = "Up", pct = pct_negative_df$pct_positive)
    )

    #===========================================================================
    # Pooled log2FoldChange Distribution
    #===========================================================================

    pooled_median <- stats::median(long_data$log2FoldChange)
    pooled_pct_negative <- mean(long_data$log2FoldChange < 0) * 100

    #===========================================================================
    # Save Tables
    #===========================================================================

    dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        pct_negative_df,
        file.path(diagnostic_dir, "gene_pct_negative_by_compound.csv"),
        row.names = FALSE
    )

    pooled_summary <- data.frame(
        n_compounds = length(compounds),
        n_genes = length(unique(long_data$gene)),
        pooled_pct_negative = pooled_pct_negative,
        pooled_pct_positive = mean(long_data$log2FoldChange > 0) * 100,
        pooled_median_log2FoldChange = pooled_median,
        overall_mean_pct_negative_across_compounds = overall_mean_pct_negative,
        overall_mean_pct_positive_across_compounds = overall_mean_pct_positive
    )

    save_csv(
        pooled_summary,
        file.path(diagnostic_dir, "gene_expression_skew_summary.csv"),
        row.names = FALSE
    )

    #===========================================================================
    # Plot: Per-Compound % Down vs. % Up Grouped Bar Chart
    #===========================================================================

    pct_long_df$compound <- factor(pct_long_df$compound, levels = pct_negative_df$compound)
    pct_long_df$direction <- factor(pct_long_df$direction, levels = c("Down", "Up"))

    panel1 <- ggplot2::ggplot(pct_long_df, ggplot2::aes(x = compound, y = pct, fill = direction)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
        ggplot2::geom_hline(yintercept = 50, linetype = "dashed", color = "grey30") +
        ggplot2::scale_fill_manual(values = c(Down = "steelblue", Up = "darkorange")) +
        ggplot2::labs(
            title = "Gene-level log2FoldChange at top dose, by compound",
            subtitle = sprintf(
                "%% of all tested genes with negative vs. positive log2FoldChange (not just significant ones); mean down = %.1f%%, mean up = %.1f%%\ndashed line = 50%% (no skew)",
                overall_mean_pct_negative, overall_mean_pct_positive
            ),
            x = "Compound", y = "% of genes", fill = "Direction"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    #===========================================================================
    # Plot: Pooled log2FoldChange Histogram
    #===========================================================================

    panel2 <- ggplot2::ggplot(long_data, ggplot2::aes(x = log2FoldChange)) +
        ggplot2::geom_histogram(bins = 80, fill = "steelblue", color = NA) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
        ggplot2::geom_vline(xintercept = pooled_median, linetype = "dotted", color = "firebrick") +
        ggplot2::coord_cartesian(xlim = c(-2, 2)) +
        ggplot2::labs(
            title = "Pooled distribution of log2FoldChange across all compounds and genes",
            subtitle = sprintf(
                "%d compounds x %d genes; %.1f%% negative; dashed = 0, dotted = pooled median (%.4f); x-axis clipped to [-2, 2]",
                length(compounds), length(unique(long_data$gene)), pooled_pct_negative, pooled_median
            ),
            x = "log2FoldChange (vs. reference group)", y = "Count"
        ) +
        ggplot2::theme_bw()

    #===========================================================================
    # Combine and Save
    #===========================================================================

    if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- patchwork::wrap_plots(panel1, panel2, ncol = 1, heights = c(1.1, 1))
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
        combined <- cowplot::plot_grid(panel1, panel2, ncol = 1, rel_heights = c(1.1, 1))
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
        combined <- gridExtra::arrangeGrob(panel1, panel2, ncol = 1, heights = c(1.1, 1))
    } else {
        stop("One of 'patchwork', 'cowplot', or 'gridExtra' is required to combine the two panels.")
    }

    plot_width <- max(config$figure_width, 11)
    plot_height <- max(config$figure_height, 11)

    if (isTRUE(config$save_png)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Gene_expression_skew.png"),
            plot = combined, width = plot_width, height = plot_height, dpi = config$figure_dpi
        )
    }

    if (isTRUE(config$save_pdf)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Gene_expression_skew.pdf"),
            plot = combined, width = plot_width, height = plot_height
        )
    }

    message("Gene-level expression skew diagnostic complete.")

    return(
        list(
            pct_negative_by_compound = pct_negative_df,
            pct_long_by_compound = pct_long_df,
            pooled_summary = pooled_summary
        )
    )

}
