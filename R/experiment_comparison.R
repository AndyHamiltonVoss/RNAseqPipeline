################################################################################
# RNAseqPipeline
#
# File: experiment_comparison.R
#
# Description:
#   Gene-level reproducibility comparison between AMP0786 and AMP0961 for
#   the two drugs shared between them (MRL24_100nM, Rosi_1000nM), using
#   the library-adjusted, shared-Control-baseline DESeq2 contrasts already
#   written by run_deseq2() on the combined dataset (config_both.R,
#   results/compare_both/DESeq2/Tables/DESeq2_group_<drug>_<experiment>_vs_Control_
#   shrunk.csv). Mirrors R/compound_similarity.R's gene-level view:
#     - Pearson correlation of log2FoldChange (all tested genes).
#     - Direction-aware Jaccard overlap of significant (padj < alpha,
#       |log2FC| > cutoff) up- and down-regulated genes.
#     - Full pairwise Jaccard (up/down) across all 4 drug x experiment
#       entities (not just the matched same-drug pairs), clustered
#       pheatmaps -- same style as R/compound_similarity.R's
#       Compound_jaccard_up/down_heatmap.png.
#
################################################################################

#===============================================================================
# Run AMP0786 vs. AMP0961 Gene-Level Reproducibility Comparison
#===============================================================================

run_experiment_comparison <- function(alpha = 0.05, log2fc_cutoff = 1) {

    message("Running AMP0786 vs. AMP0961 gene-level reproducibility comparison...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    table_dir <- "results/compare_both/DESeq2/Tables"
    out_table_dir <- "results/compare_both/Comparison/Tables"
    out_figure_dir <- "results/compare_both/Comparison/Figures"

    dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)

    drugs <- c("MRL24_100nM", "Rosi_1000nM")

    file_for <- function(drug, experiment) {
        file.path(table_dir, paste0("DESeq2_group_", drug, "_", experiment, "_vs_Control_shrunk.csv"))
    }

    summary_rows <- list()
    per_entity <- list()

    for (drug in drugs) {

        amp0786 <- read.csv(file_for(drug, "AMP0786"))
        amp0961 <- read.csv(file_for(drug, "AMP0961"))

        merged <- merge(
            amp0786[, c("gene", "log2FoldChange", "padj")],
            amp0961[, c("gene", "log2FoldChange", "padj")],
            by = "gene", suffixes = c("_AMP0786", "_AMP0961")
        )

        correlation <- stats::cor(
            merged$log2FoldChange_AMP0786, merged$log2FoldChange_AMP0961,
            method = "pearson", use = "pairwise.complete.obs"
        )

        up_786 <- merged$gene[!is.na(merged$padj_AMP0786) & merged$padj_AMP0786 < alpha & merged$log2FoldChange_AMP0786 > log2fc_cutoff]
        down_786 <- merged$gene[!is.na(merged$padj_AMP0786) & merged$padj_AMP0786 < alpha & merged$log2FoldChange_AMP0786 < -log2fc_cutoff]
        up_961 <- merged$gene[!is.na(merged$padj_AMP0961) & merged$padj_AMP0961 < alpha & merged$log2FoldChange_AMP0961 > log2fc_cutoff]
        down_961 <- merged$gene[!is.na(merged$padj_AMP0961) & merged$padj_AMP0961 < alpha & merged$log2FoldChange_AMP0961 < -log2fc_cutoff]

        jaccard_up <- length(intersect(up_786, up_961)) / length(union(up_786, up_961))
        jaccard_down <- length(intersect(down_786, down_961)) / length(union(down_786, down_961))

        per_entity[[paste0(drug, "_AMP0786")]] <- list(up = up_786, down = down_786)
        per_entity[[paste0(drug, "_AMP0961")]] <- list(up = up_961, down = down_961)

        summary_rows[[drug]] <- data.frame(
            drug = drug,
            n_genes_compared = nrow(merged),
            pearson_correlation = correlation,
            n_up_AMP0786 = length(up_786),
            n_down_AMP0786 = length(down_786),
            n_up_AMP0961 = length(up_961),
            n_down_AMP0961 = length(down_961),
            n_up_shared = length(intersect(up_786, up_961)),
            n_down_shared = length(intersect(down_786, down_961)),
            jaccard_up = jaccard_up,
            jaccard_down = jaccard_down
        )

        save_csv(merged, file.path(out_table_dir, paste0(drug, "_merged_lfc.csv")), row.names = FALSE)

        #=======================================================================
        # Scatter Plot: log2FoldChange AMP0786 vs. AMP0961, Colored by Concordance
        #=======================================================================

        merged$concordance <- "Neither significant"
        merged$concordance[merged$gene %in% intersect(up_786, up_961)] <- "Concordant up"
        merged$concordance[merged$gene %in% intersect(down_786, down_961)] <- "Concordant down"
        merged$concordance[merged$gene %in% intersect(up_786, down_961) | merged$gene %in% intersect(down_786, up_961)] <- "Discordant"
        one_sided <- (merged$gene %in% c(up_786, down_786)) != (merged$gene %in% c(up_961, down_961))
        merged$concordance[one_sided & merged$concordance == "Neither significant"] <- "One-sided"

        merged$concordance <- factor(
            merged$concordance,
            levels = c("Concordant up", "Concordant down", "Discordant", "One-sided", "Neither significant")
        )

        plot <- ggplot2::ggplot(merged, ggplot2::aes(x = log2FoldChange_AMP0786, y = log2FoldChange_AMP0961, color = concordance)) +
            ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
            ggplot2::geom_point(alpha = 0.5, size = 0.8) +
            ggplot2::scale_color_manual(values = c(
                "Concordant up" = "firebrick", "Concordant down" = "steelblue",
                "Discordant" = "purple", "One-sided" = "darkorange", "Neither significant" = "grey80"
            )) +
            ggplot2::labs(
                title = sprintf("%s: AMP0786 vs. AMP0961 log2FoldChange", drug),
                subtitle = sprintf(
                    "n=%d genes, Pearson r=%.3f; Jaccard up=%.3f, down=%.3f",
                    nrow(merged), correlation, jaccard_up, jaccard_down
                ),
                x = "log2FoldChange (AMP0786 vs. Control)", y = "log2FoldChange (AMP0961 vs. Control)",
                color = "Significance"
            ) +
            ggplot2::theme_bw()

        ggplot2::ggsave(file.path(out_figure_dir, paste0(drug, "_lfc_scatter.png")), plot = plot, width = 9, height = 7, dpi = 300)
        ggplot2::ggsave(file.path(out_figure_dir, paste0(drug, "_lfc_scatter.pdf")), plot = plot, width = 9, height = 7)

    }

    summary_table <- do.call(rbind, summary_rows)
    rownames(summary_table) <- NULL

    save_csv(summary_table, file.path(out_table_dir, "experiment_comparison_summary.csv"), row.names = FALSE)

    #===========================================================================
    # Full Pairwise Jaccard Heatmaps (All 4 Drug x Experiment Entities)
    #===========================================================================

    entities <- names(per_entity)

    jaccard_up_matrix <- matrix(
        NA_real_, nrow = length(entities), ncol = length(entities),
        dimnames = list(entities, entities)
    )
    jaccard_down_matrix <- jaccard_up_matrix

    for (i in entities) {
        for (j in entities) {

            up_union <- length(union(per_entity[[i]]$up, per_entity[[j]]$up))
            jaccard_up_matrix[i, j] <- if (up_union == 0) NA_real_ else length(intersect(per_entity[[i]]$up, per_entity[[j]]$up)) / up_union

            down_union <- length(union(per_entity[[i]]$down, per_entity[[j]]$down))
            jaccard_down_matrix[i, j] <- if (down_union == 0) NA_real_ else length(intersect(per_entity[[i]]$down, per_entity[[j]]$down)) / down_union

        }
    }

    save_csv(jaccard_up_matrix, file.path(out_table_dir, "experiment_jaccard_up.csv"))
    save_csv(jaccard_down_matrix, file.path(out_table_dir, "experiment_jaccard_down.csv"))

    heatmap_height <- length(entities) * 1.2 + 2
    heatmap_width <- heatmap_height + 2

    pheatmap::pheatmap(
        jaccard_up_matrix,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        display_numbers = TRUE,
        number_format = "%.2f",
        main = "AMP0786 vs. AMP0961: Jaccard Overlap of Upregulated Genes",
        filename = file.path(out_figure_dir, "Experiment_jaccard_up_heatmap.png"),
        width = heatmap_width,
        height = heatmap_height,
        res = 300
    )

    pheatmap::pheatmap(
        jaccard_down_matrix,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        display_numbers = TRUE,
        number_format = "%.2f",
        main = "AMP0786 vs. AMP0961: Jaccard Overlap of Downregulated Genes",
        filename = file.path(out_figure_dir, "Experiment_jaccard_down_heatmap.png"),
        width = heatmap_width,
        height = heatmap_height,
        res = 300
    )

    message("AMP0786 vs. AMP0961 gene-level comparison complete.")

    print(summary_table)

    return(
        list(
            summary = summary_table,
            jaccard_up_matrix = jaccard_up_matrix,
            jaccard_down_matrix = jaccard_down_matrix
        )
    )

}
