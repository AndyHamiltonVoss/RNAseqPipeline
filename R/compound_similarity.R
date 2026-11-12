################################################################################
# RNAseqPipeline
#
# File: compound_similarity.R
#
# Description:
#   Drug x drug similarity, comparing which genes each compound affects at
#   its highest tested dose. Two complementary views, both built from the
#   same per-compound top-dose DESeq2 (apeglm-shrunk) results:
#     - Jaccard overlap of significant (padj < alpha, |log2FC| > cutoff)
#       genes, direction-aware (a gene moving in opposite directions in two
#       compounds is not counted as shared).
#     - Pearson correlation of full-transcriptome log2FoldChange, capturing
#       sub-threshold directional agreement.
#   A Venn diagram isn't usable here: with ~20 compounds there's no tractable
#   N-way Venn, so this uses a clustered similarity-matrix heatmap instead,
#   matching the pheatmap pattern already established in
#   R/control_plate_correlation.R.
#
################################################################################

#===============================================================================
# Run Compound-vs-Compound Gene Similarity
#===============================================================================

run_compound_similarity <- function(config) {

    message("Running compound-vs-compound gene similarity...")

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    table_dir <- file.path(config$output_directory, "DESeq2", "Tables")
    figure_dir <- file.path(config$output_directory, "DESeq2", "Figures")

    shrunk_files <- list.files(
        table_dir,
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

    compounds <- sort(top_dose_table$compound)

    message(sprintf("Comparing %d compounds at their highest tested dose.", length(compounds)))

    per_compound <- lapply(compounds, function(cmpd) {

        group <- top_dose_table$group[top_dose_table$compound == cmpd]

        results <- read.csv(
            file.path(
                table_dir,
                paste0("DESeq2_group_", group, "_vs_", config$reference_group, "_shrunk.csv")
            )
        )

        classified <- classify_genes(
            results, "gene", "log2FoldChange", "padj",
            config$alpha, config$log2fc_cutoff
        )

        list(
            up = classified$up,
            down = classified$down,
            lfc = stats::setNames(results$log2FoldChange, results$gene)
        )

    })

    names(per_compound) <- compounds

    #===========================================================================
    # Jaccard Similarity (direction-aware) of Significant Genes
    #===========================================================================

    jaccard_matrix <- matrix(
        NA_real_,
        nrow = length(compounds), ncol = length(compounds),
        dimnames = list(compounds, compounds)
    )

    for (i in compounds) {
        for (j in compounds) {

            shared <- length(intersect(per_compound[[i]]$up, per_compound[[j]]$up)) +
                length(intersect(per_compound[[i]]$down, per_compound[[j]]$down))

            union_size <- length(union(
                union(per_compound[[i]]$up, per_compound[[i]]$down),
                union(per_compound[[j]]$up, per_compound[[j]]$down)
            ))

            jaccard_matrix[i, j] <- if (union_size == 0) NA_real_ else shared / union_size

        }
    }

    #===========================================================================
    # Jaccard Similarity, Split by Direction (Up-Only / Down-Only)
    #===========================================================================

    jaccard_up_matrix <- matrix(
        NA_real_,
        nrow = length(compounds), ncol = length(compounds),
        dimnames = list(compounds, compounds)
    )

    jaccard_down_matrix <- jaccard_up_matrix

    for (i in compounds) {
        for (j in compounds) {

            up_union <- length(union(per_compound[[i]]$up, per_compound[[j]]$up))

            jaccard_up_matrix[i, j] <- if (up_union == 0) {
                NA_real_
            } else {
                length(intersect(per_compound[[i]]$up, per_compound[[j]]$up)) / up_union
            }

            down_union <- length(union(per_compound[[i]]$down, per_compound[[j]]$down))

            jaccard_down_matrix[i, j] <- if (down_union == 0) {
                NA_real_
            } else {
                length(intersect(per_compound[[i]]$down, per_compound[[j]]$down)) / down_union
            }

        }
    }

    #===========================================================================
    # Pearson Correlation of Full-Transcriptome log2FoldChange
    #===========================================================================

    lfc_matrix <- Reduce(
        function(a, b) merge(a, b, by = "gene", all = TRUE),
        lapply(compounds, function(cmpd) {
            df <- data.frame(
                gene = names(per_compound[[cmpd]]$lfc),
                stringsAsFactors = FALSE
            )
            df[[cmpd]] <- as.numeric(per_compound[[cmpd]]$lfc)
            df
        })
    )

    rownames(lfc_matrix) <- lfc_matrix$gene
    lfc_matrix$gene <- NULL
    lfc_matrix <- lfc_matrix[, compounds, drop = FALSE]

    correlation_matrix <- stats::cor(
        as.matrix(lfc_matrix),
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    #===========================================================================
    # Save Tables
    #===========================================================================

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(jaccard_matrix, file.path(table_dir, "compound_jaccard_similarity.csv"))

    save_csv(jaccard_up_matrix, file.path(table_dir, "compound_jaccard_up.csv"))

    save_csv(jaccard_down_matrix, file.path(table_dir, "compound_jaccard_down.csv"))

    save_csv(correlation_matrix, file.path(table_dir, "compound_lfc_correlation.csv"))

    #===========================================================================
    # Save Heatmaps
    #===========================================================================

    heatmap_height <- max(config$figure_height, length(compounds) * 0.5)
    heatmap_width <- heatmap_height + 2

    if (isTRUE(config$save_png)) {

        pheatmap::pheatmap(
            jaccard_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: Jaccard Overlap of Significant Genes (Top Dose)",
            filename = file.path(figure_dir, "Compound_jaccard_similarity_heatmap.png"),
            width = heatmap_width,
            height = heatmap_height,
            res = config$figure_dpi
        )

        pheatmap::pheatmap(
            correlation_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: log2FoldChange Correlation (Top Dose)",
            filename = file.path(figure_dir, "Compound_LFC_correlation_heatmap.png"),
            width = heatmap_width,
            height = heatmap_height,
            res = config$figure_dpi
        )

        pheatmap::pheatmap(
            jaccard_up_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: Jaccard Overlap of Upregulated Genes (Top Dose)",
            filename = file.path(figure_dir, "Compound_jaccard_up_heatmap.png"),
            width = heatmap_width,
            height = heatmap_height,
            res = config$figure_dpi
        )

        pheatmap::pheatmap(
            jaccard_down_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: Jaccard Overlap of Downregulated Genes (Top Dose)",
            filename = file.path(figure_dir, "Compound_jaccard_down_heatmap.png"),
            width = heatmap_width,
            height = heatmap_height,
            res = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        pheatmap::pheatmap(
            jaccard_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: Jaccard Overlap of Significant Genes (Top Dose)",
            filename = file.path(figure_dir, "Compound_jaccard_similarity_heatmap.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

        pheatmap::pheatmap(
            correlation_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: log2FoldChange Correlation (Top Dose)",
            filename = file.path(figure_dir, "Compound_LFC_correlation_heatmap.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

        pheatmap::pheatmap(
            jaccard_up_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: Jaccard Overlap of Upregulated Genes (Top Dose)",
            filename = file.path(figure_dir, "Compound_jaccard_up_heatmap.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

        pheatmap::pheatmap(
            jaccard_down_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: Jaccard Overlap of Downregulated Genes (Top Dose)",
            filename = file.path(figure_dir, "Compound_jaccard_down_heatmap.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

    }

    #===========================================================================
    # Paired Comparison: Up-Overlap vs Down-Overlap, All Unique Compound Pairs
    #===========================================================================

    pair_idx <- utils::combn(compounds, 2)

    pair_comparison <- data.frame(
        compound_a = pair_idx[1, ],
        compound_b = pair_idx[2, ],
        jaccard_up = mapply(
            function(a, b) jaccard_up_matrix[a, b], pair_idx[1, ], pair_idx[2, ]
        ),
        jaccard_down = mapply(
            function(a, b) jaccard_down_matrix[a, b], pair_idx[1, ], pair_idx[2, ]
        ),
        stringsAsFactors = FALSE
    )

    pair_comparison <- pair_comparison[
        !is.na(pair_comparison$jaccard_up) & !is.na(pair_comparison$jaccard_down),
    ]

    wilcox_test <- stats::wilcox.test(
        pair_comparison$jaccard_up,
        pair_comparison$jaccard_down,
        paired = TRUE
    )

    pair_plot <- ggplot2::ggplot(
        pair_comparison,
        ggplot2::aes(x = jaccard_down, y = jaccard_up)
    ) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::coord_equal() +
        ggplot2::labs(
            title = "Upregulated vs. Downregulated Gene Overlap, All Compound Pairs",
            subtitle = sprintf(
                "Wilcoxon signed-rank p = %.3g; mean up = %.3f, mean down = %.3f (n = %d pairs)",
                wilcox_test$p.value, mean(pair_comparison$jaccard_up),
                mean(pair_comparison$jaccard_down), nrow(pair_comparison)
            ),
            x = "Jaccard overlap (downregulated genes)",
            y = "Jaccard overlap (upregulated genes)"
        ) +
        ggplot2::theme_bw()

    save_csv(
        pair_comparison,
        file.path(table_dir, "compound_jaccard_up_vs_down_pairs.csv"),
        row.names = FALSE
    )

    pair_plot_width <- max(config$figure_width, 11)
    pair_plot_height <- max(config$figure_height, 8)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "Compound_jaccard_up_vs_down_pairs.png"),
            plot = pair_plot,
            width = pair_plot_width,
            height = pair_plot_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "Compound_jaccard_up_vs_down_pairs.pdf"),
            plot = pair_plot,
            width = pair_plot_width,
            height = pair_plot_height
        )

    }

    message("Compound-vs-compound gene similarity complete.")

    return(
        list(
            top_dose_table = top_dose_table,
            jaccard_up_matrix = jaccard_up_matrix,
            jaccard_down_matrix = jaccard_down_matrix,
            pair_comparison = pair_comparison,
            wilcox_test = wilcox_test,
            jaccard_matrix = jaccard_matrix,
            correlation_matrix = correlation_matrix
        )
    )

}
