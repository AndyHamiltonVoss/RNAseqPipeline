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

    }

    message("Compound-vs-compound gene similarity complete.")

    return(
        list(
            top_dose_table = top_dose_table,
            jaccard_matrix = jaccard_matrix,
            correlation_matrix = correlation_matrix
        )
    )

}
