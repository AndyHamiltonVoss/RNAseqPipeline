################################################################################
# RNAseqPipeline
#
# File: qc.R
#
# Description:
#   Sample-level quality control: PCA on variance-stabilized counts.
#
################################################################################

#===============================================================================
# Run PCA-Based Quality Control
#===============================================================================

run_qc <- function(counts, metadata, config) {

    message("Running QC / PCA...")

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    counts <- counts[, rownames(metadata), drop = FALSE]

    vst_data <- tryCatch(
        DESeq2::vst(round(as.matrix(counts)), blind = TRUE),
        error = function(e) {
            warning(
                "DESeq2::vst() failed (", conditionMessage(e), "); ",
                "falling back to log2(counts + 1)."
            )
            log2(as.matrix(counts) + 1)
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

    scores <- merge(
        scores,
        metadata,
        by = "sample"
    )

    output_dir <- file.path(config$output_directory, "QC")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        scores,
        file.path(table_dir, "pca_scores.csv"),
        row.names = FALSE
    )

    save_csv(
        data.frame(
            PC = paste0("PC", seq_along(variance_explained)),
            percent_variance = variance_explained
        ),
        file.path(table_dir, "pca_variance_explained.csv"),
        row.names = FALSE
    )

    pca_plot <- ggplot2::ggplot(
        scores,
        ggplot2::aes(x = PC1, y = PC2, color = group)
    ) +
        ggplot2::geom_point(size = 2.5, alpha = 0.85) +
        ggplot2::stat_ellipse(type = "norm", level = 0.68, linewidth = 0.4) +
        ggplot2::labs(
            title = "PCA of Samples",
            x = sprintf("PC1 (%.1f%% variance)", variance_explained[1]),
            y = sprintf("PC2 (%.1f%% variance)", variance_explained[2]),
            color = "Group"
        ) +
        ggplot2::theme_bw()

    scree_data <- data.frame(
        PC = factor(seq_along(variance_explained), levels = seq_along(variance_explained)),
        percent_variance = variance_explained
    )

    scree_plot <- ggplot2::ggplot(
        utils::head(scree_data, 10),
        ggplot2::aes(x = PC, y = percent_variance)
    ) +
        ggplot2::geom_col(fill = "#1b9e77") +
        ggplot2::labs(
            title = "PCA Scree Plot",
            x = "Principal Component",
            y = "Percent Variance Explained"
        ) +
        ggplot2::theme_bw()

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "PCA_plot.png"),
            plot = pca_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

        ggplot2::ggsave(
            file.path(figure_dir, "PCA_scree_plot.png"),
            plot = scree_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "PCA_plot.pdf"),
            plot = pca_plot,
            width = config$figure_width,
            height = config$figure_height
        )

        ggplot2::ggsave(
            file.path(figure_dir, "PCA_scree_plot.pdf"),
            plot = scree_plot,
            width = config$figure_width,
            height = config$figure_height
        )

    }

    #===============================================================================
    # Heatmap of Top Variable Genes
    #===============================================================================

    n_heatmap <- min(config$heatmap_top_genes, nrow(vst_data))

    heatmap_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_heatmap)]

    heatmap_matrix <- vst_data[heatmap_genes, , drop = FALSE]

    save_csv(
        data.frame(
            gene = heatmap_genes,
            variance = gene_vars[heatmap_genes]
        ),
        file.path(table_dir, "heatmap_top_genes.csv"),
        row.names = FALSE
    )

    annotation_col <- data.frame(
        Group = metadata[colnames(heatmap_matrix), "group"],
        row.names = colnames(heatmap_matrix)
    )

    heatmap_height <- max(config$figure_height, 3 + n_heatmap * 0.15)
    heatmap_width <- max(config$figure_width, ncol(heatmap_matrix) * 0.1) + 2

    if (isTRUE(config$save_png)) {

        pheatmap::pheatmap(
            heatmap_matrix,
            scale = "row",
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            annotation_col = annotation_col,
            show_colnames = FALSE,
            show_rownames = n_heatmap <= 50,
            main = "Top Variable Genes",
            filename = file.path(figure_dir, "Heatmap_top_genes.png"),
            width = heatmap_width,
            height = heatmap_height,
            res = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        pheatmap::pheatmap(
            heatmap_matrix,
            scale = "row",
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            annotation_col = annotation_col,
            show_colnames = FALSE,
            show_rownames = n_heatmap <= 50,
            main = "Top Variable Genes",
            filename = file.path(figure_dir, "Heatmap_top_genes.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

    }

    message("QC / PCA complete.")

    return(
        list(
            pca = pca,
            scores = scores,
            variance_explained = variance_explained,
            plot = pca_plot,
            scree_plot = scree_plot,
            heatmap_matrix = heatmap_matrix
        )
    )

}
