################################################################################
# RNAseqPipeline
#
# File: umap.R
#
# Description:
#   Sample-level UMAP: nonlinear 2D embedding on variance-stabilized counts,
#   same top-variable-gene input as the PCA in qc.R.
#
################################################################################

#===============================================================================
# Run UMAP
#===============================================================================

run_umap <- function(counts, metadata, config) {

    message("Running UMAP...")

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.")
    }

    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
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

    n_samples <- ncol(vst_data)

    n_neighbors <- min(config$umap_n_neighbors, n_samples - 1)

    set.seed(config$random_seed)

    embedding <- uwot::umap(
        t(vst_data[top_genes, , drop = FALSE]),
        n_neighbors = n_neighbors,
        min_dist = config$umap_min_dist,
        n_components = 2
    )

    scores <- data.frame(
        sample = colnames(vst_data),
        UMAP1 = embedding[, 1],
        UMAP2 = embedding[, 2]
    )

    scores <- merge(
        scores,
        metadata,
        by = "sample"
    )

    output_dir <- file.path(config$output_directory, "QC", "UMAP")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        scores,
        file.path(table_dir, "umap_coordinates.csv"),
        row.names = FALSE
    )

    n_groups <- length(unique(scores$group))

    umap_plot <- ggplot2::ggplot(
        scores,
        ggplot2::aes(x = UMAP1, y = UMAP2, color = group)
    ) +
        ggplot2::geom_point(size = 2.5, alpha = 0.85) +
        ggplot2::labs(
            title = "UMAP of Samples",
            x = "UMAP1",
            y = "UMAP2",
            color = "Group"
        ) +
        ggplot2::theme_bw()

    # A legend with too many entries pushes the plot panel itself off the
    # canvas, so it is dropped above a group-count threshold rather than
    # rendering an unreadable (or invisible) plot.
    if (n_groups > config$umap_max_legend_groups) {

        message(
            sprintf(
                "%d groups exceeds umap_max_legend_groups (%d); legend suppressed.",
                n_groups,
                config$umap_max_legend_groups
            )
        )

        umap_plot <- umap_plot + ggplot2::theme(legend.position = "none")

    }

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "UMAP_plot.png"),
            plot = umap_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "UMAP_plot.pdf"),
            plot = umap_plot,
            width = config$figure_width,
            height = config$figure_height
        )

    }

    message("UMAP complete.")

    return(
        list(
            scores = scores,
            plot = umap_plot
        )
    )

}
