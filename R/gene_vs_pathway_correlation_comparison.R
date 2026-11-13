################################################################################
# RNAseqPipeline
#
# File: gene_vs_pathway_correlation_comparison.R
#
# Description:
#   Per-compound-pair comparison of gene-level vs. pathway-level similarity.
#   Joins the already-computed similarity matrices from R/compound_similarity.R
#   (gene level: Pearson correlation of log2FoldChange, direction-split
#   Jaccard overlap of significant genes) and R/gsva_compound_correlation.R
#   (pathway level: Pearson correlation of GSVA diff_score, direction-split
#   Jaccard overlap of significant pathways) to identify which compound pairs
#   lose similarity when aggregated from genes into pathways.
#
################################################################################

#===============================================================================
# Run Gene-Level vs. Pathway-Level Compound Similarity Comparison
#===============================================================================

run_gene_vs_pathway_comparison <- function(config) {

    message("Running gene-level vs. pathway-level compound similarity comparison...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    if (!requireNamespace("ggrepel", quietly = TRUE)) {
        stop("Package 'ggrepel' is required.")
    }

    deseq2_table_dir <- file.path(config$output_directory, "DESeq2", "Tables")
    gsva_table_dir <- file.path(config$output_directory, "GSVA", "Tables")
    figure_dir <- file.path(config$output_directory, "GSVA", "Figures")

    #===========================================================================
    # Read Gene-Level and Pathway-Level Similarity Matrices
    #===========================================================================

    read_matrix <- function(path) {
        m <- read.csv(path, row.names = 1, check.names = FALSE)
        as.matrix(m)
    }

    gene_corr <- read_matrix(file.path(deseq2_table_dir, "compound_lfc_correlation.csv"))
    gene_jaccard_up <- read_matrix(file.path(deseq2_table_dir, "compound_jaccard_up.csv"))
    gene_jaccard_down <- read_matrix(file.path(deseq2_table_dir, "compound_jaccard_down.csv"))

    pathway_corr <- read_matrix(file.path(gsva_table_dir, "compound_gsva_pathway_correlation.csv"))
    pathway_jaccard_up <- read_matrix(file.path(gsva_table_dir, "compound_gsva_jaccard_up.csv"))
    pathway_jaccard_down <- read_matrix(file.path(gsva_table_dir, "compound_gsva_jaccard_down.csv"))

    #===========================================================================
    # Restrict to the Common Compound Set
    #===========================================================================

    compounds <- Reduce(
        intersect,
        list(
            rownames(gene_corr), rownames(gene_jaccard_up), rownames(gene_jaccard_down),
            rownames(pathway_corr), rownames(pathway_jaccard_up), rownames(pathway_jaccard_down)
        )
    )

    compounds <- sort(compounds)

    message(sprintf("Comparing %d compounds common to both gene-level and pathway-level results.", length(compounds)))

    #===========================================================================
    # Build the Long Pairwise Comparison Table
    #===========================================================================

    pair_idx <- utils::combn(compounds, 2)

    get_pair <- function(matrix, a, b) matrix[a, b]

    comparison <- data.frame(
        compound_a = pair_idx[1, ],
        compound_b = pair_idx[2, ],
        gene_corr = mapply(get_pair, pair_idx[1, ], pair_idx[2, ], MoreArgs = list(matrix = gene_corr)),
        pathway_corr = mapply(get_pair, pair_idx[1, ], pair_idx[2, ], MoreArgs = list(matrix = pathway_corr)),
        gene_jaccard_up = mapply(get_pair, pair_idx[1, ], pair_idx[2, ], MoreArgs = list(matrix = gene_jaccard_up)),
        pathway_jaccard_up = mapply(get_pair, pair_idx[1, ], pair_idx[2, ], MoreArgs = list(matrix = pathway_jaccard_up)),
        gene_jaccard_down = mapply(get_pair, pair_idx[1, ], pair_idx[2, ], MoreArgs = list(matrix = gene_jaccard_down)),
        pathway_jaccard_down = mapply(get_pair, pair_idx[1, ], pair_idx[2, ], MoreArgs = list(matrix = pathway_jaccard_down)),
        stringsAsFactors = FALSE
    )

    comparison <- comparison[stats::complete.cases(comparison), ]

    comparison$corr_delta <- comparison$gene_corr - comparison$pathway_corr
    comparison$jaccard_up_delta <- comparison$gene_jaccard_up - comparison$pathway_jaccard_up
    comparison$jaccard_down_delta <- comparison$gene_jaccard_down - comparison$pathway_jaccard_down

    comparison <- comparison[order(-comparison$corr_delta), ]

    #===========================================================================
    # Save Table
    #===========================================================================

    dir.create(gsva_table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        comparison,
        file.path(gsva_table_dir, "gene_vs_pathway_similarity_comparison.csv"),
        row.names = FALSE
    )

    #===========================================================================
    # Scatter Plots: Gene-Level vs. Pathway-Level, One Per Metric
    #===========================================================================

    make_comparison_plot <- function(data, x_col, y_col, delta_col, title, x_label, y_label, n_label = 8) {

        label_data <- data[order(-abs(data[[delta_col]])), ][seq_len(min(n_label, nrow(data))), ]
        label_data$pair_label <- paste(label_data$compound_a, label_data$compound_b, sep = " / ")

        plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
            ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
            ggplot2::geom_point(alpha = 0.6) +
            ggrepel::geom_text_repel(
                data = label_data,
                ggplot2::aes(label = pair_label),
                size = 3, max.overlaps = Inf, segment.size = 0.3
            ) +
            ggplot2::labs(title = title, x = x_label, y = y_label) +
            ggplot2::theme_bw()

        plot

    }

    corr_plot <- make_comparison_plot(
        comparison, "gene_corr", "pathway_corr", "corr_delta",
        "Gene-Level vs. Pathway-Level Correlation, All Compound Pairs",
        "Pearson correlation (gene log2FoldChange)",
        "Pearson correlation (pathway diff_score)"
    )

    jaccard_up_plot <- make_comparison_plot(
        comparison, "gene_jaccard_up", "pathway_jaccard_up", "jaccard_up_delta",
        "Gene-Level vs. Pathway-Level Jaccard Overlap (Upregulated), All Compound Pairs",
        "Jaccard overlap (upregulated genes)",
        "Jaccard overlap (upregulated pathways)"
    )

    jaccard_down_plot <- make_comparison_plot(
        comparison, "gene_jaccard_down", "pathway_jaccard_down", "jaccard_down_delta",
        "Gene-Level vs. Pathway-Level Jaccard Overlap (Downregulated), All Compound Pairs",
        "Jaccard overlap (downregulated genes)",
        "Jaccard overlap (downregulated pathways)"
    )

    plot_width <- max(config$figure_width, 11)
    plot_height <- max(config$figure_height, 8)

    save_pair_plot <- function(plot, basename) {

        if (isTRUE(config$save_png)) {
            ggplot2::ggsave(
                file.path(figure_dir, paste0(basename, ".png")),
                plot = plot, width = plot_width, height = plot_height, dpi = config$figure_dpi
            )
        }

        if (isTRUE(config$save_pdf)) {
            ggplot2::ggsave(
                file.path(figure_dir, paste0(basename, ".pdf")),
                plot = plot, width = plot_width, height = plot_height
            )
        }

    }

    save_pair_plot(corr_plot, "Gene_vs_Pathway_correlation_comparison")
    save_pair_plot(jaccard_up_plot, "Gene_vs_Pathway_jaccard_up_comparison")
    save_pair_plot(jaccard_down_plot, "Gene_vs_Pathway_jaccard_down_comparison")

    message("Gene-level vs. pathway-level compound similarity comparison complete.")

    return(
        list(
            comparison = comparison
        )
    )

}
