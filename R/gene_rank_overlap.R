################################################################################
# RNAseqPipeline
#
# File: gene_rank_overlap.R
#
# Description:
#   Cross-compound consistency of the top-ranked up- and down-regulated
#   genes, by RANK (top-N) rather than by significance threshold.
#   R/compound_similarity.R already found, using a significance cutoff
#   (padj < alpha & |log2FoldChange| > cutoff), that up-regulated gene
#   sets overlap much more across compounds than down-regulated ones
#   ("up-Jaccard >> down-Jaccard"). This asks the same question a
#   different way: a rank-based cutoff doesn't depend on statistical
#   power, so it isn't itself distorted by the weak-agonist noise problem
#   quantified separately this session (some compounds, e.g. SR1664, have
#   up to ~23% of all genes essentially tied near log2FoldChange = 0;
#   strong agonists like Rosi have under 3%).
#
#   There's no single "right" N for a top-N cutoff -- how far into the
#   ranked list a compound's real signal extends before hitting its own
#   noise floor differs by compound (exactly what the tie-fraction check
#   found). So rather than picking one N, this sweeps a range and plots
#   mean pairwise Jaccard (up and down, separately) vs. N -- the TREND
#   across N is the primary, threshold-agnostic result. A single
#   representative N still gets a full compound x compound Jaccard
#   heatmap, for a concrete, inspectable artifact.
#
################################################################################

#===============================================================================
# Build the Top-Dose Gene x Compound log2FoldChange Matrix
#===============================================================================
#
# Same top-dose discovery pattern used throughout this session
# (parse_dose_group() + which.max(dose_nM)), reading each compound's
# apeglm-shrunk DESeq2 log2FoldChange -- same data source as
# R/gene_expression_skew.R.
#===============================================================================

build_top_dose_gene_lfc_matrix <- function(config) {

    deseq2_table_dir <- file.path(config$output_directory, "DESeq2", "Tables")

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

    compounds <- sort(top_dose_table$compound)

    per_compound_lfc <- lapply(compounds, function(cmpd) {
        group <- top_dose_table$group[top_dose_table$compound == cmpd]
        f <- file.path(deseq2_table_dir, paste0("DESeq2_group_", group, "_vs_", config$reference_group, "_shrunk.csv"))
        res <- utils::read.csv(f)
        stats::setNames(res$log2FoldChange, res$gene)
    })

    names(per_compound_lfc) <- compounds

    lfc_matrix <- Reduce(
        function(a, b) merge(a, b, by = "gene", all = TRUE),
        lapply(compounds, function(cmpd) {
            df <- data.frame(gene = names(per_compound_lfc[[cmpd]]), stringsAsFactors = FALSE)
            df[[cmpd]] <- as.numeric(per_compound_lfc[[cmpd]])
            df
        })
    )

    genes <- lfc_matrix$gene
    lfc_matrix$gene <- NULL
    lfc_matrix <- as.matrix(lfc_matrix[, compounds, drop = FALSE])
    rownames(lfc_matrix) <- genes

    lfc_matrix

}

#===============================================================================
# Compound x Compound Jaccard Overlap of Top-N Ranked Genes
#===============================================================================

compute_top_n_jaccard <- function(lfc_matrix, n, direction = c("up", "down")) {

    direction <- match.arg(direction)

    compounds <- colnames(lfc_matrix)

    top_n_genes <- lapply(compounds, function(cmpd) {

        values <- lfc_matrix[, cmpd]
        values <- values[!is.na(values)]

        ordered <- if (direction == "up") {
            sort(values, decreasing = TRUE)
        } else {
            sort(values, decreasing = FALSE)
        }

        names(ordered)[seq_len(min(n, length(ordered)))]

    })

    names(top_n_genes) <- compounds

    jaccard_matrix <- matrix(
        NA_real_,
        nrow = length(compounds), ncol = length(compounds),
        dimnames = list(compounds, compounds)
    )

    for (i in compounds) {
        for (j in compounds) {

            union_n <- length(union(top_n_genes[[i]], top_n_genes[[j]]))

            jaccard_matrix[i, j] <- if (union_n == 0) {
                NA_real_
            } else {
                length(intersect(top_n_genes[[i]], top_n_genes[[j]])) / union_n
            }

        }
    }

    jaccard_matrix

}

#===============================================================================
# Run Gene Rank Overlap Diagnostic
#===============================================================================

run_gene_rank_overlap <- function(config) {

    message("Running gene rank overlap diagnostic...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    lfc_matrix <- build_top_dose_gene_lfc_matrix(config)
    compounds <- colnames(lfc_matrix)

    diagnostic_dir <- file.path(config$output_directory, "gene_diagnostic")
    dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

    #===========================================================================
    # Sweep N: Mean Pairwise Jaccard (Up vs. Down) vs. N
    #===========================================================================

    n_sweep <- config$gene_rank_overlap_n_sweep

    mean_offdiag <- function(m) {
        diag(m) <- NA
        mean(m, na.rm = TRUE)
    }

    sweep_rows <- list()

    for (n in n_sweep) {

        up_matrix <- compute_top_n_jaccard(lfc_matrix, n, "up")
        down_matrix <- compute_top_n_jaccard(lfc_matrix, n, "down")

        sweep_rows[[length(sweep_rows) + 1]] <- data.frame(
            n = n, direction = "Up", mean_jaccard = mean_offdiag(up_matrix), stringsAsFactors = FALSE
        )

        sweep_rows[[length(sweep_rows) + 1]] <- data.frame(
            n = n, direction = "Down", mean_jaccard = mean_offdiag(down_matrix), stringsAsFactors = FALSE
        )

    }

    sweep_table <- do.call(rbind, sweep_rows)

    save_csv(
        sweep_table,
        file.path(diagnostic_dir, "gene_rank_overlap_jaccard_vs_n.csv"),
        row.names = FALSE
    )

    sweep_plot <- ggplot2::ggplot(
        sweep_table,
        ggplot2::aes(x = n, y = mean_jaccard, color = direction)
    ) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_color_manual(values = c(Up = "darkorange", Down = "steelblue")) +
        ggplot2::labs(
            title = "Cross-Compound Jaccard Overlap of Top-N Ranked Genes vs. N",
            subtitle = "Top-dose log2FoldChange ranking; mean pairwise Jaccard across all 20 compounds, by direction",
            x = "N (top-ranked genes, log10 scale)",
            y = "Mean pairwise Jaccard overlap",
            color = "Direction"
        ) +
        ggplot2::theme_bw()

    if (isTRUE(config$save_png)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Gene_rank_overlap_jaccard_vs_n.png"),
            plot = sweep_plot, width = config$figure_width, height = config$figure_height, dpi = config$figure_dpi
        )
    }

    if (isTRUE(config$save_pdf)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Gene_rank_overlap_jaccard_vs_n.pdf"),
            plot = sweep_plot, width = config$figure_width, height = config$figure_height
        )
    }

    #===========================================================================
    # Representative N: Full Compound x Compound Jaccard Heatmaps
    #===========================================================================

    rep_n <- config$gene_rank_overlap_n

    up_matrix <- compute_top_n_jaccard(lfc_matrix, rep_n, "up")
    down_matrix <- compute_top_n_jaccard(lfc_matrix, rep_n, "down")

    save_csv(up_matrix, file.path(diagnostic_dir, paste0("gene_rank_overlap_up_jaccard_n", rep_n, ".csv")))
    save_csv(down_matrix, file.path(diagnostic_dir, paste0("gene_rank_overlap_down_jaccard_n", rep_n, ".csv")))

    heatmap_height <- max(config$figure_height, length(compounds) * 0.5)
    heatmap_width <- heatmap_height + 2

    can_cluster_up <- !anyNA(up_matrix)
    can_cluster_down <- !anyNA(down_matrix)

    save_matrix_heatmap <- function(matrix_to_plot, can_cluster, title, filename_base) {

        for (ext in c("png", "pdf")) {

            save_this_ext <- if (ext == "png") isTRUE(config$save_png) else isTRUE(config$save_pdf)

            if (!save_this_ext) {
                next
            }

            args <- list(
                matrix_to_plot,
                cluster_rows = can_cluster && isTRUE(config$cluster_rows),
                cluster_cols = can_cluster && isTRUE(config$cluster_columns),
                display_numbers = TRUE,
                number_format = "%.2f",
                main = title,
                filename = file.path(diagnostic_dir, paste0(filename_base, ".", ext)),
                width = heatmap_width,
                height = heatmap_height
            )

            if (ext == "png") {
                args$res <- config$figure_dpi
            }

            do.call(pheatmap::pheatmap, args)

        }

    }

    save_matrix_heatmap(
        up_matrix, can_cluster_up,
        sprintf("Compound Similarity: Jaccard Overlap of Top %d Up-Ranked Genes (Top Dose)", rep_n),
        paste0("Gene_rank_overlap_up_heatmap_n", rep_n, "_heatmap")
    )

    save_matrix_heatmap(
        down_matrix, can_cluster_down,
        sprintf("Compound Similarity: Jaccard Overlap of Top %d Down-Ranked Genes (Top Dose)", rep_n),
        paste0("Gene_rank_overlap_down_heatmap_n", rep_n, "_heatmap")
    )

    message("Gene rank overlap diagnostic complete.")

    list(
        sweep_table = sweep_table,
        sweep_plot = sweep_plot,
        up_matrix_representative = up_matrix,
        down_matrix_representative = down_matrix,
        representative_n = rep_n
    )

}
