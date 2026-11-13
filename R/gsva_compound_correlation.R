################################################################################
# RNAseqPipeline
#
# File: gsva_compound_correlation.R
#
# Description:
#   Compound x compound similarity at the PATHWAY level, read from the
#   per-contrast limma results already written by run_gsva() (R/gsva.R).
#   Two complementary views, mirroring R/compound_similarity.R's gene-level
#   analysis:
#     - Pearson correlation of full-pathway-set diff_score.
#     - Jaccard overlap of significant (adj.P.Val < alpha) pathways,
#       direction-aware and split by direction (up-only / down-only), plus
#       a paired Wilcoxon test for the up/down asymmetry.
#
#   Every GSVA limma result file scores the same fixed set of KEGG pathways,
#   so (unlike the gene-level case) there's no ragged/missing-gene handling
#   needed -- this is a plain per-pathway join across compounds.
#
################################################################################

#===============================================================================
# Run Compound-vs-Compound GSVA Pathway Correlation
#===============================================================================

run_gsva_compound_correlation <- function(config) {

    message("Running compound-vs-compound GSVA pathway correlation...")

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    deseq2_table_dir <- file.path(config$output_directory, "DESeq2", "Tables")
    gsva_table_dir <- file.path(config$output_directory, "GSVA", "Tables")
    figure_dir <- file.path(config$output_directory, "GSVA", "Figures")

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

    #===========================================================================
    # Read Each Eligible Compound's GSVA Pathway diff_score Profile
    #===========================================================================

    gsva_file_for <- function(group) {
        file.path(
            gsva_table_dir,
            paste0("GSVA_KEGG_", group, "_vs_", config$reference_group, "_limma.csv")
        )
    }

    top_dose_table$has_gsva <- vapply(
        top_dose_table$group,
        function(group) file.exists(gsva_file_for(group)),
        logical(1)
    )

    excluded_compounds <- top_dose_table$compound[!top_dose_table$has_gsva]

    message(sprintf(
        "%d of %d compounds excluded (no GSVA limma result, e.g. matched config$exclude_groups_pattern): %s",
        length(excluded_compounds), nrow(top_dose_table),
        paste(excluded_compounds, collapse = ", ")
    ))

    eligible_table <- top_dose_table[top_dose_table$has_gsva, ]

    if (nrow(eligible_table) < 2) {
        stop("Fewer than 2 compounds have GSVA results; nothing to correlate.")
    }

    compounds <- sort(eligible_table$compound)

    message(sprintf("Correlating %d compounds' GSVA pathway profiles at their highest tested dose.", length(compounds)))

    #===========================================================================
    # Read Each Compound's Pathway Profile, Classify Up/Down Significant Sets
    #===========================================================================

    per_compound <- lapply(compounds, function(cmpd) {

        group <- eligible_table$group[eligible_table$compound == cmpd]
        results <- read.csv(gsva_file_for(group))

        sig <- !is.na(results$adj.P.Val) & results$adj.P.Val < config$alpha

        list(
            up = results$pathway[sig & results$diff_score > 0],
            down = results$pathway[sig & results$diff_score < 0],
            diff_score = stats::setNames(results$diff_score, results$pathway)
        )

    })

    names(per_compound) <- compounds

    #===========================================================================
    # Build the Pathway x Compound diff_score Matrix, Pearson Correlation
    #===========================================================================

    diff_score_matrix <- Reduce(
        function(a, b) merge(a, b, by = "pathway", all = TRUE),
        lapply(compounds, function(cmpd) {
            df <- data.frame(pathway = names(per_compound[[cmpd]]$diff_score), stringsAsFactors = FALSE)
            df[[cmpd]] <- as.numeric(per_compound[[cmpd]]$diff_score)
            df
        })
    )

    rownames(diff_score_matrix) <- diff_score_matrix$pathway
    diff_score_matrix$pathway <- NULL
    diff_score_matrix <- diff_score_matrix[, compounds, drop = FALSE]

    correlation_matrix <- stats::cor(
        as.matrix(diff_score_matrix),
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    #===========================================================================
    # Jaccard Overlap of Significant Pathways, Split by Direction
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
    # Save Tables and Heatmap
    #===========================================================================

    dir.create(gsva_table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(correlation_matrix, file.path(gsva_table_dir, "compound_gsva_pathway_correlation.csv"))

    save_csv(jaccard_up_matrix, file.path(gsva_table_dir, "compound_gsva_jaccard_up.csv"))

    save_csv(jaccard_down_matrix, file.path(gsva_table_dir, "compound_gsva_jaccard_down.csv"))

    heatmap_height <- max(config$figure_height, length(compounds) * 0.5)
    heatmap_width <- heatmap_height + 2

    if (isTRUE(config$save_png)) {

        pheatmap::pheatmap(
            correlation_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: GSVA KEGG Pathway diff_score Correlation (Top Dose)",
            filename = file.path(figure_dir, "Compound_GSVA_pathway_correlation_heatmap.png"),
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
            main = "Compound Similarity: Jaccard Overlap of Upregulated Pathways (Top Dose)",
            filename = file.path(figure_dir, "Compound_GSVA_jaccard_up_heatmap.png"),
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
            main = "Compound Similarity: Jaccard Overlap of Downregulated Pathways (Top Dose)",
            filename = file.path(figure_dir, "Compound_GSVA_jaccard_down_heatmap.png"),
            width = heatmap_width,
            height = heatmap_height,
            res = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        pheatmap::pheatmap(
            correlation_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: GSVA KEGG Pathway diff_score Correlation (Top Dose)",
            filename = file.path(figure_dir, "Compound_GSVA_pathway_correlation_heatmap.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

        pheatmap::pheatmap(
            jaccard_up_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: Jaccard Overlap of Upregulated Pathways (Top Dose)",
            filename = file.path(figure_dir, "Compound_GSVA_jaccard_up_heatmap.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

        pheatmap::pheatmap(
            jaccard_down_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "Compound Similarity: Jaccard Overlap of Downregulated Pathways (Top Dose)",
            filename = file.path(figure_dir, "Compound_GSVA_jaccard_down_heatmap.pdf"),
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
            title = "Upregulated vs. Downregulated Pathway Overlap, All Compound Pairs",
            subtitle = sprintf(
                "Wilcoxon signed-rank p = %.3g; mean up = %.3f, mean down = %.3f (n = %d pairs)",
                wilcox_test$p.value, mean(pair_comparison$jaccard_up),
                mean(pair_comparison$jaccard_down), nrow(pair_comparison)
            ),
            x = "Jaccard overlap (downregulated pathways)",
            y = "Jaccard overlap (upregulated pathways)"
        ) +
        ggplot2::theme_bw()

    save_csv(
        pair_comparison,
        file.path(gsva_table_dir, "compound_gsva_jaccard_up_vs_down_pairs.csv"),
        row.names = FALSE
    )

    pair_plot_width <- max(config$figure_width, 11)
    pair_plot_height <- max(config$figure_height, 8)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "Compound_GSVA_jaccard_up_vs_down_pairs.png"),
            plot = pair_plot,
            width = pair_plot_width,
            height = pair_plot_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "Compound_GSVA_jaccard_up_vs_down_pairs.pdf"),
            plot = pair_plot,
            width = pair_plot_width,
            height = pair_plot_height
        )

    }

    message("Compound-vs-compound GSVA pathway correlation complete.")

    return(
        list(
            top_dose_table = eligible_table,
            correlation_matrix = correlation_matrix,
            jaccard_up_matrix = jaccard_up_matrix,
            jaccard_down_matrix = jaccard_down_matrix,
            pair_comparison = pair_comparison,
            wilcox_test = wilcox_test
        )
    )

}
