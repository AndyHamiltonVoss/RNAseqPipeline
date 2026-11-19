################################################################################
# RNAseqPipeline
#
# File: gsea_compound_correlation.R
#
# Description:
#   Compound x compound similarity at the pathway level, read from the
#   per-contrast GSEA GO/KEGG results already written by run_gsea()
#   (R/gsea.R) -- any of four sibling directories via the gsea_dir
#   argument: "GSEA" (stat-ranked, significant-only, the default),
#   "GSEA_LFC" (shrunken-log2FoldChange-ranked, significant-only),
#   "GSEA_Dense" / "GSEA_LFC_Dense" (same two rankings, but EVERY tested
#   pathway, written when config$run_gsea_dense = TRUE -- see below).
#   Mirrors R/gsva_compound_correlation.R's design: each compound's TOP
#   TESTED DOSE only (parse_dose_group()), up to three complementary
#   views:
#     - Pearson correlation of the per-pathway NES profile (GSEA's NES in
#       place of GSVA limma's diff_score).
#     - Jaccard overlap of significant (p.adjust < alpha) pathway ID sets,
#       split by direction (NES > 0 = up, NES < 0 = down).
#     - Leading-edge (core_enrichment) GENE overlap for pathways
#       significant in BOTH compounds of a pair -- a finer-grained check
#       of whether the same genes are actually driving a shared
#       significant pathway, not just whether the pathway ID matches.
#       Gated on config$run_gsea_leading_edge_correlation; skipped for
#       "*_Dense" directories (redundant -- same significant subset as
#       the non-dense sibling).
#
#   Run separately for GO and KEGG -- different pathway universes (GO term
#   IDs vs. KEGG pathway IDs), so a single combined matrix wouldn't make
#   sense. Unlike run_gsva_compound_correlation(), this does not also run
#   the paired up/down Wilcoxon test -- that wasn't asked for here.
#
#   *** CAVEAT on the NES Pearson correlation when gsea_dir is "GSEA" or
#   "GSEA_LFC" (the significant-only directories, unlike GSVA's dense,
#   every-pathway-scored matrix): run_gsea() filters each compound's
#   result file down to its own significant terms (pvalueCutoff =
#   config$alpha) there. So a compound PAIR's correlation is computed only
#   over the pathways significant in BOTH (pairwise.complete.obs on the
#   merged matrix), not the full tested pathway universe. Two failure
#   modes follow from that:
#     - Compounds with MANY significant terms tend to share a large,
#       broadly-co-regulated "generic stress/metabolism" subset with
#       extreme, bimodal NES (~+/-2 to 2.8) -- trivially correlates near 1
#       regardless of true biological similarity (observed: Dar vs. Rosi
#       r=0.998 on 599 shared GO terms, but Dar vs. the much less similar
#       SR1664 still r=0.987 on 465 shared terms).
#     - Compounds with FEW significant terms (e.g. AV2, 16 GO terms total)
#       give tiny, unstable intersections -- e.g. AV2 vs. Dar lands at
#       r=-0.999 off just 13 points, not a real anti-correlation.
#   Pointed at "GSEA_Dense"/"GSEA_LFC_Dense" instead, this problem goes
#   away: every compound has an NES for every tested pathway, so the
#   correlation is apples-to-apples like GSVA's. The Jaccard overlap does
#   not have this distortion either way (set overlap isn't skewed by
#   intersection size the same way) -- use a Dense directory for a sound
#   NES correlation, or treat a significant-only directory's correlation
#   as a rough, lower-confidence companion to its Jaccard.
#
################################################################################

#===============================================================================
# Run Compound-vs-Compound GSEA Correlation for One Pathway Family (GO/KEGG)
#===============================================================================

run_gsea_pathway_family_correlation <- function(config, pathway_family, gsea_dir = "GSEA") {

    table_dir <- file.path(config$output_directory, gsea_dir, pathway_family, "Tables")
    figure_dir <- file.path(config$output_directory, gsea_dir, pathway_family, "Figures")

    # run_gsea()'s dense pass (gsea_dir = "GSEA_Dense"/"GSEA_LFC_Dense")
    # only creates its Tables dir, not Figures (no dotplot is written
    # there) -- create figure_dir here too so the heatmaps below have
    # somewhere to write.
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    #===========================================================================
    # Discover Each Compound's Top Dose from the GSEA Result Files on Disk
    #===========================================================================

    result_files <- list.files(
        table_dir,
        pattern = paste0("^GSEA_", pathway_family, "_.*\\.csv$")
    )

    group_names <- sub(
        "\\.csv$", "",
        sub(paste0("^GSEA_", pathway_family, "_"), "", result_files)
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

    message(sprintf(
        "Correlating %d compounds' %s %s NES profiles at their highest tested dose.",
        length(compounds), gsea_dir, pathway_family
    ))

    #===========================================================================
    # Read Each Compound's NES Profile, Build the Pathway x Compound Matrix
    #===========================================================================

    gsea_file_for <- function(group) {
        file.path(table_dir, paste0("GSEA_", pathway_family, "_", group, ".csv"))
    }

    per_compound <- lapply(compounds, function(cmpd) {

        group <- top_dose_table$group[top_dose_table$compound == cmpd]
        results <- utils::read.csv(gsea_file_for(group))

        sig <- !is.na(results$p.adjust) & results$p.adjust < config$alpha

        list(
            nes = stats::setNames(results$NES, results$ID),
            up = results$ID[sig & results$NES > 0],
            down = results$ID[sig & results$NES < 0],
            core_enrichment = stats::setNames(results$core_enrichment, results$ID),
            description = stats::setNames(results$Description, results$ID)
        )

    })

    names(per_compound) <- compounds

    nes_matrix <- Reduce(
        function(a, b) merge(a, b, by = "ID", all = TRUE),
        lapply(compounds, function(cmpd) {
            df <- data.frame(ID = names(per_compound[[cmpd]]$nes), stringsAsFactors = FALSE)
            df[[cmpd]] <- as.numeric(per_compound[[cmpd]]$nes)
            df
        })
    )

    rownames(nes_matrix) <- nes_matrix$ID
    nes_matrix$ID <- NULL
    nes_matrix <- nes_matrix[, compounds, drop = FALSE]

    correlation_matrix <- stats::cor(
        as.matrix(nes_matrix),
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
    # Save Tables and Heatmaps
    #===========================================================================

    save_csv(
        correlation_matrix,
        file.path(table_dir, paste0("compound_gsea_", pathway_family, "_nes_correlation.csv"))
    )

    save_csv(
        jaccard_up_matrix,
        file.path(table_dir, paste0("compound_gsea_", pathway_family, "_jaccard_up.csv"))
    )

    save_csv(
        jaccard_down_matrix,
        file.path(table_dir, paste0("compound_gsea_", pathway_family, "_jaccard_down.csv"))
    )

    heatmap_height <- max(config$figure_height, length(compounds) * 0.5)
    heatmap_width <- heatmap_height + 2

    # Compound pairs can share zero significant pathways in a direction
    # (NA Jaccard) or zero overlapping pathways at all (NA correlation),
    # which hclust() can't cluster on. Fall back to no clustering whenever
    # any NA is present, checked separately per matrix.
    can_cluster_correlation <- !anyNA(correlation_matrix)
    can_cluster_jaccard_up <- !anyNA(jaccard_up_matrix)
    can_cluster_jaccard_down <- !anyNA(jaccard_down_matrix)

    if (!can_cluster_correlation) {
        message(
            "Some compound pairs share no ", pathway_family,
            " pathways in common (NA correlation); disabling correlation heatmap clustering."
        )
    }

    if (!can_cluster_jaccard_up || !can_cluster_jaccard_down) {
        message(
            "Some compound pairs have no significant up- or down-NES ", pathway_family,
            " pathways at all (NA Jaccard); disabling Jaccard heatmap clustering where this occurs."
        )
    }

    save_matrix_heatmap <- function(matrix_to_plot, can_cluster, title, filename_suffix) {

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
                filename = file.path(
                    figure_dir,
                    paste0("Compound_GSEA_", pathway_family, "_", filename_suffix, ".", ext)
                ),
                width = heatmap_width,
                height = heatmap_height
            )

            if (ext == "png") {
                args$res <- config$figure_dpi
            }

            do.call(pheatmap::pheatmap, args)

        }

    }

    gsea_label <- switch(
        gsea_dir,
        "GSEA" = "GSEA",
        "GSEA_LFC" = "GSEA (LFC-Ranked)",
        "GSEA_Dense" = "GSEA (Dense, All Tested Pathways)",
        "GSEA_LFC_Dense" = "GSEA (LFC-Ranked, Dense)"
    )

    save_matrix_heatmap(
        correlation_matrix, can_cluster_correlation,
        paste0("Compound Similarity: ", gsea_label, " ", pathway_family, " NES Correlation (Top Dose)"),
        "nes_correlation_heatmap"
    )

    save_matrix_heatmap(
        jaccard_up_matrix, can_cluster_jaccard_up,
        paste0("Compound Similarity: Upregulated ", gsea_label, " ", pathway_family, " Jaccard Overlap (Top Dose)"),
        "jaccard_up_heatmap"
    )

    save_matrix_heatmap(
        jaccard_down_matrix, can_cluster_jaccard_down,
        paste0("Compound Similarity: Downregulated ", gsea_label, " ", pathway_family, " Jaccard Overlap (Top Dose)"),
        "jaccard_down_heatmap"
    )

    #===========================================================================
    # Pathway x Compound NES Heatmap (Top Variable Pathways)
    #===========================================================================
    #
    # The nes_matrix built above (pathways x compounds) is otherwise only
    # used to compute correlation_matrix -- plot it directly too, top-N
    # most-variable pathways by NES across compounds, mirroring
    # save_gsva_heatmap() in R/gsva.R. Most useful (no/few NA cells) on a
    # "*_Dense" directory, where every compound has an NES for every
    # tested pathway; on a significant-only directory most pathways are
    # NA for most compounds, so this still renders but is sparser.
    #===========================================================================

    pathway_vars <- apply(nes_matrix, 1, stats::var, na.rm = TRUE)
    pathway_vars <- pathway_vars[!is.na(pathway_vars) & is.finite(pathway_vars)]

    n_nes_heatmap <- min(config$gsea_nes_heatmap_top_pathways, length(pathway_vars))

    if (n_nes_heatmap > 0) {

        heatmap_ids <- names(sort(pathway_vars, decreasing = TRUE))[seq_len(n_nes_heatmap)]
        nes_heatmap_matrix <- as.matrix(nes_matrix[heatmap_ids, , drop = FALSE])

        description_lookup <- character(0)
        for (cmpd in compounds) {
            new_ids <- setdiff(names(per_compound[[cmpd]]$description), names(description_lookup))
            if (length(new_ids) > 0) {
                description_lookup[new_ids] <- per_compound[[cmpd]]$description[new_ids]
            }
        }

        rownames(nes_heatmap_matrix) <- ifelse(
            heatmap_ids %in% names(description_lookup),
            description_lookup[heatmap_ids],
            heatmap_ids
        )

        save_csv(
            nes_heatmap_matrix,
            file.path(table_dir, paste0("compound_gsea_", pathway_family, "_nes_top_variable_pathways.csv"))
        )

        max_abs_nes <- max(abs(nes_heatmap_matrix), na.rm = TRUE)
        nes_breaks <- seq(-max_abs_nes, max_abs_nes, length.out = 101)
        nes_palette <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

        nes_heatmap_height <- max(config$figure_height, 3 + n_nes_heatmap * 0.3)
        nes_heatmap_width <- max(config$figure_width, length(compounds) * 0.4) + 4

        for (ext in c("png", "pdf")) {

            save_this_ext <- if (ext == "png") isTRUE(config$save_png) else isTRUE(config$save_pdf)

            if (!save_this_ext) {
                next
            }

            args <- list(
                nes_heatmap_matrix,
                color = nes_palette,
                breaks = nes_breaks,
                cluster_rows = isTRUE(config$cluster_rows),
                cluster_cols = isTRUE(config$cluster_columns),
                na_col = "grey80",
                main = paste0("Top Variable ", gsea_label, " ", pathway_family, " Pathway NES (Top Dose)"),
                filename = file.path(
                    figure_dir,
                    paste0("Heatmap_GSEA_", pathway_family, "_nes_top_variable_pathways.", ext)
                ),
                width = nes_heatmap_width,
                height = nes_heatmap_height
            )

            if (ext == "png") {
                args$res <- config$figure_dpi
            }

            do.call(pheatmap::pheatmap, args)

        }

    }

    #===========================================================================
    # Pathway-Level Compound PCA
    #===========================================================================
    #
    # stats::prcomp() on the compound x pathway NES matrix (transposed
    # nes_matrix), mirroring R/qc.R's sample-level PCA (same
    # pca_center/pca_scale config) but at the compound level and on
    # pathway NES scores instead of gene expression. Most meaningful on a
    # "*_Dense" directory (complete NES profile, no NA pathways); on a
    # significant-only directory, NA cells (pathway not significant for
    # that compound) are imputed to 0 ("no enrichment signal") since
    # prcomp() cannot handle NA input.
    #===========================================================================

    if (length(compounds) >= 4) {

        pca_input <- as.matrix(nes_matrix)
        pca_input[is.na(pca_input)] <- 0

        pca <- stats::prcomp(
            t(pca_input),
            center = isTRUE(config$pca_center),
            scale. = isTRUE(config$pca_scale)
        )

        variance_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100

        pca_scores <- as.data.frame(pca$x)
        pca_scores$compound <- rownames(pca_scores)

        save_csv(
            pca_scores,
            file.path(table_dir, paste0("compound_gsea_", pathway_family, "_pathway_pca.csv")),
            row.names = FALSE
        )

        save_csv(
            data.frame(
                PC = paste0("PC", seq_along(variance_explained)),
                percent_variance = variance_explained
            ),
            file.path(table_dir, paste0("compound_gsea_", pathway_family, "_pathway_pca_variance_explained.csv")),
            row.names = FALSE
        )

        if (requireNamespace("ggplot2", quietly = TRUE)) {

            pca_plot <- ggplot2::ggplot(pca_scores, ggplot2::aes(x = PC1, y = PC2)) +
                ggplot2::geom_point(size = 3, alpha = 0.85, color = "steelblue") +
                ggplot2::labs(
                    title = paste0("Compound PCA: ", gsea_label, " ", pathway_family, " Pathway NES Profile (Top Dose)"),
                    x = sprintf("PC1 (%.1f%% variance)", variance_explained[1]),
                    y = sprintf("PC2 (%.1f%% variance)", variance_explained[2])
                ) +
                ggplot2::theme_bw()

            if (requireNamespace("ggrepel", quietly = TRUE)) {
                pca_plot <- pca_plot + ggrepel::geom_text_repel(ggplot2::aes(label = compound), size = 3.2)
            } else {
                pca_plot <- pca_plot + ggplot2::geom_text(ggplot2::aes(label = compound), size = 3.2, vjust = -1)
            }

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
                    title = paste0(gsea_label, " ", pathway_family, " Pathway PCA Scree Plot"),
                    x = "Principal Component",
                    y = "Percent Variance Explained"
                ) +
                ggplot2::theme_bw()

            if (isTRUE(config$save_png)) {

                ggplot2::ggsave(
                    file.path(figure_dir, paste0("Compound_GSEA_", pathway_family, "_pathway_pca.png")),
                    plot = pca_plot, width = config$figure_width, height = config$figure_height, dpi = config$figure_dpi
                )

                ggplot2::ggsave(
                    file.path(figure_dir, paste0("Compound_GSEA_", pathway_family, "_pathway_pca_scree.png")),
                    plot = scree_plot, width = config$figure_width, height = config$figure_height, dpi = config$figure_dpi
                )

            }

            if (isTRUE(config$save_pdf)) {

                ggplot2::ggsave(
                    file.path(figure_dir, paste0("Compound_GSEA_", pathway_family, "_pathway_pca.pdf")),
                    plot = pca_plot, width = config$figure_width, height = config$figure_height
                )

                ggplot2::ggsave(
                    file.path(figure_dir, paste0("Compound_GSEA_", pathway_family, "_pathway_pca_scree.pdf")),
                    plot = scree_plot, width = config$figure_width, height = config$figure_height
                )

            }

        }

    } else {

        pca_scores <- NULL

    }

    #===========================================================================
    # Leading-Edge (core_enrichment) Gene Overlap for Shared Significant
    # Pathways
    #===========================================================================
    #
    # Two compounds can both call the same pathway ID significant without
    # the same genes actually driving it. For every pathway ID significant
    # in BOTH compounds of a pair (already computed above as the Jaccard
    # intersection), split clusterProfiler's "/"-delimited core_enrichment
    # string and take the gene-level Jaccard overlap -- a finer-grained,
    # independent signal on top of the pathway-ID Jaccard, not a subset of
    # it. Skipped for "*_Dense" directories: a dense directory's
    # significant subset is identical to its non-dense sibling's, so this
    # would just duplicate the same files.
    #===========================================================================

    leading_edge_up_matrix <- NULL
    leading_edge_down_matrix <- NULL

    if (isTRUE(config$run_gsea_leading_edge_correlation) && !grepl("_Dense$", gsea_dir)) {

        gene_set_jaccard <- function(a, b) {

            genes_a <- unique(strsplit(a, "/", fixed = TRUE)[[1]])
            genes_b <- unique(strsplit(b, "/", fixed = TRUE)[[1]])

            union_n <- length(union(genes_a, genes_b))

            list(
                jaccard = if (union_n == 0) NA_real_ else length(intersect(genes_a, genes_b)) / union_n,
                n_a = length(genes_a),
                n_b = length(genes_b)
            )

        }

        pair_idx <- utils::combn(compounds, 2)
        detail_rows <- list()

        for (k in seq_len(ncol(pair_idx))) {

            a <- pair_idx[1, k]
            b <- pair_idx[2, k]

            for (direction in c("up", "down")) {

                shared <- intersect(per_compound[[a]][[direction]], per_compound[[b]][[direction]])

                for (pid in shared) {

                    g <- gene_set_jaccard(per_compound[[a]]$core_enrichment[[pid]], per_compound[[b]]$core_enrichment[[pid]])

                    detail_rows[[length(detail_rows) + 1]] <- data.frame(
                        compound_a = a, compound_b = b, direction = direction,
                        pathway_ID = pid, Description = per_compound[[a]]$description[[pid]],
                        leading_edge_jaccard = g$jaccard, n_core_genes_a = g$n_a, n_core_genes_b = g$n_b,
                        stringsAsFactors = FALSE
                    )

                }

            }

        }

        leading_edge_detail <- if (length(detail_rows) > 0) do.call(rbind, detail_rows) else data.frame(
            compound_a = character(0), compound_b = character(0), direction = character(0),
            pathway_ID = character(0), Description = character(0),
            leading_edge_jaccard = numeric(0), n_core_genes_a = integer(0), n_core_genes_b = integer(0)
        )

        save_csv(
            leading_edge_detail,
            file.path(table_dir, paste0("compound_gsea_", pathway_family, "_leading_edge_overlap.csv")),
            row.names = FALSE
        )

        leading_edge_jaccard_matrix <- function(direction) {

            m <- matrix(NA_real_, nrow = length(compounds), ncol = length(compounds), dimnames = list(compounds, compounds))

            for (i in compounds) {
                for (j in compounds) {

                    if (i == j) {
                        m[i, j] <- if (length(per_compound[[i]][[direction]]) == 0) NA_real_ else 1
                        next
                    }

                    rows <- leading_edge_detail[
                        leading_edge_detail$direction == direction &
                            ((leading_edge_detail$compound_a == i & leading_edge_detail$compound_b == j) |
                                (leading_edge_detail$compound_a == j & leading_edge_detail$compound_b == i)),
                    ]

                    if (nrow(rows) > 0) {
                        m[i, j] <- mean(rows$leading_edge_jaccard, na.rm = TRUE)
                    }

                }
            }

            m

        }

        leading_edge_up_matrix <- leading_edge_jaccard_matrix("up")
        leading_edge_down_matrix <- leading_edge_jaccard_matrix("down")

        save_csv(leading_edge_up_matrix, file.path(table_dir, paste0("compound_gsea_", pathway_family, "_leading_edge_jaccard_up.csv")))
        save_csv(leading_edge_down_matrix, file.path(table_dir, paste0("compound_gsea_", pathway_family, "_leading_edge_jaccard_down.csv")))

        save_matrix_heatmap(
            leading_edge_up_matrix, !anyNA(leading_edge_up_matrix),
            paste0("Compound Similarity: ", gsea_label, " ", pathway_family, " Leading-Edge Gene Jaccard, Upregulated (Top Dose)"),
            "leading_edge_jaccard_up_heatmap"
        )

        save_matrix_heatmap(
            leading_edge_down_matrix, !anyNA(leading_edge_down_matrix),
            paste0("Compound Similarity: ", gsea_label, " ", pathway_family, " Leading-Edge Gene Jaccard, Downregulated (Top Dose)"),
            "leading_edge_jaccard_down_heatmap"
        )

    }

    list(
        top_dose_table = top_dose_table,
        correlation_matrix = correlation_matrix,
        jaccard_up_matrix = jaccard_up_matrix,
        jaccard_down_matrix = jaccard_down_matrix,
        pca_scores = pca_scores,
        leading_edge_up_matrix = leading_edge_up_matrix,
        leading_edge_down_matrix = leading_edge_down_matrix
    )

}

#===============================================================================
# Run Compound-vs-Compound GSEA Pathway Correlation (GO + KEGG)
#===============================================================================

run_gsea_compound_correlation <- function(config, gsea_dir = c("GSEA", "GSEA_LFC", "GSEA_Dense", "GSEA_LFC_Dense")) {

    gsea_dir <- match.arg(gsea_dir)

    message("Running compound-vs-compound ", gsea_dir, " pathway correlation...")

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    go_result <- run_gsea_pathway_family_correlation(config, "GO", gsea_dir = gsea_dir)
    kegg_result <- run_gsea_pathway_family_correlation(config, "KEGG", gsea_dir = gsea_dir)

    message("Compound-vs-compound ", gsea_dir, " pathway correlation complete.")

    list(GO = go_result, KEGG = kegg_result)

}
