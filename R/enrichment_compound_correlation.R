################################################################################
# RNAseqPipeline
#
# File: enrichment_compound_correlation.R
#
# Description:
#   Compound x compound similarity at the pathway level, read from the
#   per-contrast, direction-split GO/KEGG over-representation results
#   already written by run_enrichment() (R/enrichment.R):
#   "<GO|KEGG>_<group>_up.csv" and "..._down.csv". Mirrors
#   R/gsea_compound_correlation.R's design (Pearson correlation of a
#   per-pathway score across compounds, each compound's TOP TESTED DOSE
#   only via parse_dose_group()), scored as -log10(p.adjust). Two
#   complementary views per family/direction, mirroring
#   R/gsva_compound_correlation.R:
#     - Pearson correlation of -log10(p.adjust) across the pathway IDs
#       either compound has.
#     - Jaccard overlap of the (already-significant, since every row in
#       an ORA direction file already passed run_enrichment()'s
#       pvalueCutoff) pathway ID sets.
#
#   Two SEPARATE sets of results per pathway family: one over only the
#   up-regulated-gene enrichment results, one over only the
#   down-regulated-gene enrichment results -- not a single combined/signed
#   score. Run separately for GO and KEGG too -- different pathway
#   universes.
#
################################################################################

#===============================================================================
# Run Compound-vs-Compound Enrichment Correlation for One Family + Direction
#===============================================================================

run_enrichment_direction_correlation <- function(config, pathway_family, direction) {

    table_dir <- file.path(config$output_directory, pathway_family, "Tables")
    figure_dir <- file.path(config$output_directory, pathway_family, "Figures")

    #===========================================================================
    # Discover Each Compound's Top Dose from This Direction's Files on Disk
    #===========================================================================

    direction_files <- list.files(
        table_dir,
        pattern = paste0("^", pathway_family, "_.*_", direction, "\\.csv$")
    )

    group_names <- sub(
        paste0("_", direction, "\\.csv$"), "",
        sub(paste0("^", pathway_family, "_"), "", direction_files)
    )

    parsed <- lapply(group_names, parse_dose_group)

    keep <- !vapply(parsed, is.null, logical(1))

    if (!any(keep)) {
        warning(
            "No parseable dose-response groups found for ", pathway_family,
            " (", direction, ") enrichment; skipping."
        )
        return(NULL)
    }

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

    if (length(compounds) < 2) {
        warning(
            "Fewer than 2 compounds have ", pathway_family,
            " (", direction, ") enrichment results; skipping."
        )
        return(NULL)
    }

    message(sprintf(
        "Correlating %d compounds' %s-regulated %s enrichment profiles at their highest tested dose.",
        length(compounds), direction, pathway_family
    ))

    #===========================================================================
    # Read Each Compound's -log10(p.adjust) Pathway Score for This Direction
    #===========================================================================

    direction_file_for <- function(group) {
        file.path(table_dir, paste0(pathway_family, "_", group, "_", direction, ".csv"))
    }

    per_compound_scores <- lapply(compounds, function(cmpd) {

        group <- top_dose_table$group[top_dose_table$compound == cmpd]
        f <- direction_file_for(group)

        if (!file.exists(f)) {
            return(stats::setNames(numeric(0), character(0)))
        }

        results <- utils::read.csv(f)

        if (nrow(results) == 0) {
            return(stats::setNames(numeric(0), character(0)))
        }

        scores <- -log10(results$p.adjust)

        stats::setNames(scores, results$ID)

    })

    names(per_compound_scores) <- compounds

    #===========================================================================
    # Build the Pathway x Compound Matrix, Pearson Correlation
    #===========================================================================

    score_matrix <- Reduce(
        function(a, b) merge(a, b, by = "ID", all = TRUE),
        lapply(compounds, function(cmpd) {
            scores <- per_compound_scores[[cmpd]]
            df <- data.frame(ID = names(scores), stringsAsFactors = FALSE)
            df[[cmpd]] <- as.numeric(scores)
            df
        })
    )

    rownames(score_matrix) <- score_matrix$ID
    score_matrix$ID <- NULL
    score_matrix <- score_matrix[, compounds, drop = FALSE]

    correlation_matrix <- stats::cor(
        as.matrix(score_matrix),
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    #===========================================================================
    # Jaccard Overlap of Significant Pathway ID Sets
    #===========================================================================
    #
    # Every ID in a compound's direction file already passed
    # run_enrichment()'s significance cutoff (pvalueCutoff = config$alpha),
    # so no extra filtering is needed here -- just set overlap.
    #===========================================================================

    pathway_sets <- lapply(per_compound_scores, names)

    jaccard_matrix <- matrix(
        NA_real_,
        nrow = length(compounds), ncol = length(compounds),
        dimnames = list(compounds, compounds)
    )

    for (i in compounds) {
        for (j in compounds) {

            set_union <- length(union(pathway_sets[[i]], pathway_sets[[j]]))

            jaccard_matrix[i, j] <- if (set_union == 0) {
                NA_real_
            } else {
                length(intersect(pathway_sets[[i]], pathway_sets[[j]])) / set_union
            }

        }
    }

    #===========================================================================
    # Save Tables and Heatmaps
    #===========================================================================

    direction_label <- tools::toTitleCase(direction)

    save_csv(
        correlation_matrix,
        file.path(table_dir, paste0("compound_", pathway_family, "_", direction, "_enrichment_correlation.csv"))
    )

    save_csv(
        jaccard_matrix,
        file.path(table_dir, paste0("compound_", pathway_family, "_", direction, "_enrichment_jaccard.csv"))
    )

    heatmap_height <- max(config$figure_height, length(compounds) * 0.5)
    heatmap_width <- heatmap_height + 2

    # ORA only reports significant pathways (unlike GSVA/GSEA's dense
    # per-pathway scores), so two compounds can share zero terms -- that
    # pair's pairwise.complete.obs correlation (and 0/0 Jaccard union) is
    # NA, which hclust() can't cluster on. Fall back to no clustering
    # whenever any NA is present, checked separately per matrix.
    can_cluster_correlation <- !anyNA(correlation_matrix)
    can_cluster_jaccard <- !anyNA(jaccard_matrix)

    if (!can_cluster_correlation) {
        message(
            "Some compound pairs share no significant ", direction, "-regulated ",
            pathway_family, " terms (NA correlation); disabling correlation heatmap clustering."
        )
    }

    if (!can_cluster_jaccard) {
        message(
            "Some compound pairs have no ", direction, "-regulated ",
            pathway_family, " terms at all (NA Jaccard); disabling Jaccard heatmap clustering."
        )
    }

    save_matrix_heatmap <- function(matrix_to_plot, can_cluster, metric_label, filename_suffix) {

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
                main = paste0(
                    "Compound Similarity: ", direction_label, "-Regulated ",
                    pathway_family, " ", metric_label, " (Top Dose)"
                ),
                filename = file.path(
                    figure_dir,
                    paste0("Compound_", pathway_family, "_", direction, "_", filename_suffix, ".", ext)
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

    save_matrix_heatmap(correlation_matrix, can_cluster_correlation, "Enrichment Correlation", "enrichment_correlation_heatmap")
    save_matrix_heatmap(jaccard_matrix, can_cluster_jaccard, "Enrichment Jaccard Overlap", "enrichment_jaccard_heatmap")

    list(
        top_dose_table = top_dose_table,
        correlation_matrix = correlation_matrix,
        jaccard_matrix = jaccard_matrix
    )

}

#===============================================================================
# Run Compound-vs-Compound Enrichment Correlation for One Pathway Family
#===============================================================================

run_enrichment_pathway_family_correlation <- function(config, pathway_family) {

    list(
        up = run_enrichment_direction_correlation(config, pathway_family, "up"),
        down = run_enrichment_direction_correlation(config, pathway_family, "down")
    )

}

#===============================================================================
# Run Compound-vs-Compound Enrichment Correlation (GO + KEGG, Up + Down)
#===============================================================================

run_enrichment_compound_correlation <- function(config) {

    message("Running compound-vs-compound GO/KEGG enrichment correlation...")

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    go_result <- if (isTRUE(config$run_go)) {
        run_enrichment_pathway_family_correlation(config, "GO")
    } else {
        NULL
    }

    kegg_result <- if (isTRUE(config$run_kegg)) {
        run_enrichment_pathway_family_correlation(config, "KEGG")
    } else {
        NULL
    }

    message("Compound-vs-compound GO/KEGG enrichment correlation complete.")

    list(GO = go_result, KEGG = kegg_result)

}
