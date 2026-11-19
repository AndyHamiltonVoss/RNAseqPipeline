################################################################################
# RNAseqPipeline
#
# File: pathway_consensus_correlation.R
#
# Description:
#   Right now GSVA (R/gsva_compound_correlation.R), GSEA
#   (R/gsea_compound_correlation.R), and GO/KEGG over-representation
#   (R/enrichment_compound_correlation.R) each give a separate, siloed
#   answer to "how similar are compound A and B at the pathway level" --
#   never combined. This module reads their ALREADY-WRITTEN Jaccard CSVs
#   from disk (no recomputation) and combines them into one consensus view
#   per pathway family (GO/KEGG) x direction (up/down), plus a
#   method-agreement diagnostic.
#
#   Only Jaccard matrices are combined (common [0,1] scale across all
#   three methods). Each method's Pearson/NES correlation matrix is on a
#   different, not-directly-comparable [-1,1] scale and is deliberately
#   NOT blended in here -- it stays inspectable individually in its own
#   output directory.
#
#   GSVA only scores KEGG gene sets (no GO), so the KEGG consensus can
#   combine up to 3 methods (GSVA + GSEA + ORA) but the GO consensus can
#   only combine 2 (GSEA + ORA). This is handled generically by checking
#   file.exists() per candidate source, not hardcoded -- a method's Jaccard
#   file missing (not run, or genuinely absent for GO) just lowers that
#   cell's method count rather than erroring.
#
#   Defaults to the primary stat-ranked "GSEA" directory (not "GSEA_LFC")
#   as its GSEA source, via the gsea_dir argument, to keep the consensus
#   well-defined; an LFC-based consensus variant could be added later as a
#   second, separately-named run without changing this design.
#
################################################################################

#===============================================================================
# Discover Which Methods Have a Jaccard File for This Family/Direction
#===============================================================================

discover_consensus_jaccard_sources <- function(config, pathway_family, direction, gsea_dir = "GSEA") {

    candidates <- list(
        GSVA = if (pathway_family == "KEGG") {
            file.path(config$output_directory, "GSVA", "Tables", paste0("compound_gsva_jaccard_", direction, ".csv"))
        } else {
            NA_character_
        },
        GSEA = file.path(
            config$output_directory, gsea_dir, pathway_family, "Tables",
            paste0("compound_gsea_", pathway_family, "_jaccard_", direction, ".csv")
        ),
        ORA = file.path(
            config$output_directory, pathway_family, "Tables",
            paste0("compound_", pathway_family, "_", direction, "_enrichment_jaccard.csv")
        )
    )

    candidates <- candidates[!is.na(candidates)]

    candidates[vapply(candidates, file.exists, logical(1))]

}

#===============================================================================
# Combine Per-Method Jaccard Matrices into a Consensus + Agreement Matrix
#===============================================================================

compute_consensus_matrices <- function(method_matrices, min_methods = 2) {

    common <- Reduce(intersect, lapply(method_matrices, rownames))

    method_matrices <- lapply(method_matrices, function(m) m[common, common, drop = FALSE])

    n <- length(common)
    k <- length(method_matrices)

    consensus <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(common, common))
    agreement <- consensus
    n_methods <- matrix(0L, nrow = n, ncol = n, dimnames = list(common, common))

    for (i in common) {
        for (j in common) {

            values <- vapply(method_matrices, function(m) m[i, j], numeric(1))
            values <- values[!is.na(values)]

            n_methods[i, j] <- length(values)

            if (length(values) >= min_methods) {
                consensus[i, j] <- mean(values)
                agreement[i, j] <- stats::sd(values)
            }

        }
    }

    list(consensus = consensus, agreement = agreement, n_methods = n_methods)

}

#===============================================================================
# Run Cross-Method Pathway Consensus Correlation for One Family (GO/KEGG)
#===============================================================================

run_pathway_family_consensus_correlation <- function(config, pathway_family, gsea_dir = "GSEA") {

    consensus_root <- if (gsea_dir == "GSEA_LFC") "Pathway_Consensus_LFC" else "Pathway_Consensus"

    table_dir <- file.path(config$output_directory, consensus_root, pathway_family, "Tables")
    figure_dir <- file.path(config$output_directory, consensus_root, pathway_family, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    methods_used_rows <- list()
    direction_results <- list()

    for (direction in c("up", "down")) {

        sources <- discover_consensus_jaccard_sources(config, pathway_family, direction, gsea_dir)

        if (length(sources) == 0) {
            warning(
                "No method Jaccard files found for ", pathway_family, " (", direction,
                "); skipping this direction's consensus."
            )
            next
        }

        message(sprintf(
            "Consensus for %s (%s): combining %d method(s) -- %s",
            pathway_family, direction, length(sources), paste(names(sources), collapse = ", ")
        ))

        method_matrices <- lapply(sources, function(f) as.matrix(utils::read.csv(f, row.names = 1, check.names = FALSE)))

        result <- compute_consensus_matrices(method_matrices, min_methods = config$pathway_consensus_min_methods)

        save_csv(result$consensus, file.path(table_dir, paste0("compound_pathway_consensus_", pathway_family, "_", direction, "_jaccard.csv")))
        save_csv(result$agreement, file.path(table_dir, paste0("compound_pathway_consensus_", pathway_family, "_", direction, "_agreement_sd.csv")))
        save_csv(result$n_methods, file.path(table_dir, paste0("compound_pathway_consensus_", pathway_family, "_", direction, "_n_methods.csv")))

        methods_used_rows[[length(methods_used_rows) + 1]] <- data.frame(
            direction = direction,
            methods = paste(names(sources), collapse = ", "),
            n_compounds = nrow(result$consensus),
            stringsAsFactors = FALSE
        )

        heatmap_height <- max(config$figure_height, nrow(result$consensus) * 0.5)
        heatmap_width <- heatmap_height + 2

        can_cluster_consensus <- !anyNA(result$consensus)
        can_cluster_agreement <- !anyNA(result$agreement)

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
                        paste0("Compound_Pathway_Consensus_", pathway_family, "_", direction, "_", filename_suffix, ".", ext)
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

        direction_label <- tools::toTitleCase(direction)

        save_matrix_heatmap(
            result$consensus, can_cluster_consensus,
            paste0("Compound Similarity: ", pathway_family, " ", direction_label, " Consensus Jaccard (", paste(names(sources), collapse = "+"), ")"),
            "jaccard_heatmap"
        )

        save_matrix_heatmap(
            result$agreement, can_cluster_agreement,
            paste0("Method Disagreement: ", pathway_family, " ", direction_label, " Jaccard SD Across Methods"),
            "agreement_heatmap"
        )

        direction_results[[direction]] <- result

    }

    methods_used_table <- if (length(methods_used_rows) > 0) do.call(rbind, methods_used_rows) else NULL

    if (!is.null(methods_used_table)) {
        save_csv(
            methods_used_table,
            file.path(table_dir, paste0("compound_pathway_consensus_", pathway_family, "_methods_used.csv")),
            row.names = FALSE
        )
    }

    list(up = direction_results[["up"]], down = direction_results[["down"]], methods_used = methods_used_table)

}

#===============================================================================
# Run Cross-Method Pathway Consensus Correlation (GO + KEGG)
#===============================================================================

run_pathway_consensus_correlation <- function(config, gsea_dir = "GSEA") {

    message("Running cross-method pathway-level compound-similarity consensus...")

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    go_result <- if (isTRUE(config$run_go)) {
        run_pathway_family_consensus_correlation(config, "GO", gsea_dir = gsea_dir)
    } else {
        NULL
    }

    kegg_result <- if (isTRUE(config$run_kegg)) {
        run_pathway_family_consensus_correlation(config, "KEGG", gsea_dir = gsea_dir)
    } else {
        NULL
    }

    message("Cross-method pathway-level compound-similarity consensus complete.")

    list(GO = go_result, KEGG = kegg_result)

}
