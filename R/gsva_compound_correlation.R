################################################################################
# RNAseqPipeline
#
# File: gsva_compound_correlation.R
#
# Description:
#   Compound x compound similarity at the PATHWAY level: Pearson correlation
#   of each compound's GSVA KEGG pathway diff_score profile (vs Control, at
#   its own top dose), read from the per-contrast limma results already
#   written by run_gsva() (R/gsva.R). Direct pathway-level analog of
#   R/compound_similarity.R's gene-level log2FoldChange correlation.
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
    # Build the Pathway x Compound diff_score Matrix
    #===========================================================================

    diff_score_matrix <- Reduce(
        function(a, b) merge(a, b, by = "pathway", all = TRUE),
        lapply(compounds, function(cmpd) {
            group <- eligible_table$group[eligible_table$compound == cmpd]
            results <- read.csv(gsva_file_for(group))
            df <- data.frame(pathway = results$pathway, stringsAsFactors = FALSE)
            df[[cmpd]] <- results$diff_score
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
    # Save Table and Heatmap
    #===========================================================================

    dir.create(gsva_table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(correlation_matrix, file.path(gsva_table_dir, "compound_gsva_pathway_correlation.csv"))

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

    }

    message("Compound-vs-compound GSVA pathway correlation complete.")

    return(
        list(
            top_dose_table = eligible_table,
            correlation_matrix = correlation_matrix
        )
    )

}
