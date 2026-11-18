################################################################################
# RNAseqPipeline
#
# File: gsea_compound_correlation.R
#
# Description:
#   Compound x compound similarity at the pathway level, read from the
#   per-contrast GSEA GO/KEGG results already written by run_gsea()
#   (R/gsea.R). Mirrors R/gsva_compound_correlation.R's design (Pearson
#   correlation of a per-pathway score across compounds), swapped onto
#   GSEA's NES column instead of GSVA limma's diff_score. Run separately
#   for GO and KEGG -- different pathway universes (GO term IDs vs. KEGG
#   pathway IDs), so a single combined matrix wouldn't make sense -- and,
#   same convention as run_gsva_compound_correlation(), using each
#   compound's TOP TESTED DOSE only (parse_dose_group()).
#
#   Unlike run_gsva_compound_correlation(), this does not also compute the
#   Jaccard up/down overlap or the paired Wilcoxon test -- those weren't
#   asked for here, just the correlation.
#
################################################################################

#===============================================================================
# Run Compound-vs-Compound GSEA Correlation for One Pathway Family (GO/KEGG)
#===============================================================================

run_gsea_pathway_family_correlation <- function(config, pathway_family) {

    table_dir <- file.path(config$output_directory, "GSEA", pathway_family, "Tables")
    figure_dir <- file.path(config$output_directory, "GSEA", pathway_family, "Figures")

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
        "Correlating %d compounds' GSEA %s NES profiles at their highest tested dose.",
        length(compounds), pathway_family
    ))

    #===========================================================================
    # Read Each Compound's NES Profile, Build the Pathway x Compound Matrix
    #===========================================================================

    gsea_file_for <- function(group) {
        file.path(table_dir, paste0("GSEA_", pathway_family, "_", group, ".csv"))
    }

    per_compound_nes <- lapply(compounds, function(cmpd) {
        group <- top_dose_table$group[top_dose_table$compound == cmpd]
        results <- utils::read.csv(gsea_file_for(group))
        stats::setNames(results$NES, results$ID)
    })

    names(per_compound_nes) <- compounds

    nes_matrix <- Reduce(
        function(a, b) merge(a, b, by = "ID", all = TRUE),
        lapply(compounds, function(cmpd) {
            df <- data.frame(ID = names(per_compound_nes[[cmpd]]), stringsAsFactors = FALSE)
            df[[cmpd]] <- as.numeric(per_compound_nes[[cmpd]])
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
    # Save Table and Heatmap
    #===========================================================================

    save_csv(
        correlation_matrix,
        file.path(table_dir, paste0("compound_gsea_", pathway_family, "_nes_correlation.csv"))
    )

    heatmap_height <- max(config$figure_height, length(compounds) * 0.5)
    heatmap_width <- heatmap_height + 2

    if (isTRUE(config$save_png)) {

        pheatmap::pheatmap(
            correlation_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = paste0("Compound Similarity: GSEA ", pathway_family, " NES Correlation (Top Dose)"),
            filename = file.path(figure_dir, paste0("Compound_GSEA_", pathway_family, "_nes_correlation_heatmap.png")),
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
            main = paste0("Compound Similarity: GSEA ", pathway_family, " NES Correlation (Top Dose)"),
            filename = file.path(figure_dir, paste0("Compound_GSEA_", pathway_family, "_nes_correlation_heatmap.pdf")),
            width = heatmap_width,
            height = heatmap_height
        )

    }

    list(top_dose_table = top_dose_table, correlation_matrix = correlation_matrix)

}

#===============================================================================
# Run Compound-vs-Compound GSEA Pathway Correlation (GO + KEGG)
#===============================================================================

run_gsea_compound_correlation <- function(config) {

    message("Running compound-vs-compound GSEA pathway correlation...")

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    go_result <- run_gsea_pathway_family_correlation(config, "GO")
    kegg_result <- run_gsea_pathway_family_correlation(config, "KEGG")

    message("Compound-vs-compound GSEA pathway correlation complete.")

    list(GO = go_result, KEGG = kegg_result)

}
