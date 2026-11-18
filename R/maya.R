################################################################################
# RNAseqPipeline
#
# File: maya.R
#
# Description:
#   MAYA (Multi-modal quantification of pathway Activity; Landais & Vallot,
#   Nature Communications 2023, github.com/One-Biosciences/MAYA)
#   pathway-activity scoring on KEGG gene sets, reusing the Ensembl-ID gene
#   set builder from R/gsva.R (build_kegg_gene_sets()) so MAYA scores the
#   same pathway family GSVA does.
#
#   *** EXPERIMENTAL / METHODOLOGICAL CAVEAT ***
#   MAYA's per-pathway-PCA + igraph/leidenbase "activation mode" detection
#   was designed and validated on single-cell RNA-seq (thousands of cells,
#   sparse/dropout noise). AMP0961 is bulk RNA-seq (~480-600 samples after
#   gene filtering / T007 exclusion) -- a different scale and noise regime
#   than MAYA was built for. Its output here is an exploratory cross-check
#   against GSVA, not a validated primary pathway-activity measure. The
#   optional differential-activity test (run_maya_limma(), gated on
#   config$run_maya_differential) was left off by default until the
#   activity-matrix heatmap was visually sanity-checked against GSVA's full
#   KEGG heatmap on the real AMP0961 data (both show the same kind of
#   large-block, non-degenerate structure); config_AMP0961.R now enables it,
#   but its p-values should still be read as exploratory, not with the same
#   confidence as GSVA's limma results.
#
################################################################################

#===============================================================================
# Save MAYA Heatmap of Top-Variable Activation Modes
#===============================================================================
#
# Mirrors save_gsva_heatmap() in R/gsva.R, retargeted at MAYA's "modes"
# (a pathway can decompose into more than one activation mode) instead of
# one row per pathway.
#===============================================================================

save_maya_heatmap <- function(activity_matrix, metadata, config, figure_dir) {

    mode_vars <- apply(activity_matrix, 1, stats::var)

    n_heatmap <- min(config$maya_heatmap_top_modes, nrow(activity_matrix))

    heatmap_modes <- names(sort(mode_vars, decreasing = TRUE))[seq_len(n_heatmap)]

    heatmap_matrix <- activity_matrix[heatmap_modes, , drop = FALSE]

    annotation_col <- data.frame(
        Group = metadata[colnames(heatmap_matrix), "group"],
        row.names = colnames(heatmap_matrix)
    )

    heatmap_height <- max(config$figure_height, 3 + n_heatmap * 0.2)
    heatmap_width <- max(config$figure_width, ncol(heatmap_matrix) * 0.1) + 4

    if (isTRUE(config$save_png)) {

        pheatmap::pheatmap(
            heatmap_matrix,
            scale = "row",
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            annotation_col = annotation_col,
            show_colnames = FALSE,
            main = "Top Variable KEGG Pathway Activation Modes (MAYA)",
            filename = file.path(figure_dir, "Heatmap_MAYA_KEGG_top_modes.png"),
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
            main = "Top Variable KEGG Pathway Activation Modes (MAYA)",
            filename = file.path(figure_dir, "Heatmap_MAYA_KEGG_top_modes.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

    }

    invisible(TRUE)

}

#===============================================================================
# Run limma on MAYA Activity Scores (Differential Pathway Activity per Contrast)
#===============================================================================
#
# Mirrors run_gsva_limma() in R/gsva.R (same design-formula/library-
# adjustment/relevel logic), operating on the MAYA activity matrix instead
# of GSVA scores. No volcano-plot labeling machinery -- that's tied to
# run_gsva_pairwise_limma(), which MAYA has no analog for here.
#===============================================================================

run_maya_limma <- function(activity_matrix, metadata, config, table_dir) {

    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required.")
    }

    metadata <- relevel_group(metadata, config$reference_group)

    activity_matrix <- activity_matrix[, rownames(metadata), drop = FALSE]

    design_formula <- adjust_design_formula_for_library(config$design_formula, metadata)

    design <- stats::model.matrix(
        stats::as.formula(design_formula),
        data = metadata
    )

    fit <- limma::lmFit(activity_matrix, design)
    fit <- limma::eBayes(fit)

    coef_names <- colnames(design)

    contrast_names <- coef_names[grepl("^group", coef_names)]

    results_list <- list()

    for (coef_name in contrast_names) {

        contrast_label <- sub("^group", "", coef_name)

        res_df <- limma::topTable(
            fit,
            coef = coef_name,
            number = Inf,
            sort.by = "P"
        )

        res_df$mode <- rownames(res_df)

        res_df <- res_df[
            ,
            c("mode", setdiff(colnames(res_df), "mode"))
        ]

        colnames(res_df)[colnames(res_df) == "logFC"] <- "diff_score"

        save_csv(
            res_df,
            file.path(
                table_dir,
                paste0(
                    "MAYA_KEGG_", contrast_label,
                    "_vs_", config$reference_group, "_limma.csv"
                )
            ),
            row.names = FALSE
        )

        results_list[[contrast_label]] <- res_df

    }

    results_list

}

#===============================================================================
# Run MAYA
#===============================================================================

run_maya <- function(counts, metadata, config, exclude_pattern = config$exclude_groups_pattern) {

    message("Running MAYA (KEGG gene sets, EXPERIMENTAL on bulk data -- see file header)...")

    if (!requireNamespace("MAYA", quietly = TRUE)) {
        stop(
            "Package 'MAYA' is required but not installed. Install with:\n",
            '  remotes::install_github("one-biosciences/maya")'
        )
    }

    org_pkg <- get_organism_db(config$organism)

    kegg_organism <- get_kegg_organism_code(config$organism)

    expr_mat <- as.matrix(counts)

    expr_mat <- expr_mat[, rownames(metadata), drop = FALSE]

    gene_sets <- build_kegg_gene_sets(
        kegg_organism = kegg_organism, org_pkg = org_pkg, universe_gene_ids = rownames(expr_mat)
    )

    # Match GSVA's size filtering over the same KEGG/Ensembl gene-set family
    # (R/gsva.R passes these same bounds to GSVA::gsvaParam) so MAYA scores
    # a comparable set of pathways rather than a differently-filtered one.
    set_sizes <- lengths(gene_sets)

    gene_sets <- gene_sets[
        set_sizes >= config$gsva_min_gene_set_size & set_sizes <= config$gsva_max_gene_set_size
    ]

    message("Scoring ", length(gene_sets), " candidate KEGG pathways across ", ncol(expr_mat), " samples...")

    maya_result <- MAYA::MAYA_pathway_analysis(
        expr_mat = expr_mat,
        modules_list = gene_sets,
        is_logcpm = FALSE,
        compute_umap = FALSE
    )

    activity_matrix <- maya_result$activity_matrix

    output_dir <- file.path(config$output_directory, "MAYA")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")
    object_dir <- file.path(output_dir, "Objects")

    for (d in c(table_dir, figure_dir, object_dir)) {
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }

    save_csv(
        data.frame(mode = rownames(activity_matrix), as.data.frame(as.matrix(activity_matrix))),
        file.path(table_dir, "MAYA_KEGG_activity_matrix.csv"),
        row.names = FALSE
    )

    saveRDS(maya_result, file.path(object_dir, "maya_result.rds"))
    saveRDS(activity_matrix, file.path(object_dir, "maya_activity_matrix.rds"))

    # Match GSVA/enrichment: drop T007-pattern co-treatment groups from the
    # heatmap and differential test, but score all samples above.
    heatmap_metadata <- metadata
    heatmap_matrix <- activity_matrix

    if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {

        excluded_groups <- unique(metadata$group[grepl(exclude_pattern, metadata$group)])

        if (length(excluded_groups) > 0) {

            message("Excluding groups matching '", exclude_pattern, "': ", paste(excluded_groups, collapse = ", "))

            heatmap_metadata <- metadata[!metadata$group %in% excluded_groups, , drop = FALSE]
            heatmap_matrix <- activity_matrix[, rownames(heatmap_metadata), drop = FALSE]

        }

    }

    save_maya_heatmap(heatmap_matrix, heatmap_metadata, config, figure_dir)

    diff_results <- NULL

    if (isTRUE(config$run_maya_differential)) {

        diff_results <- run_maya_limma(heatmap_matrix, heatmap_metadata, config, table_dir)

    }

    message("MAYA complete.")

    return(
        list(
            activity_matrix = activity_matrix,
            gene_sets = gene_sets,
            maya_object = maya_result,
            diff_results = diff_results
        )
    )

}
