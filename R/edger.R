################################################################################
# RNAseqPipeline
#
# File: edger.R
#
# Description:
#   Minimal working edgeR module.
#
################################################################################

#===============================================================================
# Run edgeR Differential Expression Analysis
#===============================================================================

run_edger <- function(counts, metadata, config) {

    message("Running edgeR...")

    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package 'edgeR' is required.")
    }

    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required for edgeR MD plots.")
    }

    metadata <- relevel_group(metadata, config$reference_group)

    counts <- counts[, rownames(metadata), drop = FALSE]

    design <- stats::model.matrix(
        stats::as.formula(config$design_formula),
        data = metadata
    )

    dge <- edgeR::DGEList(counts = round(as.matrix(counts)))
    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(dge, design)

    fit <- edgeR::glmQLFit(dge, design)

    output_dir <- file.path(config$output_directory, "edgeR")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")
    object_dir <- file.path(output_dir, "Objects")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

    normalized_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

    write.csv(
        normalized_counts,
        file.path(table_dir, "edgeR_normalized_counts.csv"),
        quote = FALSE
    )

    saveRDS(
        dge,
        file.path(object_dir, "dge_object.rds")
    )

    saveRDS(
        fit,
        file.path(object_dir, "glmQLFit_object.rds")
    )

    coef_names <- colnames(design)

    contrast_names <- coef_names[
        grepl("^group", coef_names)
    ]

    results_list <- list()

    for (coef_name in contrast_names) {

        contrast_label <- sub("^group", "", coef_name)

        message("Processing edgeR contrast: ", contrast_label)

        qlf <- edgeR::glmQLFTest(fit, coef = coef_name)

        res_df <- edgeR::topTags(
            qlf,
            n = Inf,
            sort.by = "none"
        )$table

        res_df$gene <- rownames(res_df)

        res_df <- res_df[
            ,
            c(
                "gene",
                setdiff(colnames(res_df), "gene")
            )
        ]

        write.csv(
            res_df,
            file.path(
                table_dir,
                paste0(
                    "edgeR_", contrast_label,
                    "_vs_", config$reference_group,
                    "_raw.csv"
                )
            ),
            row.names = FALSE,
            quote = FALSE
        )

        results_list[[coef_name]] <- res_df

        save_edger_md_plot(
            qlf,
            contrast_label,
            figure_dir,
            config
        )

    }

    saveRDS(
        results_list,
        file.path(object_dir, "edger_results_list.rds")
    )

    message("edgeR analysis complete.")

    return(
        list(
            dge = dge,
            fit = fit,
            results = results_list,
            normalized_counts = normalized_counts
        )
    )
}


#===============================================================================
# Save edgeR MD Plot
#===============================================================================

save_edger_md_plot <- function(
    qlf,
    contrast_name,
    figure_dir,
    config
) {

    safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", contrast_name)

    if (isTRUE(config$save_png)) {

        png(
            filename = file.path(
                figure_dir,
                paste0("MD_plot_", safe_name, ".png")
            ),
            width = config$figure_width,
            height = config$figure_height,
            units = "in",
            res = config$figure_dpi
        )

        limma::plotMD(
            qlf,
            main = contrast_name
        )

        dev.off()
    }

    if (isTRUE(config$save_pdf)) {

        pdf(
            file = file.path(
                figure_dir,
                paste0("MD_plot_", safe_name, ".pdf")
            ),
            width = config$figure_width,
            height = config$figure_height
        )

        limma::plotMD(
            qlf,
            main = contrast_name
        )

        dev.off()
    }

    invisible(TRUE)
}
