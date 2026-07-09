################################################################################
# RNAseqPipeline
#
# File: limma.R
#
# Description:
#   Minimal working limma-voom module.
#
################################################################################

#===============================================================================
# Run limma-voom Differential Expression Analysis
#===============================================================================

run_limma <- function(counts, metadata, config) {

    message("Running limma-voom...")

    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required.")
    }

    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package 'edgeR' is required for limma-voom normalization.")
    }

    metadata <- relevel_group(metadata, config$reference_group)

    counts <- counts[, rownames(metadata), drop = FALSE]

    design <- stats::model.matrix(
        stats::as.formula(config$design_formula),
        data = metadata
    )

    dge <- edgeR::DGEList(counts = round(as.matrix(counts)))
    dge <- edgeR::calcNormFactors(dge)

    output_dir <- file.path(config$output_directory, "limma")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")
    object_dir <- file.path(output_dir, "Objects")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

    voom_fit <- limma::voom(
        dge,
        design,
        plot = FALSE
    )

    write.csv(
        voom_fit$E,
        file.path(table_dir, "limma_voom_normalized_counts.csv"),
        quote = FALSE
    )

    fit <- limma::lmFit(voom_fit, design)
    fit <- limma::eBayes(fit)

    saveRDS(
        voom_fit,
        file.path(object_dir, "voom_object.rds")
    )

    saveRDS(
        fit,
        file.path(object_dir, "eBayes_fit_object.rds")
    )

    coef_names <- colnames(design)

    contrast_names <- coef_names[
        grepl("^group", coef_names)
    ]

    results_list <- list()

    for (coef_name in contrast_names) {

        contrast_label <- sub("^group", "", coef_name)

        message("Processing limma contrast: ", contrast_label)

        res_df <- limma::topTable(
            fit,
            coef = coef_name,
            number = Inf,
            sort.by = "none"
        )

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
                    "limma_", contrast_label,
                    "_vs_", config$reference_group,
                    "_raw.csv"
                )
            ),
            row.names = FALSE,
            quote = FALSE
        )

        results_list[[coef_name]] <- res_df

        save_limma_md_plot(
            fit,
            coef_name,
            contrast_label,
            figure_dir,
            config
        )

    }

    saveRDS(
        results_list,
        file.path(object_dir, "limma_results_list.rds")
    )

    message("limma-voom analysis complete.")

    return(
        list(
            voom = voom_fit,
            fit = fit,
            results = results_list,
            normalized_counts = voom_fit$E
        )
    )
}


#===============================================================================
# Save limma MD Plot
#===============================================================================

save_limma_md_plot <- function(
    fit,
    coef_name,
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
            fit,
            coef = coef_name,
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
            fit,
            coef = coef_name,
            main = contrast_name
        )

        dev.off()
    }

    invisible(TRUE)
}
