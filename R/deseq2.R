################################################################################
# RNAseqPipeline
#
# File: deseq2.R
#
# Description:
#   Minimal working DESeq2 module.
#
################################################################################

#===============================================================================
# Run DESeq2 Differential Expression Analysis
#===============================================================================

run_deseq2 <- function(counts, metadata, config) {

    message("Running DESeq2...")

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.")
    }

    if (!requireNamespace("apeglm", quietly = TRUE) &&
        isTRUE(config$lfc_shrinkage)) {
        warning("Package 'apeglm' not found. LFC shrinkage will be skipped.")
        config$lfc_shrinkage <- FALSE
    }

    metadata <- relevel_group(metadata, config$reference_group)

    counts <- counts[, rownames(metadata), drop = FALSE]

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(as.matrix(counts)),
        colData = metadata,
        design = stats::as.formula(config$design_formula)
    )

    dds <- DESeq2::DESeq(dds)

    output_dir <- file.path(config$output_directory, "DESeq2")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")
    object_dir <- file.path(output_dir, "Objects")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

    normalized_counts <- DESeq2::counts(dds, normalized = TRUE)

    write.csv(
        normalized_counts,
        file.path(table_dir, "DESeq2_normalized_counts.csv"),
        quote = FALSE
    )

    saveRDS(
        dds,
        file.path(object_dir, "dds_object.rds")
    )

    result_names <- DESeq2::resultsNames(dds)

    contrast_names <- result_names[
        grepl("^group_.*_vs_", result_names)
    ]

    results_list <- list()
    shrunk_results_list <- list()

    for (coef_name in contrast_names) {

        message("Processing DESeq2 contrast: ", coef_name)

        res_raw <- DESeq2::results(
            dds,
            name = coef_name,
            alpha = config$alpha
        )

        res_raw_df <- as.data.frame(res_raw)
        res_raw_df$gene <- rownames(res_raw_df)

        res_raw_df <- res_raw_df[
            ,
            c(
                "gene",
                setdiff(colnames(res_raw_df), "gene")
            )
        ]

        write.csv(
            res_raw_df,
            file.path(
                table_dir,
                paste0("DESeq2_", coef_name, "_raw.csv")
            ),
            row.names = FALSE,
            quote = FALSE
        )

        results_list[[coef_name]] <- res_raw_df

        if (isTRUE(config$lfc_shrinkage)) {

            res_shrunk <- DESeq2::lfcShrink(
                dds,
                coef = coef_name,
                type = config$shrinkage_method
            )

            res_shrunk_df <- as.data.frame(res_shrunk)
            res_shrunk_df$gene <- rownames(res_shrunk_df)

            res_shrunk_df <- res_shrunk_df[
                ,
                c(
                    "gene",
                    setdiff(colnames(res_shrunk_df), "gene")
                )
            ]

            write.csv(
                res_shrunk_df,
                file.path(
                    table_dir,
                    paste0("DESeq2_", coef_name, "_shrunk.csv")
                ),
                row.names = FALSE,
                quote = FALSE
            )

            shrunk_results_list[[coef_name]] <- res_shrunk_df

            save_deseq2_ma_plot(
                res_raw,
                res_shrunk,
                coef_name,
                figure_dir,
                config
            )

        } else {

            save_deseq2_ma_plot(
                res_raw,
                NULL,
                coef_name,
                figure_dir,
                config
            )

        }
    }

    saveRDS(
        results_list,
        file.path(object_dir, "deseq2_raw_results_list.rds")
    )

    saveRDS(
        shrunk_results_list,
        file.path(object_dir, "deseq2_shrunk_results_list.rds")
    )

    message("DESeq2 analysis complete.")

    return(
        list(
            dds = dds,
            raw_results = results_list,
            shrunk_results = shrunk_results_list,
            normalized_counts = normalized_counts
        )
    )
}


#===============================================================================
# Save DESeq2 MA Plot
#===============================================================================

save_deseq2_ma_plot <- function(
    res_raw,
    res_shrunk = NULL,
    contrast_name,
    figure_dir,
    config
) {

    safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", contrast_name)

    if (isTRUE(config$save_png)) {

        png(
            filename = file.path(
                figure_dir,
                paste0("MA_plot_", safe_name, ".png")
            ),
            width = config$figure_width,
            height = config$figure_height,
            units = "in",
            res = config$figure_dpi
        )

        if (!is.null(res_shrunk)) {

            par(mfrow = c(1, 2))

            DESeq2::plotMA(
                res_raw,
                main = paste("Raw:", contrast_name)
            )

            DESeq2::plotMA(
                res_shrunk,
                main = paste("Shrunk:", contrast_name)
            )

        } else {

            DESeq2::plotMA(
                res_raw,
                main = contrast_name
            )

        }

        dev.off()
    }

    if (isTRUE(config$save_pdf)) {

        pdf(
            file = file.path(
                figure_dir,
                paste0("MA_plot_", safe_name, ".pdf")
            ),
            width = config$figure_width,
            height = config$figure_height
        )

        if (!is.null(res_shrunk)) {

            par(mfrow = c(1, 2))

            DESeq2::plotMA(
                res_raw,
                main = paste("Raw:", contrast_name)
            )

            DESeq2::plotMA(
                res_shrunk,
                main = paste("Shrunk:", contrast_name)
            )

        } else {

            DESeq2::plotMA(
                res_raw,
                main = contrast_name
            )

        }

        dev.off()
    }

    invisible(TRUE)
}
