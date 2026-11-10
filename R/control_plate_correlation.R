################################################################################
# RNAseqPipeline
#
# File: control_plate_correlation.R
#
# Description:
#   Inter-plate reproducibility QC: Pearson correlation of Control-well
#   expression profiles (log2 DESeq2-normalized counts), comparing each
#   plate/library's Controls against every other plate's. A plate whose
#   Controls correlate poorly with the rest likely carries a batch effect
#   distinct from real biology.
#
################################################################################

#===============================================================================
# Run Inter-Plate Control-Well Correlation
#===============================================================================

run_control_plate_correlation <- function(normalized_counts, metadata, config) {

    message("Running inter-plate Control-well correlation...")

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    control_metadata <- metadata[metadata$group == config$reference_group, , drop = FALSE]

    if (nrow(control_metadata) == 0) {
        stop("No samples found for reference_group: ", config$reference_group)
    }

    if (!("library" %in% names(control_metadata)) ||
            length(unique(control_metadata$library)) < 2) {

        message(
            "Only one plate/library present; nothing to inter-plate compare. Skipping."
        )

        return(invisible(NULL))

    }

    log_counts <- log2(normalized_counts[, control_metadata$sample, drop = FALSE] + 1)

    sample_cor <- stats::cor(log_counts, method = "pearson")

    libraries <- sort(unique(control_metadata$library))

    plate_profiles <- vapply(
        libraries,
        function(lib) {
            samples_in_plate <- control_metadata$sample[control_metadata$library == lib]
            rowMeans(log_counts[, samples_in_plate, drop = FALSE])
        },
        numeric(nrow(log_counts))
    )

    colnames(plate_profiles) <- libraries

    plate_cor <- stats::cor(plate_profiles, method = "pearson")

    output_dir <- file.path(config$output_directory, "QC", "ControlCorrelation")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(plate_cor, file.path(table_dir, "plate_correlation_matrix.csv"))

    save_csv(sample_cor, file.path(table_dir, "sample_correlation_matrix.csv"))

    plate_height <- max(config$figure_height, length(libraries) * 0.5)
    plate_width <- plate_height + 2

    sample_annotation <- data.frame(
        Library = control_metadata$library,
        row.names = control_metadata$sample
    )

    sample_height <- max(config$figure_height, nrow(control_metadata) * 0.2)
    sample_width <- sample_height + 2

    if (isTRUE(config$save_png)) {

        pheatmap::pheatmap(
            plate_cor,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            main = "Inter-Plate Control Correlation (log2 normalized counts)",
            filename = file.path(figure_dir, "Plate_correlation_heatmap.png"),
            width = plate_width,
            height = plate_height,
            res = config$figure_dpi
        )

        pheatmap::pheatmap(
            sample_cor,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            annotation_col = sample_annotation,
            show_colnames = FALSE,
            show_rownames = FALSE,
            main = "Control-Sample Correlation (log2 normalized counts)",
            filename = file.path(figure_dir, "Sample_correlation_heatmap.png"),
            width = sample_width,
            height = sample_height,
            res = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        pheatmap::pheatmap(
            plate_cor,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            main = "Inter-Plate Control Correlation (log2 normalized counts)",
            filename = file.path(figure_dir, "Plate_correlation_heatmap.pdf"),
            width = plate_width,
            height = plate_height
        )

        pheatmap::pheatmap(
            sample_cor,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            annotation_col = sample_annotation,
            show_colnames = FALSE,
            show_rownames = FALSE,
            main = "Control-Sample Correlation (log2 normalized counts)",
            filename = file.path(figure_dir, "Sample_correlation_heatmap.pdf"),
            width = sample_width,
            height = sample_height
        )

    }

    message("Inter-plate Control-well correlation complete.")

    return(
        list(
            plate_correlation = plate_cor,
            sample_correlation = sample_cor
        )
    )

}
