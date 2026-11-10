################################################################################
# RNAseqPipeline
#
# File: ercc_qc.R
#
# Description:
#   ERCC spike-in technical reproducibility QC: percent coefficient of
#   variation (%CV) of depth-normalized ERCC counts, per plate/library.
#   ERCC controls are added in a fixed amount per well, so their recovered
#   counts (normalized for each sample's sequencing depth) are a standard
#   signal of technical consistency, independent of biology.
#
################################################################################

#===============================================================================
# Run ERCC Spike-In %CV by Plate
#===============================================================================

run_ercc_cv <- function(counts, metadata, config) {

    message("Running ERCC spike-in %CV by plate...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    if (!requireNamespace("scales", quietly = TRUE)) {
        stop("Package 'scales' is required.")
    }

    counts <- counts[, rownames(metadata), drop = FALSE]

    ercc_mask <- startsWith(rownames(counts), config$ercc_prefix)

    if (!any(ercc_mask)) {
        stop("No rows match ercc_prefix: ", config$ercc_prefix)
    }

    total_counts <- colSums(counts)
    ercc_counts <- colSums(counts[ercc_mask, , drop = FALSE])

    sample_data <- data.frame(
        sample = colnames(counts),
        total_ercc_counts = as.integer(ercc_counts),
        total_counts = as.integer(total_counts),
        ercc_percent = ercc_counts / total_counts * 100
    )

    sample_data <- merge(sample_data, metadata, by = "sample")

    if (!("library" %in% names(sample_data))) {
        sample_data$library <- config$project_name
    }

    cv_summary <- stats::aggregate(
        ercc_percent ~ library,
        data = sample_data,
        FUN = function(x) c(mean = mean(x), sd = stats::sd(x), n = length(x))
    )

    cv_summary <- do.call(data.frame, cv_summary)

    names(cv_summary) <- c("library", "mean_ercc_percent", "sd_ercc_percent", "n_samples")

    cv_summary$cv_percent <- cv_summary$sd_ercc_percent / cv_summary$mean_ercc_percent * 100

    output_dir <- file.path(config$output_directory, "QC", "ERCC")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        sample_data,
        file.path(table_dir, "ercc_sample_data.csv"),
        row.names = FALSE
    )

    save_csv(
        cv_summary,
        file.path(table_dir, "ercc_cv_by_plate.csv"),
        row.names = FALSE
    )

    cv_plot <- ggplot2::ggplot(
        cv_summary,
        ggplot2::aes(x = library, y = cv_percent)
    ) +
        ggplot2::geom_col(fill = "#447899") +
        ggplot2::labs(
            title = "ERCC Spike-In %CV by Plate",
            x = "Library / Plate",
            y = "%CV of Depth-Normalized ERCC Counts"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    percent_plot <- ggplot2::ggplot(
        sample_data,
        ggplot2::aes(x = library, y = ercc_percent)
    ) +
        ggplot2::geom_boxplot(outlier.shape = NA, color = "#447899", fill = NA) +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.5, color = "#447899", size = 1.5) +
        ggplot2::labs(
            title = "Per-Sample ERCC Percent by Plate",
            x = "Library / Plate",
            y = "ERCC Counts (% of Total Counts)"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    depth_plot <- ggplot2::ggplot(
        sample_data,
        ggplot2::aes(x = total_counts, y = total_ercc_counts, color = library)
    ) +
        ggplot2::geom_point(alpha = 0.7, size = 2) +
        ggplot2::scale_x_continuous(labels = scales::comma) +
        ggplot2::scale_y_continuous(labels = scales::comma) +
        ggplot2::labs(
            title = "Per-Sample ERCC Count vs. Read Depth",
            x = "Total Counts (read depth)",
            y = "Total ERCC Counts",
            color = "Library / Plate"
        ) +
        ggplot2::theme_bw()

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "ERCC_CV_by_plate.png"),
            plot = cv_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

        ggplot2::ggsave(
            file.path(figure_dir, "ERCC_percent_by_plate.png"),
            plot = percent_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

        ggplot2::ggsave(
            file.path(figure_dir, "ERCC_vs_depth.png"),
            plot = depth_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "ERCC_CV_by_plate.pdf"),
            plot = cv_plot,
            width = config$figure_width,
            height = config$figure_height
        )

        ggplot2::ggsave(
            file.path(figure_dir, "ERCC_percent_by_plate.pdf"),
            plot = percent_plot,
            width = config$figure_width,
            height = config$figure_height
        )

        ggplot2::ggsave(
            file.path(figure_dir, "ERCC_vs_depth.pdf"),
            plot = depth_plot,
            width = config$figure_width,
            height = config$figure_height
        )

    }

    message("ERCC spike-in %CV by plate complete.")

    return(
        list(
            sample_data = sample_data,
            cv_summary = cv_summary,
            cv_plot = cv_plot,
            percent_plot = percent_plot,
            depth_plot = depth_plot
        )
    )

}
