################################################################################
# RNAseqPipeline
#
# File: ercc_dose_heatmap.R
#
# Description:
#   Heatmap of log10(nUMI/nERCC) by drug x dose. ERCC spike-ins are added in
#   a fixed amount per well, so a low ratio of endogenous transcript UMIs to
#   spike-in UMIs flags conditions where comparatively little real RNA was
#   recovered. Reuses parse_dose_group() (R/dose_response.R) to assign each
#   sample to a (compound, dose) cell, which naturally excludes Control (no
#   dose) and co-treatment groups (two compound+dose pairs).
#
################################################################################

#===============================================================================
# Format a Dose (in nM) for Display
#===============================================================================

format_dose_label <- function(dose_nM) {

    ifelse(
        dose_nM >= 1000,
        paste0(dose_nM / 1000, "uM"),
        paste0(dose_nM, "nM")
    )

}

#===============================================================================
# Run nUMI / nERCC Ratio Heatmap by Drug x Dose
#===============================================================================

run_ercc_umi_ratio_heatmap <- function(counts, metadata, config) {

    message("Running nUMI/nERCC ratio heatmap by drug x dose...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    counts <- counts[, rownames(metadata), drop = FALSE]

    endogenous_mask <- startsWith(rownames(counts), config$saturation_gene_prefix)
    ercc_mask <- startsWith(rownames(counts), config$ercc_prefix)

    n_umi <- colSums(counts[endogenous_mask, , drop = FALSE])
    n_ercc <- colSums(counts[ercc_mask, , drop = FALSE])

    sample_data <- data.frame(
        sample = colnames(counts),
        n_umi = as.integer(n_umi),
        n_ercc = as.integer(n_ercc),
        log_ratio = log10(n_umi / n_ercc)
    )

    sample_data <- merge(sample_data, metadata, by = "sample")

    parsed <- lapply(as.character(sample_data$group), parse_dose_group)

    keep <- !vapply(parsed, is.null, logical(1))

    message(sprintf(
        "%d of %d samples excluded (no single-compound dose parsed, e.g. Control or co-treatment).",
        sum(!keep), length(keep)
    ))

    filtered <- sample_data[keep, ]
    filtered$compound <- vapply(parsed[keep], `[[`, character(1), "compound")
    filtered$dose_nM <- vapply(parsed[keep], `[[`, numeric(1), "dose_nM")

    if (nrow(filtered) == 0) {
        stop("No samples with a single-compound dose were found; nothing to plot.")
    }

    cell_data <- stats::aggregate(
        cbind(log_ratio, n_umi, n_ercc) ~ compound + dose_nM,
        data = filtered,
        FUN = function(x) c(mean = mean(x), n = length(x))
    )

    cell_data <- do.call(data.frame, cell_data)

    cell_data <- cell_data[, c("compound", "dose_nM", "log_ratio.mean", "log_ratio.n")]

    names(cell_data) <- c("compound", "dose_nM", "mean_log_ratio", "n_samples")

    cell_data <- cell_data[order(cell_data$compound, cell_data$dose_nM), ]

    cell_data$dose_rank <- stats::ave(
        cell_data$dose_nM,
        cell_data$compound,
        FUN = function(x) rank(x, ties.method = "first")
    )

    cell_data$dose_rank <- factor(
        paste("Dose", cell_data$dose_rank),
        levels = paste("Dose", sort(unique(cell_data$dose_rank)))
    )

    cell_data$dose_label <- format_dose_label(cell_data$dose_nM)

    output_dir <- file.path(config$output_directory, "QC", "ERCC")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        cell_data,
        file.path(table_dir, "ercc_umi_ratio_by_dose.csv"),
        row.names = FALSE
    )

    heatmap_plot <- ggplot2::ggplot(
        cell_data,
        ggplot2::aes(x = dose_rank, y = compound, fill = mean_log_ratio)
    ) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = dose_label), size = 2.8) +
        ggplot2::scale_fill_viridis_c(name = "log10(nUMI/nERCC)") +
        ggplot2::labs(
            title = "nUMI/nERCC Ratio by Drug and Dose",
            x = "Dose (ordinal rank, lowest to highest per compound)",
            y = "Compound"
        ) +
        ggplot2::theme_bw()

    n_compounds <- length(unique(cell_data$compound))
    n_doses <- nlevels(cell_data$dose_rank)

    plot_width <- max(config$figure_width, n_doses * 1.2)
    plot_height <- max(config$figure_height, n_compounds * 0.4)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "ERCC_UMI_ratio_heatmap.png"),
            plot = heatmap_plot,
            width = plot_width,
            height = plot_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "ERCC_UMI_ratio_heatmap.pdf"),
            plot = heatmap_plot,
            width = plot_width,
            height = plot_height
        )

    }

    message("nUMI/nERCC ratio heatmap complete.")

    return(
        list(
            cell_data = cell_data,
            plot = heatmap_plot
        )
    )

}
