################################################################################
# RNAseqPipeline
#
# File: saturation.R
#
# Description:
#   Gene detection saturation (rarefaction) curves, per individual sample,
#   computed from non-deduplicated read counts against the true total
#   sequencing depth reported in the provider's per-sample QC table (not
#   just gene-assigned reads). UMI-deduplicated counts are not used here:
#   they discard the PCR/optical duplicate structure needed to ask "how
#   many genes would be detected at a given raw sequencing depth".
#
################################################################################

#===============================================================================
# Expected Number of Detected Genes at a Given Downsampled Read Depth
#===============================================================================

expected_detected_genes <- function(counts_vec, depth) {

    total <- sum(counts_vec)

    if (depth <= 0) {
        return(0)
    }

    if (depth >= total) {
        return(sum(counts_vec > 0))
    }

    nonzero <- counts_vec[counts_vec > 0]

    if (length(nonzero) == 0) {
        return(0)
    }

    # Group by unique count value (many genes share the same count,
    # especially at low counts) so lchoose() is evaluated once per
    # distinct value rather than once per gene.
    tbl <- table(nonzero)

    c_vals <- as.numeric(names(tbl))
    freq <- as.numeric(tbl)

    unseen <- total - c_vals

    log_p_absent <- ifelse(
        unseen < depth,
        -Inf,
        lchoose(unseen, depth) - lchoose(total, depth)
    )

    sum((1 - exp(log_p_absent)) * freq)

}

#===============================================================================
# Read Per-Sample Sequencing QC Table
#===============================================================================

read_saturation_qc <- function(qc_files) {

    tables <- lapply(qc_files, data.table::fread, data.table = FALSE)

    qc <- do.call(rbind, tables)

    qc <- qc[qc$sampleID != "Unknown", ]

    if (any(duplicated(qc$sampleID))) {
        stop("Duplicated sampleID values across saturation_qc_file inputs.")
    }

    rownames(qc) <- qc$sampleID

    qc

}

#===============================================================================
# Run Gene Detection Saturation Curve
#===============================================================================

run_saturation_curve <- function(metadata, config) {

    message("Running gene detection saturation curve...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    reads <- read_counts(config$saturation_reads_file)

    qc <- read_saturation_qc(config$saturation_qc_file)

    samples <- colnames(reads)

    missing_qc <- setdiff(samples, rownames(qc))

    if (length(missing_qc) > 0) {
        stop(
            "saturation_qc_file is missing QC rows for samples:\n",
            paste(head(missing_qc, 20), collapse = "\n")
        )
    }

    gene_mask <- startsWith(rownames(reads), config$saturation_gene_prefix)

    if (!any(gene_mask)) {
        stop("No genes match saturation_gene_prefix: ", config$saturation_gene_prefix)
    }

    step <- config$saturation_step

    curve_rows <- vector("list", length(samples))
    summary_rows <- vector("list", length(samples))

    names(curve_rows) <- samples
    names(summary_rows) <- samples

    for (sample in samples) {

        info <- qc[sample, ]

        total <- as.integer(info$depth)

        assigned <- as.integer(sum(reads[, sample]))

        if (assigned != as.integer(info$mapped_to_exons)) {
            stop(sprintf(
                "%s: read sum (%d) differs from mapped_to_exons (%d).",
                sample, assigned, as.integer(info$mapped_to_exons)
            ))
        }

        if (total != assigned + as.integer(info$no_feature) + as.integer(info$not_aligned)) {
            stop(sprintf("%s: QC depth categories do not reconcile.", sample))
        }

        c_vec <- reads[gene_mask, sample]

        milestones <- Filter(function(d) d <= total, c(250000, 500000, 750000, 1000000, 2000000))

        depths <- sort(unique(c(0, total, seq(step, total, by = step), milestones)))

        y <- vapply(depths, function(d) expected_detected_genes(c_vec, d), numeric(1))

        if (any(diff(y) < -1e-5) || abs(y[length(y)] - sum(c_vec > 0)) > 1e-5) {
            stop(sprintf("%s: rarefaction monotonicity or endpoint check failed.", sample))
        }

        curve_rows[[sample]] <- data.frame(
            sample = sample,
            depth = depths,
            detected_genes = y
        )

        at1 <- if (total >= 1e6) expected_detected_genes(c_vec, 1e6) else NA_real_
        at2 <- if (total >= 2e6) expected_detected_genes(c_vec, 2e6) else NA_real_

        summary_rows[[sample]] <- data.frame(
            sample = sample,
            total_reads = total,
            gene_assigned_reads_including_spikes = assigned,
            endogenous_assigned_reads = as.integer(sum(c_vec)),
            full_depth_endogenous_genes = as.integer(sum(c_vec > 0)),
            expected_genes_at_1M = at1,
            expected_genes_at_2M = at2,
            additional_genes_1M_to_2M = at2 - at1,
            benchmark_1M = if (is.na(at1)) {
                "not_evaluable"
            } else if (at1 >= 10000) {
                "at_or_above_10000"
            } else {
                "below_10000"
            },
            provider_detected_genes = info$nb_of_detected_genes
        )

    }

    curve_data <- do.call(rbind, curve_rows)
    sample_summary <- do.call(rbind, summary_rows)

    curve_data <- merge(curve_data, metadata, by = "sample")
    sample_summary <- merge(sample_summary, metadata, by = "sample")

    if (!("library" %in% names(sample_summary))) {
        sample_summary$library <- config$project_name
        curve_data$library <- config$project_name
    }

    output_dir <- file.path(config$output_directory, "QC", "Saturation")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        curve_data,
        file.path(table_dir, "saturation_curves.csv"),
        row.names = FALSE
    )

    save_csv(
        sample_summary,
        file.path(table_dir, "sample_summary.csv"),
        row.names = FALSE
    )

    group_summary <- stats::aggregate(
        cbind(expected_genes_at_1M, full_depth_endogenous_genes) ~ group,
        data = sample_summary,
        FUN = mean,
        na.action = stats::na.pass
    )

    names(group_summary)[names(group_summary) == "expected_genes_at_1M"] <- "mean_genes_at_1M"
    names(group_summary)[names(group_summary) == "full_depth_endogenous_genes"] <- "mean_full_depth_genes"

    group_summary$samples <- as.vector(table(sample_summary$group)[group_summary$group])

    save_csv(
        group_summary,
        file.path(table_dir, "group_summary.csv"),
        row.names = FALSE
    )

    curve_data$depth_millions <- curve_data$depth / 1e6

    libraries <- sort(unique(curve_data$library))

    saturation_plot <- ggplot2::ggplot(
        curve_data,
        ggplot2::aes(x = depth_millions, y = detected_genes, group = sample)
    ) +
        ggplot2::geom_line(color = "#447899", alpha = 0.3, linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = 1, color = "#222222", linetype = "dashed", linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = 10000, color = "#b35c20", linetype = "dashed", linewidth = 0.4) +
        ggplot2::facet_wrap(~library) +
        ggplot2::labs(
            title = sprintf("Gene Detection Saturation by Library | %d samples", length(samples)),
            x = "Total reads per sample (millions; provider QC depth)",
            y = "Expected detected endogenous genes"
        ) +
        ggplot2::theme_bw()

    n_libraries <- length(libraries)

    plot_width <- max(config$figure_width, min(n_libraries, 4) * 3)
    plot_height <- max(config$figure_height, ceiling(n_libraries / 4) * 3)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "Saturation_curve.png"),
            plot = saturation_plot,
            width = plot_width,
            height = plot_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "Saturation_curve.pdf"),
            plot = saturation_plot,
            width = plot_width,
            height = plot_height
        )

    }

    message("Gene detection saturation curve complete.")

    return(
        list(
            curves = curve_data,
            sample_summary = sample_summary,
            group_summary = group_summary,
            plot = saturation_plot
        )
    )

}
