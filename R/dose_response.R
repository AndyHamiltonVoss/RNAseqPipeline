################################################################################
# RNAseqPipeline
#
# File: dose_response.R
#
# Description:
#   Cross-method dose-response comparison: counts up/downregulated genes
#   per dose for DESeq2, edgeR, and limma-voom, plus a consensus set of
#   genes called concordantly by all three, and plots the results.
#
################################################################################

#===============================================================================
# Parse Compound and Dose From a Group Name
#===============================================================================

parse_dose_group <- function(group_name) {

    match <- regmatches(
        group_name,
        regexec("^([A-Za-z0-9]+)_([0-9.]+)(pM|nM|uM|mM)$", group_name)
    )[[1]]

    if (length(match) != 4) {
        return(NULL)
    }

    unit_to_nM <- c(pM = 0.001, nM = 1, uM = 1000, mM = 1e6)

    list(
        compound = match[2],
        dose_nM = as.numeric(match[3]) * unit_to_nM[[match[4]]]
    )

}

#===============================================================================
# Build Dose Series for a Compound
#===============================================================================

get_dose_series <- function(groups, compound) {

    parsed <- lapply(groups, parse_dose_group)

    keep <- vapply(
        parsed,
        function(x) !is.null(x) && identical(x$compound, compound),
        logical(1)
    )

    series <- data.frame(
        group = groups[keep],
        dose_nM = vapply(parsed[keep], function(x) x$dose_nM, numeric(1)),
        stringsAsFactors = FALSE
    )

    series <- series[order(series$dose_nM), ]

    rownames(series) <- NULL

    return(series)

}

#===============================================================================
# Classify Regulated Genes From a Results Table
#===============================================================================

classify_genes <- function(
    results,
    gene_col,
    lfc_col,
    padj_col,
    alpha,
    log2fc_cutoff
) {

    sig <- !is.na(results[[padj_col]]) & results[[padj_col]] < alpha

    up <- results[[gene_col]][sig & results[[lfc_col]] > log2fc_cutoff]

    down <- results[[gene_col]][sig & results[[lfc_col]] < -log2fc_cutoff]

    list(up = up, down = down)

}

#===============================================================================
# Load Per-Method Results for One Contrast
#===============================================================================

load_method_results <- function(group, output_directory, reference_group) {

    deseq2_file <- file.path(
        output_directory, "DESeq2", "Tables",
        paste0("DESeq2_group_", group, "_vs_", reference_group, "_shrunk.csv")
    )

    edger_file <- file.path(
        output_directory, "edgeR", "Tables",
        paste0("edgeR_", group, "_vs_", reference_group, "_raw.csv")
    )

    limma_file <- file.path(
        output_directory, "limma", "Tables",
        paste0("limma_", group, "_vs_", reference_group, "_raw.csv")
    )

    list(
        DESeq2 = if (file.exists(deseq2_file)) read.csv(deseq2_file) else NULL,
        edgeR = if (file.exists(edger_file)) read.csv(edger_file) else NULL,
        limma = if (file.exists(limma_file)) read.csv(limma_file) else NULL
    )

}

#===============================================================================
# Summarize Gene Counts Across Methods for One Contrast
#===============================================================================

summarize_contrast_counts <- function(
    group,
    dose_nM,
    output_directory,
    reference_group,
    alpha,
    log2fc_cutoff
) {

    results <- load_method_results(group, output_directory, reference_group)

    empty <- list(up = character(0), down = character(0))

    classified <- list(
        DESeq2 = if (!is.null(results$DESeq2))
            classify_genes(results$DESeq2, "gene", "log2FoldChange", "padj", alpha, log2fc_cutoff)
        else empty,

        edgeR = if (!is.null(results$edgeR))
            classify_genes(results$edgeR, "gene", "logFC", "FDR", alpha, log2fc_cutoff)
        else empty,

        limma = if (!is.null(results$limma))
            classify_genes(results$limma, "gene", "logFC", "adj.P.Val", alpha, log2fc_cutoff)
        else empty
    )

    consensus_up <- Reduce(
        intersect,
        list(classified$DESeq2$up, classified$edgeR$up, classified$limma$up)
    )

    consensus_down <- Reduce(
        intersect,
        list(classified$DESeq2$down, classified$edgeR$down, classified$limma$down)
    )

    counts <- data.frame(
        group = group,
        dose_nM = dose_nM,
        method = c("DESeq2", "edgeR", "limma", "Consensus"),
        n_up = c(
            length(classified$DESeq2$up),
            length(classified$edgeR$up),
            length(classified$limma$up),
            length(consensus_up)
        ),
        n_down = c(
            length(classified$DESeq2$down),
            length(classified$edgeR$down),
            length(classified$limma$down),
            length(consensus_down)
        ),
        stringsAsFactors = FALSE
    )

    return(counts)

}

#===============================================================================
# Build Dose-Response Summary Table for a Compound
#===============================================================================

build_dose_response_table <- function(
    compound,
    groups,
    output_directory,
    reference_group,
    alpha,
    log2fc_cutoff
) {

    series <- get_dose_series(groups, compound)

    if (nrow(series) == 0) {
        warning("No dose series found for compound: ", compound)
        return(NULL)
    }

    counts_list <- lapply(
        seq_len(nrow(series)),
        function(i)
            summarize_contrast_counts(
                group = series$group[i],
                dose_nM = series$dose_nM[i],
                output_directory = output_directory,
                reference_group = reference_group,
                alpha = alpha,
                log2fc_cutoff = log2fc_cutoff
            )
    )

    counts <- do.call(rbind, counts_list)

    counts$compound <- compound

    return(counts)

}

#===============================================================================
# Plot Dose-Response Curves
#===============================================================================

plot_dose_response <- function(counts, direction = c("up", "down")) {

    direction <- match.arg(direction)

    y_col <- if (direction == "up") "n_up" else "n_down"

    plot_data <- counts
    plot_data$n_genes <- plot_data[[y_col]]

    method_colors <- c(
        DESeq2 = "#1b9e77",
        edgeR = "#d95f02",
        limma = "#7570b3",
        Consensus = "#000000"
    )

    line_data <- plot_data[plot_data$method != "Consensus", ]
    consensus_data <- plot_data[plot_data$method == "Consensus", ]

    plot <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = dose_nM,
            y = n_genes,
            color = method,
            group = method
        )
    ) +
        ggplot2::geom_line(
            data = line_data,
            linewidth = 0.8
        ) +
        ggplot2::geom_point(
            data = line_data,
            size = 2
        ) +
        ggplot2::geom_line(
            data = consensus_data,
            linewidth = 1,
            linetype = "dashed"
        ) +
        ggplot2::geom_point(
            data = consensus_data,
            size = 2.5,
            shape = 17
        ) +
        ggplot2::facet_wrap(~compound, scales = "free_x") +
        ggplot2::scale_x_log10() +
        ggplot2::scale_color_manual(values = method_colors) +
        ggplot2::labs(
            title = paste0(
                if (direction == "up") "Upregulated" else "Downregulated",
                " Genes vs. Dose"
            ),
            x = "Dose (nM, log10 scale)",
            y = paste0("Number of ", direction, "regulated genes"),
            color = "Method"
        ) +
        ggplot2::theme_bw()

    return(plot)

}

#===============================================================================
# Run Full Dose-Response Comparison
#===============================================================================

run_dose_response <- function(config, metadata) {

    message("Running dose-response comparison...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    if (!isTRUE(config$run_deseq2) || !isTRUE(config$run_edger) || !isTRUE(config$run_limma)) {
        warning(
            "Dose-response comparison requires DESeq2, edgeR, and limma to all be ",
            "enabled; results from disabled methods are treated as empty gene sets, ",
            "so the Consensus line will be empty."
        )
    }

    groups <- setdiff(
        levels(factor(metadata$group)),
        config$reference_group
    )

    output_dir <- file.path(config$output_directory, "Figures")

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    all_counts <- do.call(
        rbind,
        lapply(
            config$dose_response_compounds,
            function(compound)
                build_dose_response_table(
                    compound = compound,
                    groups = groups,
                    output_directory = config$output_directory,
                    reference_group = config$reference_group,
                    alpha = config$alpha,
                    log2fc_cutoff = config$log2fc_cutoff
                )
        )
    )

    if (is.null(all_counts) || nrow(all_counts) == 0) {
        warning("No dose-response series found; skipping dose-response plots.")
        return(invisible(NULL))
    }

    save_csv(
        all_counts,
        file.path(output_dir, "dose_response_gene_counts.csv"),
        row.names = FALSE
    )

    up_plot <- plot_dose_response(all_counts, direction = "up")
    down_plot <- plot_dose_response(all_counts, direction = "down")

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(output_dir, "dose_response_upregulated.png"),
            plot = up_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

        ggplot2::ggsave(
            file.path(output_dir, "dose_response_downregulated.png"),
            plot = down_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(output_dir, "dose_response_upregulated.pdf"),
            plot = up_plot,
            width = config$figure_width,
            height = config$figure_height
        )

        ggplot2::ggsave(
            file.path(output_dir, "dose_response_downregulated.pdf"),
            plot = down_plot,
            width = config$figure_width,
            height = config$figure_height
        )

    }

    message("Dose-response comparison complete.")

    return(
        list(
            counts = all_counts,
            up_plot = up_plot,
            down_plot = down_plot
        )
    )

}
