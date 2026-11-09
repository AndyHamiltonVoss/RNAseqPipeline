################################################################################
# RNAseqPipeline
#
# File: dose_concordance.R
#
# Description:
#   Dose-to-dose log2FC concordance analysis: for a compound's dose series,
#   checks whether genes that are called as consensus-regulated at only some
#   doses ("flicker" in/out of the Venn overlap) sit near the significance
#   threshold with correlated effect sizes (noise at the cutoff), or show
#   genuinely different magnitudes across doses (real dose-dependent
#   recruitment). Reuses get_dose_series() / load_method_results() from
#   R/dose_response.R and get_consensus_genes() / format_direction_label()
#   from R/compound_venn.R.
#
################################################################################

#===============================================================================
# Build a Gene x Dose log2FC Matrix (DESeq2 Shrunk log2FC vs Reference)
#===============================================================================

build_dose_log2fc_matrix <- function(series, output_directory, reference_group) {

    log2fc_list <- lapply(
        series$group,
        function(group) {

            results <- load_method_results(group, output_directory, reference_group)

            if (is.null(results$DESeq2)) {
                stop("DESeq2 results required for log2FC concordance (group: ", group, ")")
            }

            stats::setNames(results$DESeq2$log2FoldChange, results$DESeq2$gene)

        }
    )

    all_genes <- Reduce(union, lapply(log2fc_list, names))

    log2fc_matrix <- vapply(
        log2fc_list,
        function(x) unname(x[all_genes]),
        numeric(length(all_genes))
    )

    rownames(log2fc_matrix) <- all_genes

    return(log2fc_matrix)

}

#===============================================================================
# Plot Pairwise Dose-to-Dose log2FC Concordance for a Compound
#===============================================================================

plot_dose_log2fc_concordance <- function(
    compound,
    config,
    metadata,
    direction = "up",
    log2fc_cutoff = config$log2fc_cutoff
) {

    groups <- setdiff(
        levels(factor(metadata$group)),
        config$reference_group
    )

    series <- get_dose_series(groups, compound)

    if (nrow(series) < 2) {
        stop("Need at least 2 doses for compound '", compound, "' to plot concordance.")
    }

    dose_labels <- sub(paste0("^", compound, "_"), "", series$group)

    called_genes <- stats::setNames(
        lapply(series$group, function(g) get_consensus_genes(g, direction, config, log2fc_cutoff)),
        dose_labels
    )

    ever_called <- Reduce(union, called_genes)

    log2fc_matrix <- build_dose_log2fc_matrix(series, config$output_directory, config$reference_group)
    colnames(log2fc_matrix) <- dose_labels

    log2fc_matrix <- log2fc_matrix[rownames(log2fc_matrix) %in% ever_called, , drop = FALSE]

    dose_pairs <- utils::combn(dose_labels, 2, simplify = FALSE)

    plot_data <- do.call(
        rbind,
        lapply(
            dose_pairs,
            function(pair) {

                x_dose <- pair[1]
                y_dose <- pair[2]

                called_x <- rownames(log2fc_matrix) %in% called_genes[[x_dose]]
                called_y <- rownames(log2fc_matrix) %in% called_genes[[y_dose]]

                call_status <- ifelse(
                    called_x & called_y, "Called in both",
                    ifelse(called_x | called_y, "Called in one", "Called in neither")
                )

                data.frame(
                    gene = rownames(log2fc_matrix),
                    x = log2fc_matrix[, x_dose],
                    y = log2fc_matrix[, y_dose],
                    call_status = call_status,
                    pair = paste0(x_dose, " vs ", y_dose),
                    stringsAsFactors = FALSE
                )

            }
        )
    )

    plot_data$pair <- factor(plot_data$pair, levels = unique(plot_data$pair))

    status_colors <- c(
        "Called in both" = "#1b9e77",
        "Called in one" = "#d95f02",
        "Called in neither" = "grey70"
    )

    cutoff_lines <- if (identical(direction, "combined")) {
        c(-log2fc_cutoff, log2fc_cutoff)
    } else if (identical(direction, "down")) {
        -log2fc_cutoff
    } else {
        log2fc_cutoff
    }

    plot <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = x, y = y, color = call_status)
    ) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::geom_hline(yintercept = cutoff_lines, linetype = "dotted", color = "grey50") +
        ggplot2::geom_vline(xintercept = cutoff_lines, linetype = "dotted", color = "grey50") +
        ggplot2::geom_point(alpha = 0.6, size = 1.3) +
        ggplot2::facet_wrap(~pair) +
        ggplot2::scale_color_manual(values = status_colors) +
        ggplot2::labs(
            title = paste0(
                compound, ": Dose-to-Dose log2FC Concordance (",
                format_direction_label(direction), ")"
            ),
            subtitle = "DESeq2 shrunk log2FC vs DMSO; dotted lines mark the log2FC cutoff used for calls",
            x = "log2FC (dose on x-axis)",
            y = "log2FC (dose on y-axis)",
            color = "Consensus call"
        ) +
        ggplot2::theme_bw()

    return(list(plot = plot, data = plot_data))

}

#===============================================================================
# Run Dose-to-Dose log2FC Concordance Analysis for a Compound
#===============================================================================

run_dose_log2fc_concordance <- function(compound, config, metadata, direction = "up") {

    message("Running log2FC concordance analysis: ", compound, " (", direction, ")")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    result <- plot_dose_log2fc_concordance(compound, config, metadata, direction)

    figure_dir <- file.path(config$output_directory, "Figures", "DoseConcordance")

    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    filename_base <- paste0(compound, "_log2FC_concordance_", direction)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_base, ".png")),
            plot = result$plot,
            width = config$figure_width * 1.3,
            height = config$figure_height * 1.3,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_base, ".pdf")),
            plot = result$plot,
            width = config$figure_width * 1.3,
            height = config$figure_height * 1.3
        )

    }

    save_csv(
        result$data,
        file.path(figure_dir, paste0(filename_base, "_data.csv")),
        row.names = FALSE
    )

    message("log2FC concordance analysis complete.")

    return(result)

}
