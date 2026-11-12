################################################################################
# RNAseqPipeline
#
# File: de_gene_heatmap.R
#
# Description:
#   Heatmap of the number of differentially expressed genes (DESeq2,
#   apeglm-shrunk log2FoldChange) by drug x dose. Reuses parse_dose_group()
#   (R/dose_response.R) to assign each contrast to a (compound, dose) cell,
#   which naturally excludes the reference group and co-treatment contrasts,
#   classify_genes() (R/dose_response.R) for the significance/fold-change
#   call, and format_dose_label() (R/ercc_dose_heatmap.R) for the per-cell
#   dose annotation.
#
################################################################################

#===============================================================================
# Run DE Gene Count Heatmap by Drug x Dose
#===============================================================================

run_de_gene_heatmap <- function(config) {

    message("Running DE gene count heatmap by drug x dose...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    table_dir <- file.path(config$output_directory, "DESeq2", "Tables")
    figure_dir <- file.path(config$output_directory, "DESeq2", "Figures")

    shrunk_files <- list.files(
        table_dir,
        pattern = "^DESeq2_group_.*_shrunk\\.csv$"
    )

    group_names <- sub(
        paste0("_vs_", config$reference_group, "_shrunk\\.csv$"), "",
        sub("^DESeq2_group_", "", shrunk_files)
    )

    parsed <- lapply(group_names, parse_dose_group)

    keep <- !vapply(parsed, is.null, logical(1))

    message(sprintf(
        "%d of %d contrasts excluded (no single-compound dose parsed, e.g. co-treatment arms).",
        sum(!keep), length(keep)
    ))

    if (sum(keep) == 0) {
        stop("No contrasts with a single-compound dose were found; nothing to plot.")
    }

    cell_rows <- lapply(seq_along(group_names)[keep], function(i) {

        group <- group_names[i]

        results <- read.csv(
            file.path(
                table_dir,
                paste0("DESeq2_group_", group, "_vs_", config$reference_group, "_shrunk.csv")
            )
        )

        classified <- classify_genes(
            results, "gene", "log2FoldChange", "padj",
            config$alpha, config$log2fc_cutoff
        )

        data.frame(
            compound = parsed[[i]]$compound,
            dose_nM = parsed[[i]]$dose_nM,
            n_up = length(classified$up),
            n_down = length(classified$down),
            n_deg = length(classified$up) + length(classified$down)
        )

    })

    cell_data <- do.call(rbind, cell_rows)

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

    cell_data$log2_up_down_ratio <- ifelse(
        cell_data$n_up + cell_data$n_down == 0,
        NA_real_,
        log2((cell_data$n_up + 1) / (cell_data$n_down + 1))
    )

    cell_data$ratio_label <- ifelse(
        is.na(cell_data$log2_up_down_ratio), "",
        sprintf("%.2f", cell_data$log2_up_down_ratio)
    )

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        cell_data,
        file.path(table_dir, "de_gene_counts_by_dose.csv"),
        row.names = FALSE
    )

    heatmap_plot <- ggplot2::ggplot(
        cell_data,
        ggplot2::aes(x = dose_rank, y = compound, fill = n_deg)
    ) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = n_deg), size = 2.8) +
        ggplot2::scale_fill_viridis_c(name = "DE genes (n)") +
        ggplot2::labs(
            title = "Differentially Expressed Gene Count by Drug and Dose",
            subtitle = sprintf(
                "DESeq2, apeglm-shrunk LFC; alpha = %s, |log2FC| > %s",
                config$alpha, config$log2fc_cutoff
            ),
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
            file.path(figure_dir, "DE_gene_count_heatmap.png"),
            plot = heatmap_plot,
            width = plot_width,
            height = plot_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "DE_gene_count_heatmap.pdf"),
            plot = heatmap_plot,
            width = plot_width,
            height = plot_height
        )

    }

    updown_data <- rbind(
        data.frame(
            compound = cell_data$compound,
            dose_rank = cell_data$dose_rank,
            direction = "Upregulated",
            n_genes = cell_data$n_up
        ),
        data.frame(
            compound = cell_data$compound,
            dose_rank = cell_data$dose_rank,
            direction = "Downregulated",
            n_genes = cell_data$n_down
        )
    )

    updown_data$direction <- factor(
        updown_data$direction,
        levels = c("Upregulated", "Downregulated")
    )

    updown_plot <- ggplot2::ggplot(
        updown_data,
        ggplot2::aes(x = dose_rank, y = compound, fill = n_genes)
    ) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = n_genes), size = 2.8) +
        ggplot2::scale_fill_viridis_c(name = "DE genes (n)") +
        ggplot2::facet_wrap(~direction) +
        ggplot2::labs(
            title = "Up- vs. Downregulated Gene Count by Drug and Dose",
            subtitle = sprintf(
                "DESeq2, apeglm-shrunk LFC; alpha = %s, |log2FC| > %s",
                config$alpha, config$log2fc_cutoff
            ),
            x = "Dose (ordinal rank, lowest to highest per compound)",
            y = "Compound"
        ) +
        ggplot2::theme_bw()

    updown_width <- max(config$figure_width * 1.8, n_doses * 1.2 * 2)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "DE_gene_count_heatmap_up_down.png"),
            plot = updown_plot,
            width = updown_width,
            height = plot_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "DE_gene_count_heatmap_up_down.pdf"),
            plot = updown_plot,
            width = updown_width,
            height = plot_height
        )

    }

    ratio_plot <- ggplot2::ggplot(
        cell_data,
        ggplot2::aes(x = dose_rank, y = compound, fill = log2_up_down_ratio)
    ) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = ratio_label), size = 2.8) +
        ggplot2::scale_fill_gradient2(
            name = "log2(up/down)\n(+1 pseudocount)",
            low = "#2166ac", mid = "white", high = "#b2182b",
            midpoint = 0, na.value = "grey80"
        ) +
        ggplot2::labs(
            title = "Up/Down-regulated Gene Ratio by Drug and Dose",
            subtitle = sprintf(
                "log2((up+1)/(down+1)); DESeq2, apeglm-shrunk LFC; alpha = %s, |log2FC| > %s. Grey = no DE genes either direction.",
                config$alpha, config$log2fc_cutoff
            ),
            x = "Dose (ordinal rank, lowest to highest per compound)",
            y = "Compound"
        ) +
        ggplot2::theme_bw()

    ratio_plot_width <- max(plot_width, 11)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "DE_gene_up_down_ratio_heatmap.png"),
            plot = ratio_plot,
            width = ratio_plot_width,
            height = plot_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "DE_gene_up_down_ratio_heatmap.pdf"),
            plot = ratio_plot,
            width = ratio_plot_width,
            height = plot_height
        )

    }

    message("DE gene count heatmap complete.")

    return(
        list(
            cell_data = cell_data,
            plot = heatmap_plot,
            updown_plot = updown_plot,
            ratio_plot = ratio_plot
        )
    )

}
