################################################################################
# RNAseqPipeline
#
# File: pathway_downregulation_skew.R
#
# Description:
#   Diagnostic: does GSVA KEGG pathway activity (diff_score vs. the
#   reference group, from the per-contrast limma results already written by
#   run_gsva() in R/gsva.R) skew toward down-regulation at each compound's
#   top tested dose? Two views:
#     - Per-compound % of all scored pathways with negative diff_score (not
#       just the ones crossing the significance cutoff).
#     - Pooled distribution of diff_score across every compound x pathway
#       combination.
#   This was added after finding that pathway-level compound-vs-compound
#   similarity (R/gsva_compound_correlation.R) is dominated by a shared,
#   dose-dependent, largely compound-nonspecific negative shift across most
#   pathways -- this figure is the diagnostic evidence for that shift.
#
################################################################################

#===============================================================================
# Run Pathway Down-Regulation Skew Diagnostic
#===============================================================================

run_pathway_downregulation_skew <- function(config) {

    message("Running pathway down-regulation skew diagnostic...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    deseq2_table_dir <- file.path(config$output_directory, "DESeq2", "Tables")
    gsva_table_dir <- file.path(config$output_directory, "GSVA", "Tables")
    diagnostic_dir <- file.path(config$output_directory, "pathway_diagnostic")

    #===========================================================================
    # Discover Each Compound's Top Dose
    #===========================================================================

    shrunk_files <- list.files(
        deseq2_table_dir,
        pattern = "^DESeq2_group_.*_shrunk\\.csv$"
    )

    group_names <- sub(
        paste0("_vs_", config$reference_group, "_shrunk\\.csv$"), "",
        sub("^DESeq2_group_", "", shrunk_files)
    )

    parsed <- lapply(group_names, parse_dose_group)

    keep <- !vapply(parsed, is.null, logical(1))

    dose_table <- data.frame(
        group = group_names[keep],
        compound = vapply(parsed[keep], `[[`, character(1), "compound"),
        dose_nM = vapply(parsed[keep], `[[`, numeric(1), "dose_nM")
    )

    top_dose_table <- do.call(
        rbind,
        lapply(split(dose_table, dose_table$compound), function(d) {
            d[which.max(d$dose_nM), ]
        })
    )

    gsva_file_for <- function(group) {
        file.path(gsva_table_dir, paste0("GSVA_KEGG_", group, "_vs_", config$reference_group, "_limma.csv"))
    }

    top_dose_table$has_gsva <- vapply(top_dose_table$group, function(g) file.exists(gsva_file_for(g)), logical(1))
    eligible <- top_dose_table[top_dose_table$has_gsva, ]
    compounds <- sort(eligible$compound)

    message(sprintf("Checking pathway diff_score skew across %d compounds at their highest tested dose.", length(compounds)))

    #===========================================================================
    # Read diff_score for All Pathways (Not Just Significant Ones)
    #===========================================================================

    long_data <- do.call(rbind, lapply(compounds, function(cmpd) {
        group <- eligible$group[eligible$compound == cmpd]
        res <- read.csv(gsva_file_for(group))
        data.frame(compound = cmpd, pathway = res$pathway, diff_score = res$diff_score)
    }))

    #===========================================================================
    # Per-Compound % of Pathways With Negative diff_score
    #===========================================================================

    pct_negative_df <- do.call(rbind, lapply(compounds, function(cmpd) {
        sub <- long_data[long_data$compound == cmpd, ]
        data.frame(compound = cmpd, pct_negative = mean(sub$diff_score < 0) * 100)
    }))

    overall_mean_pct_negative <- mean(pct_negative_df$pct_negative)

    pct_negative_df <- pct_negative_df[order(-pct_negative_df$pct_negative), ]

    #===========================================================================
    # Pooled diff_score Distribution
    #===========================================================================

    pooled_median <- stats::median(long_data$diff_score)
    pooled_pct_negative <- mean(long_data$diff_score < 0) * 100

    #===========================================================================
    # Save Tables
    #===========================================================================

    dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        pct_negative_df,
        file.path(diagnostic_dir, "pathway_pct_negative_by_compound.csv"),
        row.names = FALSE
    )

    pooled_summary <- data.frame(
        n_compounds = length(compounds),
        n_pathways = length(unique(long_data$pathway)),
        pooled_pct_negative = pooled_pct_negative,
        pooled_median_diff_score = pooled_median,
        overall_mean_pct_negative_across_compounds = overall_mean_pct_negative
    )

    save_csv(
        pooled_summary,
        file.path(diagnostic_dir, "pathway_downregulation_skew_summary.csv"),
        row.names = FALSE
    )

    #===========================================================================
    # Plot: Per-Compound % Negative Bar Chart
    #===========================================================================

    pct_negative_df$compound <- factor(pct_negative_df$compound, levels = pct_negative_df$compound)

    panel1 <- ggplot2::ggplot(pct_negative_df, ggplot2::aes(x = compound, y = pct_negative)) +
        ggplot2::geom_col(fill = "steelblue") +
        ggplot2::geom_hline(yintercept = 50, linetype = "dashed", color = "grey30") +
        ggplot2::geom_hline(yintercept = overall_mean_pct_negative, linetype = "dotted", color = "firebrick") +
        ggplot2::annotate(
            "text", x = Inf, y = overall_mean_pct_negative, hjust = 1.05, vjust = -0.5,
            label = sprintf("mean = %.1f%%", overall_mean_pct_negative), color = "firebrick", size = 3
        ) +
        ggplot2::labs(
            title = "Most KEGG pathways trend down at top dose, across nearly every compound",
            subtitle = sprintf(
                "%% of all %d scored pathways with negative diff_score (not just significant ones); dashed line = 50%% (no skew)",
                length(unique(long_data$pathway))
            ),
            x = "Compound", y = "% of pathways with negative diff_score"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    #===========================================================================
    # Plot: Pooled diff_score Histogram
    #===========================================================================

    panel2 <- ggplot2::ggplot(long_data, ggplot2::aes(x = diff_score)) +
        ggplot2::geom_histogram(bins = 80, fill = "steelblue", color = NA) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
        ggplot2::geom_vline(xintercept = pooled_median, linetype = "dotted", color = "firebrick") +
        ggplot2::labs(
            title = "Pooled distribution of diff_score across all compounds and pathways",
            subtitle = sprintf(
                "%d compounds x %d pathways; %.1f%% negative; dashed = 0, dotted = pooled median (%.3f)",
                length(compounds), length(unique(long_data$pathway)), pooled_pct_negative, pooled_median
            ),
            x = "diff_score (pathway activity vs. reference group)", y = "Count"
        ) +
        ggplot2::theme_bw()

    #===========================================================================
    # Combine and Save
    #===========================================================================

    if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- patchwork::wrap_plots(panel1, panel2, ncol = 1, heights = c(1.1, 1))
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
        combined <- cowplot::plot_grid(panel1, panel2, ncol = 1, rel_heights = c(1.1, 1))
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
        combined <- gridExtra::arrangeGrob(panel1, panel2, ncol = 1, heights = c(1.1, 1))
    } else {
        stop("One of 'patchwork', 'cowplot', or 'gridExtra' is required to combine the two panels.")
    }

    plot_width <- max(config$figure_width, 11)
    plot_height <- max(config$figure_height, 11)

    if (isTRUE(config$save_png)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Pathway_downregulation_skew.png"),
            plot = combined, width = plot_width, height = plot_height, dpi = config$figure_dpi
        )
    }

    if (isTRUE(config$save_pdf)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Pathway_downregulation_skew.pdf"),
            plot = combined, width = plot_width, height = plot_height
        )
    }

    message("Pathway down-regulation skew diagnostic complete.")

    return(
        list(
            pct_negative_by_compound = pct_negative_df,
            pooled_summary = pooled_summary
        )
    )

}
