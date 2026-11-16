################################################################################
# RNAseqPipeline
#
# File: pathway_specificity_comparison.R
#
# Description:
#   Diagnostic: are up-regulated genes concentrated in a specific, shared
#   set of pathways across compounds, while down-regulated genes are
#   comparatively scattered/non-specific? Mines the GO and KEGG enrichment
#   results already written by run_enrichment() (R/enrichment.R) for each
#   compound's top tested dose, separately for the up- and down-regulated
#   gene sets, and asks: how many compounds share each significant term?
#   A term private to one compound reflects scattered/idiosyncratic signal;
#   a term shared across many compounds reflects convergent, specific
#   biology. This complements R/compound_similarity.R's gene-identity-level
#   finding (up-Jaccard >> down-Jaccard between compounds) with the same
#   question asked at the pathway-annotation level.
#
################################################################################

#===============================================================================
# Run Pathway Specificity Comparison (Up vs. Down)
#===============================================================================

run_pathway_specificity_comparison <- function(config) {

    message("Running pathway specificity comparison (up vs. down)...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    deseq2_table_dir <- file.path(config$output_directory, "DESeq2", "Tables")
    go_table_dir <- file.path(config$output_directory, "GO", "Tables")
    kegg_table_dir <- file.path(config$output_directory, "KEGG", "Tables")
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

    exclude_pattern <- config$exclude_groups_pattern

    if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {
        dose_table <- dose_table[!grepl(exclude_pattern, dose_table$group), ]
    }

    top_dose_table <- do.call(
        rbind,
        lapply(split(dose_table, dose_table$compound), function(d) {
            d[which.max(d$dose_nM), ]
        })
    )

    compounds <- sort(top_dose_table$compound)

    message(sprintf("Comparing pathway term specificity across %d compounds at their highest tested dose.", length(compounds)))

    #===========================================================================
    # Read Significant Term Sets (Up / Down) Per Compound, Per Source
    #===========================================================================

    read_terms <- function(path) {
        if (!file.exists(path)) {
            return(character(0))
        }
        result <- tryCatch(read.csv(path), error = function(e) NULL)
        if (is.null(result) || nrow(result) == 0 || !("Description" %in% names(result))) {
            return(character(0))
        }
        unique(result$Description)
    }

    build_term_sets <- function(table_dir, source_label) {
        sets <- lapply(compounds, function(cmpd) {
            group <- top_dose_table$group[top_dose_table$compound == cmpd]
            list(
                up = read_terms(file.path(table_dir, paste0(source_label, "_", group, "_up.csv"))),
                down = read_terms(file.path(table_dir, paste0(source_label, "_", group, "_down.csv")))
            )
        })
        names(sets) <- compounds
        sets
    }

    kegg_sets <- build_term_sets(kegg_table_dir, "KEGG")
    go_sets <- build_term_sets(go_table_dir, "GO")

    #===========================================================================
    # Build Term-Frequency Tables (How Many Compounds Share Each Term)
    #===========================================================================

    build_frequency_df <- function(term_sets, source_label) {

        up_freq <- table(unlist(lapply(term_sets, `[[`, "up")))
        down_freq <- table(unlist(lapply(term_sets, `[[`, "down")))

        rbind(
            if (length(up_freq) > 0) {
                data.frame(source = source_label, direction = "Up", term = names(up_freq), n_compounds = as.integer(up_freq), stringsAsFactors = FALSE)
            },
            if (length(down_freq) > 0) {
                data.frame(source = source_label, direction = "Down", term = names(down_freq), n_compounds = as.integer(down_freq), stringsAsFactors = FALSE)
            }
        )

    }

    kegg_frequency_df <- build_frequency_df(kegg_sets, "KEGG")
    go_frequency_df <- build_frequency_df(go_sets, "GO")

    #===========================================================================
    # Specificity Metrics Per Source x Direction
    #===========================================================================

    summarize_direction <- function(freq_df, source_label, direction_label) {

        sub <- freq_df[freq_df$direction == direction_label, ]

        if (nrow(sub) == 0) {
            return(data.frame(
                source = source_label, direction = direction_label,
                n_unique_terms = 0L, total_calls = 0L, pct_private = NA_real_,
                top10_concentration_pct = NA_real_, mean_frequency = NA_real_, median_frequency = NA_real_
            ))
        }

        freq_sorted <- sort(sub$n_compounds, decreasing = TRUE)

        data.frame(
            source = source_label,
            direction = direction_label,
            n_unique_terms = nrow(sub),
            total_calls = sum(sub$n_compounds),
            pct_private = mean(sub$n_compounds == 1) * 100,
            top10_concentration_pct = sum(head(freq_sorted, 10)) / sum(freq_sorted) * 100,
            mean_frequency = mean(sub$n_compounds),
            median_frequency = stats::median(sub$n_compounds)
        )

    }

    specificity_summary <- rbind(
        summarize_direction(kegg_frequency_df, "KEGG", "Up"),
        summarize_direction(kegg_frequency_df, "KEGG", "Down"),
        summarize_direction(go_frequency_df, "GO", "Up"),
        summarize_direction(go_frequency_df, "GO", "Down")
    )

    #===========================================================================
    # Save Tables
    #===========================================================================

    dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(specificity_summary, file.path(diagnostic_dir, "pathway_specificity_summary.csv"), row.names = FALSE)
    save_csv(kegg_frequency_df, file.path(diagnostic_dir, "pathway_specificity_kegg_term_frequency.csv"), row.names = FALSE)
    save_csv(go_frequency_df, file.path(diagnostic_dir, "pathway_specificity_go_term_frequency.csv"), row.names = FALSE)

    #===========================================================================
    # Plot: Term-Sharing Frequency Distribution, KEGG and GO
    #===========================================================================

    make_sharing_plot <- function(freq_df, source_label, n_compounds) {

        freq_df$direction <- factor(freq_df$direction, levels = c("Down", "Up"))

        ggplot2::ggplot(freq_df, ggplot2::aes(x = n_compounds, fill = direction)) +
            ggplot2::geom_bar(position = ggplot2::position_dodge2(preserve = "single", padding = 0.1)) +
            ggplot2::scale_fill_manual(values = c(Down = "steelblue", Up = "darkorange")) +
            ggplot2::scale_x_continuous(breaks = seq_len(n_compounds)) +
            ggplot2::labs(
                title = sprintf("%s term sharing: how many compounds call each term significant?", source_label),
                x = "Number of compounds the term is significant in",
                y = "Number of unique terms",
                fill = "Direction"
            ) +
            ggplot2::theme_bw()

    }

    kegg_plot <- make_sharing_plot(kegg_frequency_df, "KEGG", length(compounds))
    go_plot <- make_sharing_plot(go_frequency_df, "GO", length(compounds))

    if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- patchwork::wrap_plots(kegg_plot, go_plot, ncol = 1)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
        combined <- cowplot::plot_grid(kegg_plot, go_plot, ncol = 1)
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
        combined <- gridExtra::arrangeGrob(kegg_plot, go_plot, ncol = 1)
    } else {
        stop("One of 'patchwork', 'cowplot', or 'gridExtra' is required to combine the two panels.")
    }

    plot_width <- max(config$figure_width, 11)
    plot_height <- max(config$figure_height, 11)

    if (isTRUE(config$save_png)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Pathway_specificity_term_sharing.png"),
            plot = combined, width = plot_width, height = plot_height, dpi = config$figure_dpi
        )
    }

    if (isTRUE(config$save_pdf)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Pathway_specificity_term_sharing.pdf"),
            plot = combined, width = plot_width, height = plot_height
        )
    }

    message("Pathway specificity comparison complete.")

    return(
        list(
            specificity_summary = specificity_summary,
            kegg_frequency = kegg_frequency_df,
            go_frequency = go_frequency_df
        )
    )

}
