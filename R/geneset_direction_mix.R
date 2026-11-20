################################################################################
# RNAseqPipeline
#
# File: geneset_direction_mix.R
#
# Description:
#   For each compound and each KEGG pathway, checks whether that
#   pathway's FIXED gene membership contains both significantly
#   up-regulated and significantly down-regulated genes at the same time
#   ("Mixed"), only one direction ("Up-only"/"Down-only"), or neither. If
#   a pathway's own member genes are internally split, whatever single
#   net direction GSEA/GSVA report for it (NES sign, diff_score sign) is
#   just whichever side has slightly more weight, not evidence of a
#   coherent, unified biological response -- a further, independent
#   reason (alongside the gene-level up/down consistency findings in
#   R/gene_rank_overlap.R) to treat pathway-level directionality calls
#   cautiously.
#
#   KEGG only: reuses build_kegg_gene_sets() from R/gsva.R as-is (the
#   only reusable gene-set-membership builder in this codebase; GSVA
#   itself is KEGG-only here for the same reason -- gseGO() builds the
#   full GO DAG internally via clusterProfiler's unexported machinery,
#   not something worth reimplementing just for this diagnostic).
#
################################################################################

#===============================================================================
# Run Gene Set Direction Mix Diagnostic
#===============================================================================

run_geneset_direction_mix <- function(config) {

    message("Running gene set direction mix diagnostic...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    deseq2_table_dir <- file.path(config$output_directory, "DESeq2", "Tables")
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

    shrunk_file_for <- function(group) {
        file.path(deseq2_table_dir, paste0("DESeq2_group_", group, "_vs_", config$reference_group, "_shrunk.csv"))
    }

    top_dose_table$has_deseq2 <- vapply(top_dose_table$group, function(g) file.exists(shrunk_file_for(g)), logical(1))
    eligible <- top_dose_table[top_dose_table$has_deseq2, ]
    compounds <- sort(eligible$compound)

    message(sprintf("Checking KEGG gene-set direction mix across %d compounds at their highest tested dose.", length(compounds)))

    #===========================================================================
    # Build the KEGG Gene Set Universe
    #===========================================================================

    org_pkg <- get_organism_db(config$organism)
    kegg_organism <- get_kegg_organism_code(config$organism)

    first_results <- utils::read.csv(shrunk_file_for(eligible$group[1]))
    universe_gene_ids <- first_results$gene

    gene_sets <- build_kegg_gene_sets(
        kegg_organism = kegg_organism,
        org_pkg = org_pkg,
        universe_gene_ids = universe_gene_ids
    )

    gene_sets <- Filter(
        function(s) length(s) >= config$gsva_min_gene_set_size && length(s) <= config$gsva_max_gene_set_size,
        gene_sets
    )

    message(sprintf("Checking direction mix across %d KEGG pathways.", length(gene_sets)))

    #===========================================================================
    # Classify Genes Up/Down Per Compound, Cross With Gene Set Membership
    #===========================================================================

    long_rows <- list()

    for (cmpd in compounds) {

        group <- eligible$group[eligible$compound == cmpd]
        results <- utils::read.csv(shrunk_file_for(group))

        classified <- classify_genes(results, "gene", "log2FoldChange", "padj", config$alpha, config$log2fc_cutoff)

        for (pathway_name in names(gene_sets)) {

            members <- gene_sets[[pathway_name]]

            n_up <- length(intersect(members, classified$up))
            n_down <- length(intersect(members, classified$down))

            category <- if (n_up > 0 && n_down > 0) {
                "Mixed"
            } else if (n_up > 0) {
                "Up-only"
            } else if (n_down > 0) {
                "Down-only"
            } else {
                "Neither"
            }

            long_rows[[length(long_rows) + 1]] <- data.frame(
                compound = cmpd,
                pathway = pathway_name,
                n_up = n_up,
                n_down = n_down,
                n_tested = length(members),
                category = category,
                pct_up_of_sig = if (n_up + n_down > 0) n_up / (n_up + n_down) else NA_real_,
                stringsAsFactors = FALSE
            )

        }

    }

    long_table <- do.call(rbind, long_rows)

    save_csv(
        long_table,
        file.path(diagnostic_dir, "geneset_direction_mix_long.csv"),
        row.names = FALSE
    )

    #===========================================================================
    # Per-Compound Summary: % of Pathways in Each Category
    #===========================================================================

    summary_rows <- lapply(compounds, function(cmpd) {

        sub <- long_table[long_table$compound == cmpd, ]
        n_total <- nrow(sub)

        data.frame(
            compound = cmpd,
            pct_mixed = mean(sub$category == "Mixed") * 100,
            pct_up_only = mean(sub$category == "Up-only") * 100,
            pct_down_only = mean(sub$category == "Down-only") * 100,
            pct_neither = mean(sub$category == "Neither") * 100,
            n_pathways = n_total
        )

    })

    summary_table <- do.call(rbind, summary_rows)
    summary_table <- summary_table[order(-summary_table$pct_mixed), ]

    save_csv(
        summary_table,
        file.path(diagnostic_dir, "geneset_direction_mix_summary_by_compound.csv"),
        row.names = FALSE
    )

    #===========================================================================
    # Plot 1: Per-Compound Stacked Bar Chart of Category %
    #===========================================================================

    stacked_df <- rbind(
        data.frame(compound = summary_table$compound, category = "Mixed", pct = summary_table$pct_mixed),
        data.frame(compound = summary_table$compound, category = "Up-only", pct = summary_table$pct_up_only),
        data.frame(compound = summary_table$compound, category = "Down-only", pct = summary_table$pct_down_only),
        data.frame(compound = summary_table$compound, category = "Neither", pct = summary_table$pct_neither)
    )

    stacked_df$compound <- factor(stacked_df$compound, levels = summary_table$compound)
    stacked_df$category <- factor(stacked_df$category, levels = c("Neither", "Down-only", "Mixed", "Up-only"))

    category_plot <- ggplot2::ggplot(stacked_df, ggplot2::aes(x = compound, y = pct, fill = category)) +
        ggplot2::geom_col(position = "stack", width = 0.7) +
        ggplot2::scale_fill_manual(values = c(
            "Neither" = "grey80", "Down-only" = "steelblue", "Mixed" = "mediumpurple", "Up-only" = "darkorange"
        )) +
        ggplot2::labs(
            title = "KEGG Gene Sets by Internal Up/Down Direction Mix, by Compound (Top Dose)",
            subtitle = "Ordered by % Mixed (pathways with BOTH significant up and down genes) descending",
            x = "Compound", y = "% of tested KEGG pathways", fill = "Category"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    if (isTRUE(config$save_png)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Geneset_direction_mix_by_compound.png"),
            plot = category_plot, width = max(config$figure_width, 11), height = config$figure_height, dpi = config$figure_dpi
        )
    }

    if (isTRUE(config$save_pdf)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Geneset_direction_mix_by_compound.pdf"),
            plot = category_plot, width = max(config$figure_width, 11), height = config$figure_height
        )
    }

    #===========================================================================
    # Plot 2: Pooled Histogram of the Internal Split (Mixed Pathways Only)
    #===========================================================================

    mixed_rows <- long_table[long_table$category == "Mixed", ]

    split_plot <- ggplot2::ggplot(mixed_rows, ggplot2::aes(x = pct_up_of_sig)) +
        ggplot2::geom_histogram(bins = 40, fill = "mediumpurple", color = NA) +
        ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey30") +
        ggplot2::labs(
            title = "Internal Up/Down Split Within \"Mixed\" KEGG Pathways",
            subtitle = sprintf(
                "%d (compound x pathway) mixed cells pooled across %d compounds; dashed = 50/50 split",
                nrow(mixed_rows), length(compounds)
            ),
            x = "Fraction of significant genes in the pathway that are up-regulated",
            y = "Count (compound x pathway cells)"
        ) +
        ggplot2::theme_bw()

    if (isTRUE(config$save_png)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Geneset_direction_mix_split_histogram.png"),
            plot = split_plot, width = config$figure_width, height = config$figure_height, dpi = config$figure_dpi
        )
    }

    if (isTRUE(config$save_pdf)) {
        ggplot2::ggsave(
            file.path(diagnostic_dir, "Geneset_direction_mix_split_histogram.pdf"),
            plot = split_plot, width = config$figure_width, height = config$figure_height
        )
    }

    message("Gene set direction mix diagnostic complete.")

    list(
        long_table = long_table,
        summary_by_compound = summary_table,
        category_plot = category_plot,
        split_plot = split_plot
    )

}
