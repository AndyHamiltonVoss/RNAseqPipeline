################################################################################
# RNAseqPipeline
#
# File: enrichment.R
#
# Description:
#   GO and KEGG over-representation analysis on DESeq2 significant genes,
#   per contrast, using clusterProfiler.
#
################################################################################

#===============================================================================
# Resolve Organism Annotation Package / KEGG Organism Code
#===============================================================================

get_organism_db <- function(organism) {

    org_pkg <- switch(
        tolower(organism),
        human = "org.Hs.eg.db",
        mouse = "org.Mm.eg.db",
        stop("Unsupported organism for enrichment: ", organism)
    )

    if (!requireNamespace(org_pkg, quietly = TRUE)) {
        stop(
            "Package '", org_pkg, "' is required for ", organism,
            " GO/KEGG enrichment but is not installed."
        )
    }

    return(org_pkg)

}

get_kegg_organism_code <- function(organism) {

    switch(
        tolower(organism),
        human = "hsa",
        mouse = "mmu",
        stop("Unsupported organism for KEGG enrichment: ", organism)
    )

}

#===============================================================================
# Run GO Over-Representation Analysis
#===============================================================================

run_go_enrichment <- function(gene_list, universe, org_pkg, ont = "BP", alpha = 0.05) {

    if (length(gene_list) == 0) {
        return(NULL)
    }

    org_db <- get(org_pkg, envir = asNamespace(org_pkg))

    result <- tryCatch(
        clusterProfiler::enrichGO(
            gene = gene_list,
            universe = universe,
            OrgDb = org_db,
            keyType = "ENSEMBL",
            ont = ont,
            pAdjustMethod = "BH",
            pvalueCutoff = alpha,
            readable = TRUE
        ),
        error = function(e) {
            warning("GO enrichment failed: ", conditionMessage(e))
            NULL
        }
    )

    return(result)

}

#===============================================================================
# Run KEGG Over-Representation Analysis
#===============================================================================

run_kegg_enrichment <- function(gene_list, universe_entrez_ids, org_pkg, kegg_organism, alpha = 0.05) {

    if (length(gene_list) == 0) {
        return(NULL)
    }

    org_db <- get(org_pkg, envir = asNamespace(org_pkg))

    gene_ids <- suppressMessages(
        clusterProfiler::bitr(
            gene_list,
            fromType = "ENSEMBL",
            toType = "ENTREZID",
            OrgDb = org_db
        )
    )

    if (nrow(gene_ids) == 0) {
        return(NULL)
    }

    result <- tryCatch(
        clusterProfiler::enrichKEGG(
            gene = gene_ids$ENTREZID,
            universe = universe_entrez_ids,
            organism = kegg_organism,
            pAdjustMethod = "BH",
            pvalueCutoff = alpha
        ),
        error = function(e) {
            warning("KEGG enrichment failed: ", conditionMessage(e))
            NULL
        }
    )

    return(result)

}

#===============================================================================
# Save Enrichment Result (Table + Dotplot)
#===============================================================================

save_enrichment_result <- function(result, filename_base, table_dir, figure_dir, config) {

    if (is.null(result)) {
        return(invisible(FALSE))
    }

    result_df <- as.data.frame(result)

    save_csv(
        result_df,
        file.path(table_dir, paste0(filename_base, ".csv")),
        row.names = FALSE
    )

    if (nrow(result_df) == 0) {
        return(invisible(FALSE))
    }

    n_show <- min(20, nrow(result_df))

    plot_title <- gsub("_", " ", filename_base)

    dot_plot <- enrichplot::dotplot(result, showCategory = n_show) +
        ggplot2::ggtitle(plot_title)

    # Long GO/KEGG term names wrap onto multiple lines, so a fixed height
    # crowds/overlaps y-axis labels once many categories are shown; scale
    # with n_show instead.
    plot_height <- max(config$figure_height, n_show * 0.5)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_base, "_dotplot.png")),
            plot = dot_plot,
            width = config$figure_width,
            height = plot_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_base, "_dotplot.pdf")),
            plot = dot_plot,
            width = config$figure_width,
            height = plot_height
        )

    }

    invisible(TRUE)

}

#===============================================================================
# Plot Dose-Response Curves of Enrichment Input Gene Counts
#===============================================================================

plot_enrichment_gene_counts <- function(summary_table, compounds) {

    parsed <- lapply(summary_table$group, parse_dose_group)

    keep <- vapply(parsed, Negate(is.null), logical(1))

    plot_data <- summary_table[keep, ]
    plot_data$compound <- vapply(parsed[keep], function(x) x$compound, character(1))
    plot_data$dose_nM <- vapply(parsed[keep], function(x) x$dose_nM, numeric(1))

    plot_data <- plot_data[plot_data$compound %in% compounds, ]

    if (nrow(plot_data) == 0) {
        return(NULL)
    }

    plot_data$direction <- factor(plot_data$direction, levels = c("up", "down", "combined"))

    direction_colors <- c(
        up = "#d95f02",
        down = "#1b9e77",
        combined = "#7570b3"
    )

    plot <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = dose_nM,
            y = n_input_genes,
            color = direction,
            group = direction
        )
    ) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2) +
        ggplot2::facet_wrap(~compound, scales = "free_x") +
        ggplot2::scale_x_log10() +
        ggplot2::scale_color_manual(values = direction_colors) +
        ggplot2::labs(
            title = "GO/KEGG Input Gene Counts vs. Dose",
            x = "Dose (nM, log10 scale)",
            y = "Number of genes input to enrichment",
            color = "Direction"
        ) +
        ggplot2::theme_bw()

    return(plot)

}

#===============================================================================
# Classify Up/Down Genes for a Set of Sample Groups
#===============================================================================

classify_genes_for_groups <- function(group_names, output_directory, reference_group, alpha, log2fc_cutoff) {

    classified_list <- lapply(group_names, function(group) {

        deseq2_file <- file.path(
            output_directory, "DESeq2", "Tables",
            paste0("DESeq2_group_", group, "_vs_", reference_group, "_shrunk.csv")
        )

        if (!file.exists(deseq2_file)) {
            warning("No DESeq2 results found for group: ", group, "; skipping.")
            return(NULL)
        }

        classify_genes(
            read.csv(deseq2_file), "gene", "log2FoldChange", "padj",
            alpha, log2fc_cutoff
        )

    })

    Filter(Negate(is.null), classified_list)

}

#===============================================================================
# Filter a Compound's Dose Series Down to a Set of Dose Strings
#===============================================================================

filter_dose_series <- function(groups, compound, doses) {

    series <- get_dose_series(groups, compound)

    if (!is.null(doses)) {
        dose_nM_values <- vapply(doses, dose_string_to_nM, numeric(1))
        series <- series[series$dose_nM %in% dose_nM_values, ]
    }

    series

}

#===============================================================================
# Build Dose-Consistent "Core" Gene Sets for a Compound
#===============================================================================

build_core_gene_sets <- function(
    compound,
    groups,
    output_directory,
    reference_group,
    alpha,
    log2fc_cutoff,
    doses = NULL,
    exclude_doses = NULL
) {

    series <- filter_dose_series(groups, compound, doses)

    if (nrow(series) == 0) {
        warning("No dose series found for compound: ", compound, "; skipping core gene set.")
        return(NULL)
    }

    classified_list <- classify_genes_for_groups(
        series$group, output_directory, reference_group, alpha, log2fc_cutoff
    )

    if (length(classified_list) == 0) {
        return(NULL)
    }

    core_up <- Reduce(intersect, lapply(classified_list, function(x) x$up))
    core_down <- Reduce(intersect, lapply(classified_list, function(x) x$down))

    # Subtract genes already up/down at the excluded (e.g. sub-saturating)
    # doses, isolating genes newly recruited only once those doses are cleared.
    if (!is.null(exclude_doses)) {

        exclude_series <- filter_dose_series(groups, compound, exclude_doses)

        exclude_classified <- classify_genes_for_groups(
            exclude_series$group, output_directory, reference_group, alpha, log2fc_cutoff
        )

        exclude_up <- if (length(exclude_classified) > 0)
            Reduce(union, lapply(exclude_classified, function(x) x$up))
        else character(0)

        exclude_down <- if (length(exclude_classified) > 0)
            Reduce(union, lapply(exclude_classified, function(x) x$down))
        else character(0)

        core_up <- setdiff(core_up, exclude_up)
        core_down <- setdiff(core_down, exclude_down)

    }

    list(
        up = core_up,
        down = core_down,
        combined = union(core_up, core_down),
        n_doses = length(classified_list)
    )

}

#===============================================================================
# Run GO/KEGG Enrichment for One Named Gene-Set Group (Up/Down/Combined)
#===============================================================================

process_enrichment_group <- function(
    group_label,
    gene_lists,
    universe_gene_ids,
    universe_entrez_ids,
    org_pkg,
    kegg_organism,
    config,
    go_table_dir,
    go_figure_dir,
    kegg_table_dir,
    kegg_figure_dir,
    message_suffix = ""
) {

    rows <- list()

    for (direction in c("up", "down", "combined")) {

        gene_list <- gene_lists[[direction]]

        message(
            "Processing enrichment: ", group_label, " (", direction, "regulated, ",
            length(gene_list), " genes", message_suffix, ")"
        )

        n_go_terms <- NA
        n_kegg_terms <- NA

        if (isTRUE(config$run_go)) {

            go_result <- run_go_enrichment(
                gene_list, universe_gene_ids, org_pkg,
                ont = config$go_ontology, alpha = config$alpha
            )

            save_enrichment_result(
                go_result,
                paste0("GO_", group_label, "_", direction),
                go_table_dir, go_figure_dir, config
            )

            n_go_terms <- if (!is.null(go_result)) nrow(as.data.frame(go_result)) else 0

        }

        if (isTRUE(config$run_kegg)) {

            kegg_result <- run_kegg_enrichment(
                gene_list, universe_entrez_ids, org_pkg, kegg_organism, alpha = config$alpha
            )

            save_enrichment_result(
                kegg_result,
                paste0("KEGG_", group_label, "_", direction),
                kegg_table_dir, kegg_figure_dir, config
            )

            n_kegg_terms <- if (!is.null(kegg_result)) nrow(as.data.frame(kegg_result)) else 0

        }

        rows[[direction]] <- data.frame(
            group = group_label,
            direction = direction,
            n_input_genes = length(gene_list),
            n_go_terms = n_go_terms,
            n_kegg_terms = n_kegg_terms,
            stringsAsFactors = FALSE
        )

    }

    do.call(rbind, rows)

}

#===============================================================================
# Run Full GO/KEGG Enrichment Across All Contrasts
#===============================================================================

run_enrichment <- function(
    config,
    metadata,
    log2fc_cutoff = config$go_kegg_log2fc_cutoff,
    exclude_pattern = config$exclude_groups_pattern
) {

    message("Running GO/KEGG enrichment...")

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required.")
    }

    org_pkg <- get_organism_db(config$organism)

    kegg_organism <- if (isTRUE(config$run_kegg)) get_kegg_organism_code(config$organism) else NA

    groups <- setdiff(
        levels(factor(metadata$group)),
        config$reference_group
    )

    if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {
        excluded <- groups[grepl(exclude_pattern, groups)]
        if (length(excluded) > 0) {
            message("Excluding groups matching '", exclude_pattern, "': ", paste(excluded, collapse = ", "))
        }
        groups <- groups[!grepl(exclude_pattern, groups)]
    }

    go_table_dir <- file.path(config$output_directory, "GO", "Tables")
    go_figure_dir <- file.path(config$output_directory, "GO", "Figures")
    kegg_table_dir <- file.path(config$output_directory, "KEGG", "Tables")
    kegg_figure_dir <- file.path(config$output_directory, "KEGG", "Figures")

    dir.create(go_table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(go_figure_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(kegg_table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(kegg_figure_dir, recursive = TRUE, showWarnings = FALSE)

    # All contrasts share one fitted DESeq2 model, so the tested gene universe
    # is identical across groups; read it once rather than per call.
    first_group_file <- file.path(
        config$output_directory, "DESeq2", "Tables",
        paste0("DESeq2_group_", groups[1], "_vs_", config$reference_group, "_shrunk.csv")
    )

    universe_gene_ids <- read.csv(first_group_file)$gene

    universe_entrez_ids <- NULL

    if (isTRUE(config$run_kegg)) {

        org_db <- get(org_pkg, envir = asNamespace(org_pkg))

        universe_entrez_ids <- suppressMessages(
            clusterProfiler::bitr(
                universe_gene_ids,
                fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = org_db
            )
        )$ENTREZID

    }

    summary_rows <- list()

    for (group in groups) {

        deseq2_file <- file.path(
            config$output_directory, "DESeq2", "Tables",
            paste0("DESeq2_group_", group, "_vs_", config$reference_group, "_shrunk.csv")
        )

        if (!file.exists(deseq2_file)) {
            warning("No DESeq2 results found for group: ", group, "; skipping enrichment.")
            next
        }

        results_table <- read.csv(deseq2_file)

        universe <- results_table$gene

        classified <- classify_genes(
            results_table, "gene", "log2FoldChange", "padj",
            config$alpha, log2fc_cutoff
        )

        gene_lists <- list(
            up = classified$up,
            down = classified$down,
            combined = union(classified$up, classified$down)
        )

        summary_rows[[length(summary_rows) + 1]] <- process_enrichment_group(
            group_label = group,
            gene_lists = gene_lists,
            universe_gene_ids = universe,
            universe_entrez_ids = universe_entrez_ids,
            org_pkg = org_pkg,
            kegg_organism = kegg_organism,
            config = config,
            go_table_dir = go_table_dir,
            go_figure_dir = go_figure_dir,
            kegg_table_dir = kegg_table_dir,
            kegg_figure_dir = kegg_figure_dir
        )

    }

    # Per-dose "up"/"down" gene sets can churn substantially between doses even
    # when their sizes stay flat, which makes per-dose enrichment noisy. Also
    # run enrichment on the "core" set of genes consistently up/down across
    # every dose of a compound, for a more stable pathway signature.
    for (compound in config$dose_response_compounds) {

        core_sets <- build_core_gene_sets(
            compound = compound,
            groups = groups,
            output_directory = config$output_directory,
            reference_group = config$reference_group,
            alpha = config$alpha,
            log2fc_cutoff = log2fc_cutoff
        )

        if (is.null(core_sets)) next

        summary_rows[[length(summary_rows) + 1]] <- process_enrichment_group(
            group_label = paste0(compound, "_core"),
            gene_lists = core_sets,
            universe_gene_ids = universe_gene_ids,
            universe_entrez_ids = universe_entrez_ids,
            org_pkg = org_pkg,
            kegg_organism = kegg_organism,
            config = config,
            go_table_dir = go_table_dir,
            go_figure_dir = go_figure_dir,
            kegg_table_dir = kegg_table_dir,
            kegg_figure_dir = kegg_figure_dir,
            message_suffix = paste0(", consistent across ", core_sets$n_doses, " doses")
        )

    }

    # Custom dose-subset core sets (e.g. only the doses above a compound's
    # receptor binding affinity), configured explicitly since they don't apply
    # uniformly to every compound.
    for (subset_name in names(config$core_gene_set_subsets)) {

        subset_spec <- config$core_gene_set_subsets[[subset_name]]

        core_sets <- build_core_gene_sets(
            compound = subset_spec$compound,
            groups = groups,
            output_directory = config$output_directory,
            reference_group = config$reference_group,
            alpha = config$alpha,
            log2fc_cutoff = log2fc_cutoff,
            doses = subset_spec$doses,
            exclude_doses = subset_spec$exclude_doses
        )

        if (is.null(core_sets)) next

        exclude_note <- if (!is.null(subset_spec$exclude_doses))
            paste0("; excluding genes also seen at ", paste(subset_spec$exclude_doses, collapse = "/"))
        else ""

        summary_rows[[length(summary_rows) + 1]] <- process_enrichment_group(
            group_label = subset_name,
            gene_lists = core_sets,
            universe_gene_ids = universe_gene_ids,
            universe_entrez_ids = universe_entrez_ids,
            org_pkg = org_pkg,
            kegg_organism = kegg_organism,
            config = config,
            go_table_dir = go_table_dir,
            go_figure_dir = go_figure_dir,
            kegg_table_dir = kegg_table_dir,
            kegg_figure_dir = kegg_figure_dir,
            message_suffix = paste0(
                ", consistent across ", core_sets$n_doses, " doses", exclude_note
            )
        )

    }

    summary_table <- do.call(rbind, summary_rows)

    save_csv(
        summary_table,
        file.path(config$output_directory, "GO", "enrichment_summary.csv"),
        row.names = FALSE
    )

    dose_plot <- plot_enrichment_gene_counts(
        summary_table,
        compounds = config$dose_response_compounds
    )

    if (!is.null(dose_plot)) {

        figures_dir <- file.path(config$output_directory, "Figures")
        dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

        if (isTRUE(config$save_png)) {

            ggplot2::ggsave(
                file.path(figures_dir, "enrichment_gene_counts_dose_response.png"),
                plot = dose_plot,
                width = config$figure_width,
                height = config$figure_height,
                dpi = config$figure_dpi
            )

        }

        if (isTRUE(config$save_pdf)) {

            ggplot2::ggsave(
                file.path(figures_dir, "enrichment_gene_counts_dose_response.pdf"),
                plot = dose_plot,
                width = config$figure_width,
                height = config$figure_height
            )

        }

    }

    message("GO/KEGG enrichment complete.")

    return(summary_table)

}
