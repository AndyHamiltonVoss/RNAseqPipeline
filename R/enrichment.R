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

    dot_plot <- enrichplot::dotplot(result, showCategory = n_show)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_base, "_dotplot.png")),
            plot = dot_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_base, "_dotplot.pdf")),
            plot = dot_plot,
            width = config$figure_width,
            height = config$figure_height
        )

    }

    invisible(TRUE)

}

#===============================================================================
# Run Full GO/KEGG Enrichment Across All Contrasts
#===============================================================================

run_enrichment <- function(config, metadata) {

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

    go_table_dir <- file.path(config$output_directory, "GO", "Tables")
    go_figure_dir <- file.path(config$output_directory, "GO", "Figures")
    kegg_table_dir <- file.path(config$output_directory, "KEGG", "Tables")
    kegg_figure_dir <- file.path(config$output_directory, "KEGG", "Figures")

    dir.create(go_table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(go_figure_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(kegg_table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(kegg_figure_dir, recursive = TRUE, showWarnings = FALSE)

    # All contrasts share one fitted DESeq2 model, so the tested gene universe
    # is identical across groups; map it to Entrez IDs once rather than per call.
    universe_entrez_ids <- NULL

    if (isTRUE(config$run_kegg)) {

        org_db <- get(org_pkg, envir = asNamespace(org_pkg))

        first_group_file <- file.path(
            config$output_directory, "DESeq2", "Tables",
            paste0("DESeq2_group_", groups[1], "_vs_", config$reference_group, "_shrunk.csv")
        )

        universe_gene_ids <- read.csv(first_group_file)$gene

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

        for (direction in c("up", "down")) {

            gene_list <- classify_genes(
                results_table, "gene", "log2FoldChange", "padj",
                config$alpha, config$log2fc_cutoff
            )[[direction]]

            message(
                "Processing enrichment: ", group, " (", direction, "regulated, ",
                length(gene_list), " genes)"
            )

            n_go_terms <- NA
            n_kegg_terms <- NA

            if (isTRUE(config$run_go)) {

                go_result <- run_go_enrichment(
                    gene_list, universe, org_pkg,
                    ont = config$go_ontology, alpha = config$alpha
                )

                save_enrichment_result(
                    go_result,
                    paste0("GO_", group, "_", direction),
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
                    paste0("KEGG_", group, "_", direction),
                    kegg_table_dir, kegg_figure_dir, config
                )

                n_kegg_terms <- if (!is.null(kegg_result)) nrow(as.data.frame(kegg_result)) else 0

            }

            summary_rows[[length(summary_rows) + 1]] <- data.frame(
                group = group,
                direction = direction,
                n_input_genes = length(gene_list),
                n_go_terms = n_go_terms,
                n_kegg_terms = n_kegg_terms,
                stringsAsFactors = FALSE
            )

        }

    }

    summary_table <- do.call(rbind, summary_rows)

    save_csv(
        summary_table,
        file.path(config$output_directory, "GO", "enrichment_summary.csv"),
        row.names = FALSE
    )

    message("GO/KEGG enrichment complete.")

    return(summary_table)

}
