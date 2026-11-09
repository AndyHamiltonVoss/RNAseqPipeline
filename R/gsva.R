################################################################################
# RNAseqPipeline
#
# File: gsva.R
#
# Description:
#   Gene Set Variation Analysis (GSVA): per-sample KEGG pathway activity
#   scores from the VST-normalized expression matrix, plus limma differential
#   pathway activity per dose-vs-reference contrast.
#
################################################################################

#===============================================================================
# Build KEGG Gene Sets (Ensembl Gene IDs) for GSVA
#===============================================================================

build_kegg_gene_sets <- function(kegg_organism, org_pkg, universe_gene_ids) {

    gene_to_pathway <- KEGGREST::keggLink("pathway", kegg_organism)

    pathway_names <- KEGGREST::keggList("pathway", kegg_organism)
    pathway_names <- sub(" - [^-]*$", "", pathway_names)

    entrez_ids <- sub(paste0("^", kegg_organism, ":"), "", names(gene_to_pathway))
    pathway_ids <- sub("^path:", "", unname(gene_to_pathway))

    org_db <- get(org_pkg, envir = asNamespace(org_pkg))

    id_map <- suppressMessages(
        clusterProfiler::bitr(
            unique(entrez_ids), fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = org_db
        )
    )

    mapping <- merge(
        data.frame(entrez_id = entrez_ids, pathway_id = pathway_ids, stringsAsFactors = FALSE),
        id_map,
        by.x = "entrez_id", by.y = "ENTREZID"
    )

    mapping <- mapping[mapping$ENSEMBL %in% universe_gene_ids, ]

    gene_sets <- split(mapping$ENSEMBL, mapping$pathway_id)
    gene_sets <- lapply(gene_sets, unique)

    ids <- names(gene_sets)
    names(gene_sets) <- ifelse(ids %in% names(pathway_names), unname(pathway_names[ids]), ids)

    gene_sets

}

#===============================================================================
# Save GSVA Heatmap of Top-Variable Pathways
#===============================================================================

save_gsva_heatmap <- function(gsva_scores, metadata, config, figure_dir) {

    pathway_vars <- apply(gsva_scores, 1, stats::var)

    n_heatmap <- min(config$gsva_heatmap_top_pathways, nrow(gsva_scores))

    heatmap_pathways <- names(sort(pathway_vars, decreasing = TRUE))[seq_len(n_heatmap)]

    heatmap_matrix <- gsva_scores[heatmap_pathways, , drop = FALSE]

    annotation_col <- data.frame(
        Group = metadata[colnames(heatmap_matrix), "group"],
        row.names = colnames(heatmap_matrix)
    )

    heatmap_height <- max(config$figure_height, 3 + n_heatmap * 0.2)
    heatmap_width <- max(config$figure_width, ncol(heatmap_matrix) * 0.1) + 4

    if (isTRUE(config$save_png)) {

        pheatmap::pheatmap(
            heatmap_matrix,
            scale = "row",
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            annotation_col = annotation_col,
            show_colnames = FALSE,
            main = "Top Variable KEGG Pathway Activity (GSVA)",
            filename = file.path(figure_dir, "Heatmap_GSVA_KEGG_top_pathways.png"),
            width = heatmap_width,
            height = heatmap_height,
            res = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        pheatmap::pheatmap(
            heatmap_matrix,
            scale = "row",
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            annotation_col = annotation_col,
            show_colnames = FALSE,
            main = "Top Variable KEGG Pathway Activity (GSVA)",
            filename = file.path(figure_dir, "Heatmap_GSVA_KEGG_top_pathways.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

    }

    invisible(TRUE)

}

#===============================================================================
# Run limma on GSVA Scores (Differential Pathway Activity per Contrast)
#===============================================================================

run_gsva_limma <- function(gsva_scores, metadata, config, table_dir, figure_dir, label_pathways_by_group = list()) {

    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required.")
    }

    metadata <- relevel_group(metadata, config$reference_group)

    gsva_scores <- gsva_scores[, rownames(metadata), drop = FALSE]

    design <- stats::model.matrix(
        stats::as.formula(config$design_formula),
        data = metadata
    )

    fit <- limma::lmFit(gsva_scores, design)
    fit <- limma::eBayes(fit)

    coef_names <- colnames(design)

    contrast_names <- coef_names[grepl("^group", coef_names)]

    results_list <- list()

    for (coef_name in contrast_names) {

        contrast_label <- sub("^group", "", coef_name)

        res_df <- limma::topTable(
            fit,
            coef = coef_name,
            number = Inf,
            sort.by = "P"
        )

        res_df$pathway <- rownames(res_df)

        res_df <- res_df[
            ,
            c("pathway", setdiff(colnames(res_df), "pathway"))
        ]

        # limma always calls the coefficient "logFC"; here it's a difference
        # of GSVA enrichment scores, not a log fold change of expression.
        colnames(res_df)[colnames(res_df) == "logFC"] <- "diff_score"

        save_csv(
            res_df,
            file.path(
                table_dir,
                paste0(
                    "GSVA_KEGG_", contrast_label,
                    "_vs_", config$reference_group, "_limma.csv"
                )
            ),
            row.names = FALSE
        )

        save_gsva_volcano_figure(
            res_df, contrast_label, config$reference_group, config, figure_dir,
            label_pathways = label_pathways_by_group[[contrast_label]]
        )

        results_list[[coef_name]] <- res_df

    }

    results_list

}

#===============================================================================
# Pairwise Limma Comparison Between Two GSVA Groups (Direct, Not vs Reference)
#===============================================================================

run_gsva_pairwise_limma <- function(gsva_scores, metadata, group1, group2) {

    keep_samples <- rownames(metadata)[metadata$group %in% c(group1, group2)]

    pairwise_scores <- gsva_scores[, keep_samples, drop = FALSE]
    pairwise_metadata <- metadata[keep_samples, , drop = FALSE]
    pairwise_metadata$group <- factor(pairwise_metadata$group, levels = c(group2, group1))

    design <- stats::model.matrix(~ group, data = pairwise_metadata)

    fit <- limma::lmFit(pairwise_scores, design)
    fit <- limma::eBayes(fit)

    res_df <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")

    res_df$pathway <- rownames(res_df)
    res_df <- res_df[, c("pathway", setdiff(colnames(res_df), "pathway"))]

    # limma always calls the coefficient "logFC"; here it's a difference of
    # GSVA enrichment scores, not a log fold change of expression.
    colnames(res_df)[colnames(res_df) == "logFC"] <- "diff_score"

    res_df

}

#===============================================================================
# Build a GSVA Volcano Plot from a Limma topTable-style Result
#===============================================================================

build_gsva_volcano_plot <- function(result_df, group1, group2, alpha, label_pathways = NULL) {

    plot_data <- result_df
    plot_data$neg_log10_p <- -log10(plot_data$P.Value)
    plot_data$significant <- plot_data$adj.P.Val < alpha

    if (!is.null(label_pathways)) {
        label_data <- plot_data[plot_data$pathway %in% label_pathways, ]
    } else {
        n_label <- min(10, sum(plot_data$significant))
        label_data <- plot_data[order(plot_data$P.Value), ][seq_len(n_label), ]
    }

    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = diff_score, y = neg_log10_p, color = significant)
    ) +
        ggplot2::geom_point(size = 2, alpha = 0.7) +
        ggplot2::scale_color_manual(
            values = c(`TRUE` = "#d95f02", `FALSE` = "grey60"),
            labels = c(`TRUE` = paste0("padj < ", alpha), `FALSE` = "n.s."),
            name = NULL
        ) +
        ggrepel::geom_text_repel(
            data = label_data,
            ggplot2::aes(label = pathway),
            color = "black",
            size = 3,
            max.overlaps = Inf,
            box.padding = 0.6,
            point.padding = 0.3,
            force = 4,
            min.segment.length = 0,
            seed = 1234,
            show.legend = FALSE
        ) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::labs(
            title = paste0(
                "GSVA KEGG Pathway Activity: ",
                gsub("_", " ", group1), " vs. ", gsub("_", " ", group2)
            ),
            x = paste0("Difference in GSVA score (", group1, " − ", group2, ")"),
            y = expression(-log[10](p~value))
        ) +
        ggplot2::theme_bw()

}

#===============================================================================
# Save a GSVA Volcano Plot Figure (PNG/PDF)
#===============================================================================

save_gsva_volcano_figure <- function(result_df, group1, group2, config, figure_dir, alpha = config$alpha, label_pathways = NULL) {

    if (!requireNamespace("ggrepel", quietly = TRUE)) {
        stop("Package 'ggrepel' is required for GSVA volcano plots.")
    }

    comparison_label <- paste0(group1, "_vs_", group2)

    volcano_plot <- build_gsva_volcano_plot(result_df, group1, group2, alpha, label_pathways)

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0("Volcano_GSVA_KEGG_", comparison_label, ".png")),
            plot = volcano_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0("Volcano_GSVA_KEGG_", comparison_label, ".pdf")),
            plot = volcano_plot,
            width = config$figure_width,
            height = config$figure_height
        )

    }

    invisible(volcano_plot)

}

#===============================================================================
# Save GSVA Volcano Plot (+ Table) for a Pairwise Group Comparison
#===============================================================================

save_gsva_volcano <- function(pairwise_result, group1, group2, config, figure_dir, table_dir, alpha = config$alpha) {

    comparison_label <- paste0(group1, "_vs_", group2)

    save_csv(
        pairwise_result,
        file.path(table_dir, paste0("GSVA_KEGG_", comparison_label, "_limma.csv")),
        row.names = FALSE
    )

    save_gsva_volcano_figure(pairwise_result, group1, group2, config, figure_dir, alpha)

}

#===============================================================================
# Run Full GSVA (KEGG Gene Sets) with Differential Pathway Activity
#===============================================================================

run_gsva <- function(config, metadata, exclude_pattern = config$exclude_groups_pattern) {

    message("Running GSVA (KEGG gene sets)...")

    if (!requireNamespace("GSVA", quietly = TRUE)) {
        stop("Package 'GSVA' is required.")
    }

    dds_file <- file.path(config$output_directory, "DESeq2", "Objects", "dds_object.rds")

    if (!file.exists(dds_file)) {
        stop("DESeq2 dds object not found at ", dds_file, "; run DESeq2 first.")
    }

    dds <- readRDS(dds_file)

    vst_data <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))

    # Match the GO/KEGG ORA step: drop co-treatment/antagonist arms (e.g.
    # T007) from scoring and the downstream per-contrast comparisons. VST is
    # still computed from the full fit above, so normalization is unaffected.
    if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {

        excluded_groups <- unique(metadata$group[grepl(exclude_pattern, metadata$group)])

        if (length(excluded_groups) > 0) {
            message("Excluding groups matching '", exclude_pattern, "': ", paste(excluded_groups, collapse = ", "))
            metadata <- metadata[!metadata$group %in% excluded_groups, , drop = FALSE]
            vst_data <- vst_data[, rownames(metadata), drop = FALSE]
        }

    }

    org_pkg <- get_organism_db(config$organism)
    kegg_organism <- get_kegg_organism_code(config$organism)

    gene_sets <- build_kegg_gene_sets(
        kegg_organism = kegg_organism,
        org_pkg = org_pkg,
        universe_gene_ids = rownames(vst_data)
    )

    message("Scoring ", length(gene_sets), " candidate KEGG pathways across ", ncol(vst_data), " samples...")

    gsva_param <- GSVA::gsvaParam(
        exprData = vst_data,
        geneSets = gene_sets,
        kcdf = config$gsva_kcdf,
        minSize = config$gsva_min_gene_set_size,
        maxSize = config$gsva_max_gene_set_size
    )

    gsva_scores <- GSVA::gsva(gsva_param, verbose = FALSE)

    message(nrow(gsva_scores), " pathways passed the gene-set size filter.")

    output_dir <- file.path(config$output_directory, "GSVA")
    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")
    object_dir <- file.path(output_dir, "Objects")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        data.frame(pathway = rownames(gsva_scores), as.data.frame(gsva_scores)),
        file.path(table_dir, "GSVA_KEGG_scores.csv"),
        row.names = FALSE
    )

    saveRDS(gsva_scores, file.path(object_dir, "gsva_kegg_scores.rds"))

    # Run pairwise (group-vs-group) comparisons first: their significant
    # pathways are used to label the corresponding groups' vs-reference
    # volcano plots below, instead of each plot auto-labeling its own
    # top hits independently.
    pairwise_results <- list()
    label_pathways_by_group <- list()

    if (length(config$gsva_pairwise_comparisons) > 0) {

        for (comparison_name in names(config$gsva_pairwise_comparisons)) {

            comparison <- config$gsva_pairwise_comparisons[[comparison_name]]
            group1 <- comparison[1]
            group2 <- comparison[2]

            if (!all(c(group1, group2) %in% metadata$group)) {
                warning(
                    "Skipping GSVA pairwise comparison '", comparison_name,
                    "': group(s) not found in metadata (may have been excluded)."
                )
                next
            }

            pairwise_result <- run_gsva_pairwise_limma(gsva_scores, metadata, group1, group2)

            save_gsva_volcano(pairwise_result, group1, group2, config, figure_dir, table_dir)

            pairwise_results[[comparison_name]] <- pairwise_result

            significant_pathways <- pairwise_result$pathway[pairwise_result$adj.P.Val < config$alpha]

            label_pathways_by_group[[group1]] <- union(label_pathways_by_group[[group1]], significant_pathways)
            label_pathways_by_group[[group2]] <- union(label_pathways_by_group[[group2]], significant_pathways)

        }

    }

    diff_results <- run_gsva_limma(gsva_scores, metadata, config, table_dir, figure_dir, label_pathways_by_group)

    save_gsva_heatmap(gsva_scores, metadata, config, figure_dir)

    message("GSVA complete.")

    return(
        list(
            scores = gsva_scores,
            gene_sets = gene_sets,
            diff_results = diff_results,
            pairwise_results = pairwise_results
        )
    )

}
