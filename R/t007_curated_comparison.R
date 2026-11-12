################################################################################
# RNAseqPipeline
#
# File: t007_curated_comparison.R
#
# Description:
#   "T007-curated" upregulated gene sets: for each compound with a T007
#   co-treatment contrast, start from the genes upregulated by the
#   compound alone at its own top dose, then remove any gene that is also
#   upregulated when T007 is co-administered at that same dose. What's
#   left is the drug-alone-specific upregulated signature. Compounds
#   without a T007 co-treatment contrast (SR1664, SR16832, T007 itself)
#   are excluded -- curation is impossible without that data.
#
#   Reruns the compound-similarity (Jaccard + LFC correlation) and
#   ligand-group Venn analyses on these curated sets, writing to a
#   dedicated Output_AMP0961/T007_Curated/ folder so they're easy to find
#   and distinguish from the uncurated versions.
#
################################################################################

#===============================================================================
# Run T007-Curated Upregulated Gene Comparison
#===============================================================================

run_t007_curated_comparison <- function(config) {

    message("Running T007-curated upregulated gene comparison...")

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    if (!requireNamespace("VennDiagram", quietly = TRUE)) {
        stop("Package 'VennDiagram' is required.")
    }

    source_table_dir <- file.path(config$output_directory, "DESeq2", "Tables")

    output_dir <- file.path(config$output_directory, "T007_Curated")
    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")
    venn_dir <- file.path(figure_dir, "Venn")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)

    read_shrunk <- function(group) {
        read.csv(
            file.path(
                source_table_dir,
                paste0("DESeq2_group_", group, "_vs_", config$reference_group, "_shrunk.csv")
            )
        )
    }

    #===========================================================================
    # Discover Each Compound's Top Dose
    #===========================================================================

    shrunk_files <- list.files(
        source_table_dir,
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

    #===========================================================================
    # Determine Which Compounds Have a T007 Co-Treatment Contrast
    #===========================================================================

    combo_group_exists <- function(alone_group) {
        combo_group <- paste0(alone_group, "_", config$t007_co_treatment_suffix)
        file.exists(file.path(
            source_table_dir,
            paste0("DESeq2_group_", combo_group, "_vs_", config$reference_group, "_shrunk.csv")
        ))
    }

    top_dose_table$has_combo <- vapply(top_dose_table$group, combo_group_exists, logical(1))

    excluded_compounds <- top_dose_table$compound[!top_dose_table$has_combo]

    message(sprintf(
        "%d of %d compounds excluded (no T007 co-treatment contrast found): %s",
        length(excluded_compounds), nrow(top_dose_table),
        paste(excluded_compounds, collapse = ", ")
    ))

    eligible_table <- top_dose_table[top_dose_table$has_combo, ]

    if (nrow(eligible_table) == 0) {
        stop("No compounds have a T007 co-treatment contrast; nothing to curate.")
    }

    #===========================================================================
    # Build Curated Upregulated Gene Sets
    #===========================================================================

    curation <- lapply(seq_len(nrow(eligible_table)), function(i) {

        cmpd <- eligible_table$compound[i]
        alone_group <- eligible_table$group[i]
        combo_group <- paste0(alone_group, "_", config$t007_co_treatment_suffix)

        alone_results <- read_shrunk(alone_group)
        combo_results <- read_shrunk(combo_group)

        up_alone <- classify_genes(
            alone_results, "gene", "log2FoldChange", "padj",
            config$alpha, config$log2fc_cutoff
        )$up

        up_combo <- classify_genes(
            combo_results, "gene", "log2FoldChange", "padj",
            config$alpha, config$log2fc_cutoff
        )$up

        curated_up <- setdiff(up_alone, up_combo)

        list(
            curated_up = curated_up,
            lfc = stats::setNames(alone_results$log2FoldChange, alone_results$gene),
            n_up_alone = length(up_alone),
            n_up_with_t007 = length(up_combo),
            n_removed = length(intersect(up_alone, up_combo)),
            n_curated = length(curated_up)
        )

    })

    names(curation) <- eligible_table$compound

    compounds <- eligible_table$compound

    curated_up_sets <- stats::setNames(
        lapply(compounds, function(cmpd) curation[[cmpd]]$curated_up),
        compounds
    )

    #===========================================================================
    # Save Curation Tables
    #===========================================================================

    gene_rows <- lapply(compounds, function(cmpd) {
        genes <- curated_up_sets[[cmpd]]
        if (length(genes) == 0) {
            return(NULL)
        }
        data.frame(compound = cmpd, gene = genes, stringsAsFactors = FALSE)
    })

    curated_up_genes <- do.call(rbind, gene_rows)

    save_csv(
        curated_up_genes,
        file.path(table_dir, "curated_up_genes.csv"),
        row.names = FALSE
    )

    curated_up_gene_counts <- data.frame(
        compound = compounds,
        n_up_alone = vapply(compounds, function(cmpd) curation[[cmpd]]$n_up_alone, numeric(1)),
        n_up_with_t007 = vapply(compounds, function(cmpd) curation[[cmpd]]$n_up_with_t007, numeric(1)),
        n_removed = vapply(compounds, function(cmpd) curation[[cmpd]]$n_removed, numeric(1)),
        n_curated = vapply(compounds, function(cmpd) curation[[cmpd]]$n_curated, numeric(1)),
        stringsAsFactors = FALSE
    )

    save_csv(
        curated_up_gene_counts,
        file.path(table_dir, "curated_up_gene_counts.csv"),
        row.names = FALSE
    )

    #===========================================================================
    # Jaccard Similarity of Curated Upregulated Gene Sets
    #===========================================================================

    jaccard_matrix <- matrix(
        NA_real_,
        nrow = length(compounds), ncol = length(compounds),
        dimnames = list(compounds, compounds)
    )

    for (i in compounds) {
        for (j in compounds) {

            union_size <- length(union(curated_up_sets[[i]], curated_up_sets[[j]]))

            jaccard_matrix[i, j] <- if (union_size == 0) {
                NA_real_
            } else {
                length(intersect(curated_up_sets[[i]], curated_up_sets[[j]])) / union_size
            }

        }
    }

    save_csv(jaccard_matrix, file.path(table_dir, "compound_jaccard_curated.csv"))

    #===========================================================================
    # Pearson Correlation of log2FoldChange, Restricted to the Curated Gene Panel
    #===========================================================================

    curated_gene_panel <- Reduce(union, curated_up_sets)

    lfc_matrix <- Reduce(
        function(a, b) merge(a, b, by = "gene", all = TRUE),
        lapply(compounds, function(cmpd) {
            df <- data.frame(
                gene = names(curation[[cmpd]]$lfc),
                stringsAsFactors = FALSE
            )
            df[[cmpd]] <- as.numeric(curation[[cmpd]]$lfc)
            df
        })
    )

    rownames(lfc_matrix) <- lfc_matrix$gene
    lfc_matrix$gene <- NULL
    lfc_matrix <- lfc_matrix[rownames(lfc_matrix) %in% curated_gene_panel, compounds, drop = FALSE]

    correlation_matrix <- stats::cor(
        as.matrix(lfc_matrix),
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    save_csv(correlation_matrix, file.path(table_dir, "compound_lfc_correlation_curated.csv"))

    #===========================================================================
    # Save Heatmaps
    #===========================================================================

    heatmap_height <- max(config$figure_height, length(compounds) * 0.5)
    heatmap_width <- heatmap_height + 2

    if (isTRUE(config$save_png)) {

        pheatmap::pheatmap(
            jaccard_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "T007-Curated Compound Similarity: Jaccard Overlap (Upregulated, Top Dose)",
            filename = file.path(figure_dir, "Compound_jaccard_curated_heatmap.png"),
            width = heatmap_width,
            height = heatmap_height,
            res = config$figure_dpi
        )

        pheatmap::pheatmap(
            correlation_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "T007-Curated Compound Similarity: log2FoldChange Correlation (Curated Gene Panel)",
            filename = file.path(figure_dir, "Compound_LFC_correlation_curated_heatmap.png"),
            width = heatmap_width,
            height = heatmap_height,
            res = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        pheatmap::pheatmap(
            jaccard_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "T007-Curated Compound Similarity: Jaccard Overlap (Upregulated, Top Dose)",
            filename = file.path(figure_dir, "Compound_jaccard_curated_heatmap.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

        pheatmap::pheatmap(
            correlation_matrix,
            cluster_rows = isTRUE(config$cluster_rows),
            cluster_cols = isTRUE(config$cluster_columns),
            display_numbers = TRUE,
            number_format = "%.2f",
            main = "T007-Curated Compound Similarity: log2FoldChange Correlation (Curated Gene Panel)",
            filename = file.path(figure_dir, "Compound_LFC_correlation_curated_heatmap.pdf"),
            width = heatmap_width,
            height = heatmap_height
        )

    }

    #===========================================================================
    # Within-Group Venn Diagrams (Curated Upregulated Genes)
    #===========================================================================

    palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

    results <- list(
        curated_up_sets = curated_up_sets,
        curated_up_gene_counts = curated_up_gene_counts,
        jaccard_matrix = jaccard_matrix,
        correlation_matrix = correlation_matrix
    )

    group_consensus_sets <- list()
    group_union_sets <- list()

    for (group_name in names(config$venn_compound_groups)) {

        group_compounds <- config$venn_compound_groups[[group_name]]

        missing_compounds <- setdiff(group_compounds, compounds)

        if (length(missing_compounds) > 0) {
            warning(
                "Skipping group '", group_name, "' for T007-curated comparison: ",
                "compound(s) not eligible (no T007 co-treatment contrast, or not found): ",
                paste(missing_compounds, collapse = ", ")
            )
            next
        }

        if (length(group_compounds) < 2 || length(group_compounds) > 5) {
            next
        }

        gene_sets <- curated_up_sets[group_compounds]

        group_consensus_sets[[group_name]] <- Reduce(intersect, gene_sets)
        group_union_sets[[group_name]] <- Reduce(union, gene_sets)

        futile.logger::flog.threshold(
            futile.logger::ERROR,
            name = "VennDiagramLogger"
        )

        venn_grob <- VennDiagram::venn.diagram(
            x = gene_sets,
            filename = NULL,
            fill = palette[seq_along(gene_sets)],
            alpha = 0.5,
            cex = 1.1,
            cat.cex = 1.1,
            main = sprintf("%s: T007-curated upregulated genes (top dose)", group_name),
            main.cex = 1.3,
            margin = 0.1
        )

        filename_base <- paste0("Venn_", group_name, "_curated_up")

        save_venn_plot(
            venn_grob, filename_base, venn_dir, config,
            width = config$figure_width, height = config$figure_height
        )

        overlap_summary <- venn_regions_n(gene_sets)

        all_genes <- Reduce(union, gene_sets)

        gene_membership <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

        for (k in seq_along(gene_sets)) {
            gene_membership[[group_compounds[k]]] <- gene_membership$gene %in% gene_sets[[k]]
        }

        save_csv(
            overlap_summary,
            file.path(venn_dir, paste0(filename_base, "_overlap_counts.csv")),
            row.names = FALSE
        )

        save_csv(
            gene_membership,
            file.path(venn_dir, paste0(filename_base, "_gene_membership.csv")),
            row.names = FALSE
        )

        results[[paste0(group_name, "_curated_up")]] <- list(
            gene_sets = gene_sets,
            overlap_summary = overlap_summary,
            gene_membership = gene_membership
        )

    }

    #===========================================================================
    # Cross-Group Comparison (Curated Upregulated Genes)
    #===========================================================================

    core_set_types <- list(consensus = group_consensus_sets, union = group_union_sets)

    for (core_type in names(core_set_types)) {

        group_gene_sets <- core_set_types[[core_type]]

        if (length(group_gene_sets) < 2 || length(group_gene_sets) > 5) {
            next
        }

        futile.logger::flog.threshold(
            futile.logger::ERROR,
            name = "VennDiagramLogger"
        )

        venn_grob <- VennDiagram::venn.diagram(
            x = group_gene_sets,
            filename = NULL,
            fill = palette[seq_along(group_gene_sets)],
            alpha = 0.5,
            cex = 1.1,
            cat.cex = 1.1,
            main = sprintf(
                "Group Core Comparison (%s): T007-curated upregulated genes",
                core_type
            ),
            main.cex = 1.3,
            margin = 0.1
        )

        filename_base <- paste0("Venn_group_", core_type, "_curated_up")

        save_venn_plot(
            venn_grob, filename_base, venn_dir, config,
            width = max(config$figure_width, 11), height = max(config$figure_height, 7)
        )

        overlap_summary <- venn_regions_n(group_gene_sets)

        all_genes <- Reduce(union, group_gene_sets)

        gene_membership <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

        group_names_in_set <- names(group_gene_sets)

        for (k in seq_along(group_gene_sets)) {
            gene_membership[[group_names_in_set[k]]] <- gene_membership$gene %in% group_gene_sets[[k]]
        }

        save_csv(
            overlap_summary,
            file.path(venn_dir, paste0(filename_base, "_overlap_counts.csv")),
            row.names = FALSE
        )

        save_csv(
            gene_membership,
            file.path(venn_dir, paste0(filename_base, "_gene_membership.csv")),
            row.names = FALSE
        )

        results[[paste0("group_", core_type, "_curated_up")]] <- list(
            gene_sets = group_gene_sets,
            overlap_summary = overlap_summary,
            gene_membership = gene_membership
        )

    }

    message("T007-curated upregulated gene comparison complete.")

    return(results)

}
