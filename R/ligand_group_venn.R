################################################################################
# RNAseqPipeline
#
# File: ligand_group_venn.R
#
# Description:
#   Venn diagrams for specific, named groups of 2-5 compounds (e.g. a
#   pharmacologically related cluster), showing shared vs. unique DE genes
#   at each compound's own highest tested dose. Reuses the DESeq2-only
#   top-dose discovery logic already proven in R/compound_similarity.R
#   (NOT get_consensus_genes()/find_peak_dose() from R/compound_venn.R,
#   which require DESeq2 + edgeR + limma consensus and would silently
#   return empty gene sets for AMP0961, where edgeR/limma are disabled).
#   Venn geometry/rendering reuses VennDiagram::venn.diagram()'s native
#   list-of-sets interface (same call shape as build_method_venn() in
#   R/venn_diagrams.R) and save_venn_plot() for png/pdf output.
#
################################################################################

#===============================================================================
# Run Ligand-Group Venn Diagrams
#===============================================================================

run_ligand_group_venn <- function(config) {

    message("Running ligand-group Venn diagrams...")

    if (!requireNamespace("VennDiagram", quietly = TRUE)) {
        stop("Package 'VennDiagram' is required.")
    }

    if (length(config$venn_compound_groups) == 0) {
        message("No venn_compound_groups configured; nothing to do.")
        return(invisible(NULL))
    }

    table_dir <- file.path(config$output_directory, "DESeq2", "Tables")
    figure_dir <- file.path(config$output_directory, "Figures", "Venn")

    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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

    per_compound <- stats::setNames(
        lapply(top_dose_table$compound, function(cmpd) {

            group <- top_dose_table$group[top_dose_table$compound == cmpd]

            results <- read.csv(
                file.path(
                    table_dir,
                    paste0("DESeq2_group_", group, "_vs_", config$reference_group, "_shrunk.csv")
                )
            )

            classify_genes(
                results, "gene", "log2FoldChange", "padj",
                config$alpha, config$log2fc_cutoff
            )

        }),
        top_dose_table$compound
    )

    palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

    results_by_group <- list()

    group_consensus_sets <- list(up = list(), down = list())
    group_union_sets <- list(up = list(), down = list())

    for (group_name in names(config$venn_compound_groups)) {

        compounds <- config$venn_compound_groups[[group_name]]

        missing_compounds <- setdiff(compounds, names(per_compound))

        if (length(missing_compounds) > 0) {
            warning(
                "Skipping ligand group '", group_name, "': compound(s) not found ",
                "among discovered top-dose contrasts: ",
                paste(missing_compounds, collapse = ", ")
            )
            next
        }

        if (length(compounds) < 2 || length(compounds) > 5) {
            warning(
                "Skipping ligand group '", group_name, "': VennDiagram supports ",
                "2-5 sets, got ", length(compounds), "."
            )
            next
        }

        for (direction in c("up", "down")) {

            gene_sets <- stats::setNames(
                lapply(compounds, function(cmpd) per_compound[[cmpd]][[direction]]),
                compounds
            )

            if (!(group_name %in% config$venn_cross_group_exclude)) {
                group_consensus_sets[[direction]][[group_name]] <- Reduce(intersect, gene_sets)
                group_union_sets[[direction]][[group_name]] <- Reduce(union, gene_sets)
            }

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
                main = sprintf(
                    "%s: %sregulated genes (top dose)",
                    group_name, direction
                ),
                main.cex = 1.3,
                margin = 0.1
            )

            filename_base <- paste0("Venn_", group_name, "_", direction)

            save_venn_plot(
                venn_grob, filename_base, figure_dir, config,
                width = config$figure_width, height = config$figure_height
            )

            overlap_summary <- venn_regions_n(gene_sets)

            all_genes <- Reduce(union, gene_sets)

            gene_membership <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

            for (i in seq_along(gene_sets)) {
                gene_membership[[compounds[i]]] <- gene_membership$gene %in% gene_sets[[i]]
            }

            save_csv(
                overlap_summary,
                file.path(figure_dir, paste0(filename_base, "_overlap_counts.csv")),
                row.names = FALSE
            )

            save_csv(
                gene_membership,
                file.path(figure_dir, paste0(filename_base, "_gene_membership.csv")),
                row.names = FALSE
            )

            results_by_group[[paste0(group_name, "_", direction)]] <- list(
                gene_sets = gene_sets,
                overlap_summary = overlap_summary,
                gene_membership = gene_membership
            )

        }

    }

    #===========================================================================
    # Cross-Group Comparison: Each Group's Core Gene Set vs. Every Other Group
    #===========================================================================

    core_set_types <- list(consensus = group_consensus_sets, union = group_union_sets)

    core_type_labels <- c(
        consensus = "consensus: all compounds in group agree",
        union = "union: any compound in group"
    )

    for (core_type in names(core_set_types)) {

        for (direction in c("up", "down")) {

            group_gene_sets <- core_set_types[[core_type]][[direction]]

            if (length(group_gene_sets) < 2) {
                next
            }

            if (length(group_gene_sets) > 5) {
                warning(
                    "Skipping group-", core_type, " ", direction, " comparison: ",
                    "VennDiagram supports 2-5 sets, got ", length(group_gene_sets), "."
                )
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
                    "Group Core Comparison (%s): %sregulated genes",
                    core_type_labels[[core_type]], direction
                ),
                main.cex = 1.3,
                margin = 0.1
            )

            filename_base <- paste0("Venn_group_", core_type, "_", direction)

            save_venn_plot(
                venn_grob, filename_base, figure_dir, config,
                width = max(config$figure_width, 11), height = max(config$figure_height, 7)
            )

            overlap_summary <- venn_regions_n(group_gene_sets)

            all_genes <- Reduce(union, group_gene_sets)

            gene_membership <- data.frame(gene = all_genes, stringsAsFactors = FALSE)

            group_names_in_set <- names(group_gene_sets)

            for (i in seq_along(group_gene_sets)) {
                gene_membership[[group_names_in_set[i]]] <- gene_membership$gene %in% group_gene_sets[[i]]
            }

            save_csv(
                overlap_summary,
                file.path(figure_dir, paste0(filename_base, "_overlap_counts.csv")),
                row.names = FALSE
            )

            save_csv(
                gene_membership,
                file.path(figure_dir, paste0(filename_base, "_gene_membership.csv")),
                row.names = FALSE
            )

            results_by_group[[paste0("group_", core_type, "_", direction)]] <- list(
                gene_sets = group_gene_sets,
                overlap_summary = overlap_summary,
                gene_membership = gene_membership
            )

        }

    }

    message("Ligand-group Venn diagrams complete.")

    return(results_by_group)

}
