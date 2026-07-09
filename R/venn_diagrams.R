################################################################################
# RNAseqPipeline
#
# File: venn_diagrams.R
#
# Description:
#   Cross-method Venn diagrams showing overlap of up/downregulated genes
#   called by DESeq2, edgeR, and limma-voom for each dose-response contrast.
#   Reuses classify_genes() / load_method_results() from dose_response.R.
#
################################################################################

#===============================================================================
# Build a Three-Way Venn Diagram Grob
#===============================================================================

build_method_venn <- function(gene_sets, title) {

    futile.logger::flog.threshold(
        futile.logger::ERROR,
        name = "VennDiagramLogger"
    )

    grob <- VennDiagram::venn.diagram(
        x = gene_sets,
        filename = NULL,
        fill = c("#1b9e77", "#d95f02", "#7570b3"),
        alpha = 0.5,
        cex = 1.2,
        cat.cex = 1.1,
        main = title,
        main.cex = 1.3,
        margin = 0.1
    )

    return(grob)

}

#===============================================================================
# Save Venn Diagram Grob
#===============================================================================

save_venn_plot <- function(
    venn_grob,
    filename_base,
    figure_dir,
    config,
    width = config$figure_width,
    height = config$figure_height
) {

    if (isTRUE(config$save_png)) {

        png(
            filename = file.path(figure_dir, paste0(filename_base, ".png")),
            width = width,
            height = height,
            units = "in",
            res = config$figure_dpi
        )

        grid::grid.newpage()
        grid::grid.draw(venn_grob)

        dev.off()
    }

    if (isTRUE(config$save_pdf)) {

        pdf(
            file = file.path(figure_dir, paste0(filename_base, ".pdf")),
            width = width,
            height = height
        )

        grid::grid.newpage()
        grid::grid.draw(venn_grob)

        dev.off()
    }

    invisible(TRUE)
}

#===============================================================================
# Break Three Gene Sets Into Non-Overlapping Venn Regions
#===============================================================================

venn_regions <- function(gene_sets) {

    method_names <- names(gene_sets)

    a <- gene_sets[[1]]
    b <- gene_sets[[2]]
    c <- gene_sets[[3]]

    regions <- list()

    regions[[paste(method_names[1], "only")]] <- setdiff(a, union(b, c))
    regions[[paste(method_names[2], "only")]] <- setdiff(b, union(a, c))
    regions[[paste(method_names[3], "only")]] <- setdiff(c, union(a, b))
    regions[[paste(method_names[1], "&", method_names[2])]] <- setdiff(intersect(a, b), c)
    regions[[paste(method_names[1], "&", method_names[3])]] <- setdiff(intersect(a, c), b)
    regions[[paste(method_names[2], "&", method_names[3])]] <- setdiff(intersect(b, c), a)
    regions[["All three"]] <- intersect(intersect(a, b), c)

    return(regions)

}

#===============================================================================
# Run Full Venn Diagram Comparison
#===============================================================================

run_venn_diagrams <- function(config, metadata) {

    message("Running cross-method Venn diagram comparison...")

    if (!requireNamespace("VennDiagram", quietly = TRUE)) {
        stop("Package 'VennDiagram' is required.")
    }

    if (!isTRUE(config$run_deseq2) || !isTRUE(config$run_edger) || !isTRUE(config$run_limma)) {
        warning(
            "Venn diagram comparison requires DESeq2, edgeR, and limma to all be ",
            "enabled; results from disabled methods are treated as empty gene sets."
        )
    }

    groups <- setdiff(
        levels(factor(metadata$group)),
        config$reference_group
    )

    figure_dir <- file.path(config$output_directory, "Figures", "Venn")

    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    contrast_groups <- unique(
        unlist(
            lapply(
                config$dose_response_compounds,
                function(compound) get_dose_series(groups, compound)$group
            )
        )
    )

    overlap_counts <- list()
    gene_membership <- list()

    for (group in contrast_groups) {

        results <- load_method_results(group, config$output_directory, config$reference_group)

        empty <- list(up = character(0), down = character(0))

        classified <- list(
            DESeq2 = if (!is.null(results$DESeq2))
                classify_genes(results$DESeq2, "gene", "log2FoldChange", "padj", config$alpha, config$log2fc_cutoff)
            else empty,

            edgeR = if (!is.null(results$edgeR))
                classify_genes(results$edgeR, "gene", "logFC", "FDR", config$alpha, config$log2fc_cutoff)
            else empty,

            limma = if (!is.null(results$limma))
                classify_genes(results$limma, "gene", "logFC", "adj.P.Val", config$alpha, config$log2fc_cutoff)
            else empty
        )

        for (direction in c("up", "down")) {

            message("Processing Venn diagram: ", group, " (", direction, ")")

            gene_sets <- list(
                DESeq2 = classified$DESeq2[[direction]],
                edgeR = classified$edgeR[[direction]],
                limma = classified$limma[[direction]]
            )

            title <- paste0(group, " (", direction, "regulated)")

            venn_grob <- build_method_venn(gene_sets, title)

            filename_base <- paste0(
                "Venn_",
                gsub("[^A-Za-z0-9_\\-]", "_", group),
                "_", direction
            )

            save_venn_plot(venn_grob, filename_base, figure_dir, config)

            regions <- venn_regions(gene_sets)

            overlap_counts[[length(overlap_counts) + 1]] <- data.frame(
                group = group,
                direction = direction,
                region = names(regions),
                n_genes = vapply(regions, length, integer(1)),
                stringsAsFactors = FALSE
            )

            all_genes <- Reduce(union, gene_sets)

            if (length(all_genes) > 0) {

                gene_membership[[length(gene_membership) + 1]] <- data.frame(
                    group = group,
                    direction = direction,
                    gene = all_genes,
                    DESeq2 = all_genes %in% gene_sets$DESeq2,
                    edgeR = all_genes %in% gene_sets$edgeR,
                    limma = all_genes %in% gene_sets$limma,
                    stringsAsFactors = FALSE
                )

            }

        }

    }

    overlap_counts <- do.call(rbind, overlap_counts)
    gene_membership <- do.call(rbind, gene_membership)

    save_csv(
        overlap_counts,
        file.path(figure_dir, "venn_overlap_counts.csv"),
        row.names = FALSE
    )

    save_csv(
        gene_membership,
        file.path(figure_dir, "venn_gene_membership.csv"),
        row.names = FALSE
    )

    message("Venn diagram comparison complete.")

    return(
        list(
            overlap_counts = overlap_counts,
            gene_membership = gene_membership
        )
    )

}
