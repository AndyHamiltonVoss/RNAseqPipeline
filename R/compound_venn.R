################################################################################
# RNAseqPipeline
#
# File: compound_venn.R
#
# Description:
#   Venn diagram comparing the consensus regulated gene set between two
#   compounds, each taken at its own peak-regulation dose. Reuses
#   classify_genes() / load_method_results() from dose_response.R and
#   save_venn_plot() from venn_diagrams.R.
#
################################################################################

#===============================================================================
# Get Consensus Regulated Genes for One Contrast
#===============================================================================

get_consensus_genes <- function(group, direction, config) {

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

    Reduce(
        intersect,
        list(classified$DESeq2[[direction]], classified$edgeR[[direction]], classified$limma[[direction]])
    )

}

#===============================================================================
# Find the Dose With the Most Consensus Regulated Genes for a Compound
#===============================================================================

find_peak_dose <- function(compound, groups, direction, config) {

    series <- get_dose_series(groups, compound)

    if (nrow(series) == 0) {
        stop("No dose series found for compound: ", compound)
    }

    series$n_genes <- vapply(
        series$group,
        function(group) length(get_consensus_genes(group, direction, config)),
        integer(1)
    )

    peak <- series[which.max(series$n_genes), ]

    return(peak)

}

#===============================================================================
# Build a Two-Way Venn Diagram Grob
#===============================================================================

build_two_way_venn <- function(gene_sets, title) {

    futile.logger::flog.threshold(
        futile.logger::ERROR,
        name = "VennDiagramLogger"
    )

    grob <- VennDiagram::venn.diagram(
        x = gene_sets,
        filename = NULL,
        fill = c("#1b9e77", "#d95f02"),
        alpha = 0.5,
        cex = 1.3,
        cat.cex = 1.2,
        main = title,
        main.cex = 1.3,
        margin = 0.1
    )

    return(grob)

}

#===============================================================================
# Compare Two Compounds at Their Peak-Regulation Dose
#===============================================================================

run_compound_comparison <- function(
    compound_a,
    compound_b,
    config,
    metadata,
    direction = "up"
) {

    message("Running compound comparison: ", compound_a, " vs ", compound_b)

    if (!requireNamespace("VennDiagram", quietly = TRUE)) {
        stop("Package 'VennDiagram' is required.")
    }

    groups <- setdiff(
        levels(factor(metadata$group)),
        config$reference_group
    )

    peak_a <- find_peak_dose(compound_a, groups, direction, config)
    peak_b <- find_peak_dose(compound_b, groups, direction, config)

    message(
        compound_a, " peak dose: ", peak_a$group,
        " (", peak_a$n_genes, " ", direction, "regulated consensus genes)"
    )
    message(
        compound_b, " peak dose: ", peak_b$group,
        " (", peak_b$n_genes, " ", direction, "regulated consensus genes)"
    )

    genes_a <- get_consensus_genes(peak_a$group, direction, config)
    genes_b <- get_consensus_genes(peak_b$group, direction, config)

    gene_sets <- stats::setNames(
        list(genes_a, genes_b),
        c(peak_a$group, peak_b$group)
    )

    title <- paste0(
        peak_a$group, " vs ", peak_b$group,
        " (", direction, "regulated, consensus)"
    )

    venn_grob <- build_two_way_venn(gene_sets, title)

    figure_dir <- file.path(config$output_directory, "Figures", "Venn")

    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    filename_base <- paste0(
        "Venn_", compound_a, "_vs_", compound_b, "_", direction, "_peak_dose"
    )

    save_venn_plot(venn_grob, filename_base, figure_dir, config)

    only_a <- setdiff(genes_a, genes_b)
    only_b <- setdiff(genes_b, genes_a)
    shared <- intersect(genes_a, genes_b)

    overlap_summary <- data.frame(
        region = c(
            paste(peak_a$group, "only"),
            paste(peak_b$group, "only"),
            "Shared"
        ),
        n_genes = c(length(only_a), length(only_b), length(shared)),
        stringsAsFactors = FALSE
    )

    all_genes <- union(genes_a, genes_b)

    gene_membership <- data.frame(
        gene = all_genes,
        stringsAsFactors = FALSE
    )

    gene_membership[[peak_a$group]] <- gene_membership$gene %in% genes_a
    gene_membership[[peak_b$group]] <- gene_membership$gene %in% genes_b

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

    message("Compound comparison complete.")

    return(
        list(
            peak_a = peak_a,
            peak_b = peak_b,
            gene_sets = gene_sets,
            overlap_summary = overlap_summary,
            gene_membership = gene_membership
        )
    )

}
