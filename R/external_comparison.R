################################################################################
# RNAseqPipeline
#
# File: external_comparison.R
#
# Description:
#   Compares this pipeline's DESeq2 results against an externally published
#   dataset (Rayl et al. 2024, Molecular Pharmacology). Three checks:
#     1. Reprocess the published normalized count table through this
#        pipeline's own DESeq2 module, to confirm it reproduces the paper's
#        reported results (reproducibility check).
#     2. Compare regulated gene sets from this study's dose-response
#        contrasts against the published significant gene calls (Venn).
#
################################################################################

#===============================================================================
# Read a Published DESeq2 Results Tab
#===============================================================================

read_published_deseq2 <- function(xlsx_path, sheet) {

    df <- readxl::read_excel(xlsx_path, sheet = sheet, skip = 1)
    df <- as.data.frame(df)
    colnames(df)[1] <- "gene"

    # The published padj is NA for many genes DESeq2 excluded via independent
    # filtering (typically low baseMean / high dispersion). Recompute
    # BH-adjusted padj directly from the raw p-values so those genes aren't
    # dropped outright.
    df$padj <- stats::p.adjust(df$pvalue, method = "BH")

    return(df)

}

#===============================================================================
# Read the Published Normalized Count Table
#===============================================================================

read_published_counts <- function(xlsx_path, sheet = "Normalized count table") {

    df <- readxl::read_excel(xlsx_path, sheet = sheet, skip = 2)
    df <- as.data.frame(df)
    colnames(df)[1] <- "gene"

    rownames(df) <- df$gene
    df$gene <- NULL

    return(df)

}

#===============================================================================
# Map a Set of Ensembl Gene IDs to Unique Gene Symbols
#===============================================================================

ensembl_to_symbol <- function(ensembl_ids, org_db) {

    if (length(ensembl_ids) == 0) {
        return(character(0))
    }

    mapped <- suppressMessages(
        clusterProfiler::bitr(
            ensembl_ids,
            fromType = "ENSEMBL",
            toType = "SYMBOL",
            OrgDb = org_db
        )
    )

    unique(mapped$SYMBOL)

}

#===============================================================================
# Build a Two-Way Venn Grob With Wider Category-Label Spacing
#
# This study's gene sets and the published gene sets often differ by 1-2
# orders of magnitude, which crowds the smaller circle's count callouts
# against its label under build_two_way_venn()'s default spacing.
#===============================================================================

build_external_venn <- function(gene_sets, title, margin = 0.15) {

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
        cat.dist = c(0.025, 0.09),
        cat.pos = c(-30, 30),
        main = title,
        main.cex = 1.3,
        margin = margin
    )

    return(grob)

}

#===============================================================================
# Reprocess the Published Normalized Counts Through This Pipeline's DESeq2
#===============================================================================

run_published_reprocessing <- function(xlsx_path, config) {

    message("Reprocessing published normalized counts through DESeq2...")

    counts <- read_published_counts(xlsx_path)

    external_config <- config
    external_config$output_directory <- file.path(config$output_directory, "External_DESeq2")
    external_config$reference_group <- "dmso"

    create_output_directories(external_config$output_directory)

    counts <- filter_low_counts(
        counts = counts,
        minimum_counts = external_config$minimum_counts,
        minimum_samples = external_config$minimum_samples
    )

    metadata <- create_metadata(
        counts = counts,
        replicate_pattern = external_config$replicate_pattern
    )

    validate_metadata(
        metadata = metadata,
        minimum_replicates = external_config$minimum_samples
    )

    save_metadata(
        metadata = metadata,
        output_directory = external_config$output_directory
    )

    run_deseq2(
        counts = counts,
        metadata = metadata,
        config = external_config
    )

}

#===============================================================================
# Compare Reprocessed Results Against the Published Results (Reproducibility)
#===============================================================================

run_reproducibility_check <- function(xlsx_path, config, results_dir) {

    message("Running reproducibility check against published results...")

    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

    figure_dir <- file.path(config$output_directory, "Figures", "Venn_Reproducibility")
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    reprocessed_dir <- file.path(config$output_directory, "External_DESeq2", "DESeq2", "Tables")

    pairs <- list(
        Rosi = list(
            sheet = "Rosi_vs_DMSO",
            reprocessed = file.path(reprocessed_dir, "DESeq2_group_rosi_vs_dmso_raw.csv")
        ),
        MRL24 = list(
            sheet = "MRL24_vs_DMSO",
            reprocessed = file.path(reprocessed_dir, "DESeq2_group_mrl24_vs_dmso_raw.csv")
        )
    )

    correlation_summary <- list()
    overlap_summary <- list()

    for (compound in names(pairs)) {

        published <- read_published_deseq2(xlsx_path, pairs[[compound]]$sheet)
        reprocessed <- read.csv(pairs[[compound]]$reprocessed)

        merged <- merge(
            published, reprocessed,
            by = "gene", suffixes = c("_published", "_reprocessed")
        )

        r <- stats::cor(
            merged$log2FoldChange_published,
            merged$log2FoldChange_reprocessed,
            use = "complete.obs"
        )

        message(compound, ": n=", nrow(merged), " Pearson r=", round(r, 4))

        correlation_summary[[length(correlation_summary) + 1]] <- data.frame(
            compound = compound,
            n_genes = nrow(merged),
            pearson_r = r
        )

        p <- ggplot2::ggplot(
            merged,
            ggplot2::aes(x = log2FoldChange_published, y = log2FoldChange_reprocessed)
        ) +
            ggplot2::geom_point(alpha = 0.3, size = 0.8, color = "#1b9e77") +
            ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
            ggplot2::labs(
                title = paste0(compound, " vs DMSO: reprocessed vs published log2FoldChange"),
                subtitle = paste0("Pearson r = ", round(r, 4), ", n = ", nrow(merged), " genes"),
                x = "Published log2FoldChange",
                y = "Reprocessed (this pipeline) log2FoldChange"
            ) +
            ggplot2::theme_bw()

        if (isTRUE(config$save_png)) {
            ggplot2::ggsave(
                file.path(figure_dir, paste0("Scatter_", compound, "_reprocessed_vs_published.png")),
                plot = p, width = config$figure_width, height = config$figure_height, dpi = config$figure_dpi
            )
        }

        if (isTRUE(config$save_pdf)) {
            ggplot2::ggsave(
                file.path(figure_dir, paste0("Scatter_", compound, "_reprocessed_vs_published.pdf")),
                plot = p, width = config$figure_width, height = config$figure_height
            )
        }

        # Significance-based classification (padj < alpha, any effect size):
        # the published paper's own reported gene counts use this convention.
        published_classified <- classify_genes(
            published, "gene", "log2FoldChange", "padj", config$alpha, 0
        )

        reprocessed_classified <- classify_genes(
            reprocessed, "gene", "log2FoldChange", "padj", config$alpha, 0
        )

        for (direction in c("up", "down")) {

            gene_sets <- stats::setNames(
                list(reprocessed_classified[[direction]], published_classified[[direction]]),
                c("Reprocessed", "Published")
            )

            title <- paste0(compound, " reprocessed vs published (", direction, "regulated)")

            venn_grob <- build_external_venn(gene_sets, title)

            filename_base <- paste0("Venn_", compound, "_reprocessed_vs_published_", direction)

            save_venn_plot(
                venn_grob, filename_base, figure_dir, config,
                width = config$figure_width * 1.3, height = config$figure_height
            )

            shared <- intersect(reprocessed_classified[[direction]], published_classified[[direction]])

            overlap_summary[[length(overlap_summary) + 1]] <- data.frame(
                compound = compound,
                direction = direction,
                n_reprocessed = length(reprocessed_classified[[direction]]),
                n_published = length(published_classified[[direction]]),
                n_shared = length(shared),
                stringsAsFactors = FALSE
            )

        }

    }

    correlation_summary <- do.call(rbind, correlation_summary)
    overlap_summary <- do.call(rbind, overlap_summary)

    save_csv(correlation_summary, file.path(results_dir, "reproducibility_correlation_summary.csv"), row.names = FALSE)
    save_csv(overlap_summary, file.path(results_dir, "reproducibility_overlap_summary.csv"), row.names = FALSE)

    message("Reproducibility check complete.")

    return(
        list(
            correlation_summary = correlation_summary,
            overlap_summary = overlap_summary
        )
    )

}

#===============================================================================
# Compare This Study's Dose-Response Gene Sets Against the Published Gene Sets
#===============================================================================

run_dose_comparison <- function(xlsx_path, config, results_dir, doses = c("1nM", "10nM", "100nM", "1000nM")) {

    message("Comparing this study's dose-response contrasts to published gene sets...")

    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

    figure_dir <- file.path(config$output_directory, "Figures", "Venn_External")
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    org_pkg <- get_organism_db(config$organism)
    org_db <- get(org_pkg, envir = asNamespace(org_pkg))

    published <- list(
        Rosi = read_published_deseq2(xlsx_path, "Rosi_vs_DMSO"),
        MRL24 = read_published_deseq2(xlsx_path, "MRL24_vs_DMSO")
    )

    # No log2FC magnitude cutoff for this comparison: classify purely on
    # significance (padj < alpha) and direction of change, since the
    # published study's reported gene counts use this same definition.
    published_classified <- lapply(
        published,
        function(df) classify_genes(df, "gene", "log2FoldChange", "padj", config$alpha, 0)
    )

    overlap_summaries <- list()

    for (compound in names(published)) {

        for (dose in doses) {

            group <- paste0(compound, "_", dose)

            deseq2_file <- file.path(
                config$output_directory, "DESeq2", "Tables",
                paste0("DESeq2_group_", group, "_vs_", config$reference_group, "_shrunk.csv")
            )

            if (!file.exists(deseq2_file)) {
                message("Skipping ", group, ": no DESeq2 results found.")
                next
            }

            mine <- read.csv(deseq2_file)

            mine_classified <- classify_genes(
                mine, "gene", "log2FoldChange", "padj", config$alpha, 0
            )

            for (direction in c("up", "down")) {

                message("Processing external comparison: ", group, " (", direction, ")")

                mine_symbols <- ensembl_to_symbol(mine_classified[[direction]], org_db)
                published_symbols <- published_classified[[compound]][[direction]]

                gene_sets <- stats::setNames(
                    list(mine_symbols, published_symbols),
                    c(group, "Published")
                )

                title <- paste0(group, " vs Published ", compound, " (", direction, "regulated)")

                venn_grob <- build_external_venn(gene_sets, title)

                filename_base <- paste0("Venn_", group, "_vs_Published_", direction)

                save_venn_plot(
                    venn_grob, filename_base, figure_dir, config,
                    width = config$figure_width * 1.3, height = config$figure_height
                )

                shared <- intersect(mine_symbols, published_symbols)

                overlap_summaries[[length(overlap_summaries) + 1]] <- data.frame(
                    compound = compound,
                    dose = dose,
                    direction = direction,
                    n_this_study = length(mine_symbols),
                    n_published = length(published_symbols),
                    n_shared = length(shared),
                    stringsAsFactors = FALSE
                )

            }

        }

    }

    overlap_summary <- do.call(rbind, overlap_summaries)

    save_csv(overlap_summary, file.path(results_dir, "dose_comparison_overlap_summary.csv"), row.names = FALSE)

    message("Dose comparison complete.")

    return(overlap_summary)

}
