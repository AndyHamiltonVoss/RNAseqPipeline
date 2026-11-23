################################################################################
# RNAseqPipeline
#
# File: compare_exp_rosi_ora.R
#
# Description:
#   One-off, cross-dataset KEGG over-representation comparison, built
#   first for Rosiglitazone and then generalized to any compound with a
#   matching sheet in the external reference file: an external, published
#   reference (Rayl, Nemetchek, Voss, Hughes, "Agonists of the nuclear
#   receptor PPARg can produce biased signaling," Molecular Pharmacology
#   -- 24h DESeq2 results vs. DMSO, gene SYMBOLs, one sheet per compound
#   named "<compound>_vs_DMSO") against this project's own AMP0786 and
#   AMP0961 DESeq2 output, at every pure (non-combo) dose tested for that
#   compound in each dataset.
#
#   Gene selection is intentionally simple and identical across all
#   sources: every gene with padj < alpha, NO log2FoldChange magnitude
#   cutoff, NOT split by up/down direction -- unlike R/enrichment.R's
#   classify_genes()-based up/down/combined convention (which needs a
#   project-specific log2fc_cutoff that the external paper's own
#   methodology may not match), so the same rule can be applied uniformly
#   everywhere.
#
#   KEGG ORA mechanics (bitr -> enrichKEGG) mirror R/enrichment.R's
#   run_kegg_enrichment()/save_enrichment_result() exactly, generalized
#   with a `from_type` argument so the same code handles both this
#   project's own Ensembl gene IDs and the external file's gene symbols.
#
#   Standalone script, not wired into DGE_pipeline.R -- this is a
#   cross-config comparison (reads from both config.R's and
#   config_AMP0961.R's output directories plus an external file), the
#   same precedent as R/prepare_combined_experiment_data.R.
#
################################################################################

#===============================================================================
# Read the External Reference's <Compound>-vs-DMSO Sheet
#===============================================================================

read_external_compound_results <- function(path, sheet) {

    if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required.")
    }

    raw <- readxl::read_excel(path, sheet = sheet, skip = 1)

    names(raw)[1] <- "gene"

    as.data.frame(raw)

}

#===============================================================================
# KEGG ORA on Genes Passing padj < alpha Only (No LFC Cutoff, Not Split)
#===============================================================================
#
# `pvalue_cutoff` controls only how far down clusterProfiler's OWN output
# filtering goes, not which genes are selected as the foreground (that's
# always padj < alpha, per run_kegg_ora_padj_only()'s name). Pass 1 (via
# run_kegg_ora_dense()) to get every KEGG pathway clusterProfiler actually
# tested against the universe, not just the ones that cleared
# pvalueCutoff -- needed for a fair correlation across sources (see
# run_compare_exp_pathway_correlation()): a compound/dose with very few
# significant genes (e.g. AMP0961 Rosi_0.125nM) returns ZERO pathways at
# pvalueCutoff = alpha, and R/gsea_compound_correlation.R already found
# correlating only-significant-pathway subsets across sources with very
# different gene-list sizes gives unstable, misleading correlations.
#===============================================================================

run_kegg_ora_padj_only <- function(results_table, gene_col, padj_col, from_type, org_db, kegg_organism, alpha, pvalue_cutoff = alpha) {

    padj <- suppressWarnings(as.numeric(results_table[[padj_col]]))

    selected_genes <- unique(results_table[[gene_col]][!is.na(padj) & padj < alpha])

    universe_genes <- unique(results_table[[gene_col]][!is.na(results_table[[gene_col]])])

    gene_entrez <- suppressMessages(
        clusterProfiler::bitr(selected_genes, fromType = from_type, toType = "ENTREZID", OrgDb = org_db)
    )

    universe_entrez <- suppressMessages(
        clusterProfiler::bitr(universe_genes, fromType = from_type, toType = "ENTREZID", OrgDb = org_db)
    )

    result <- if (nrow(gene_entrez) == 0) {
        NULL
    } else {
        tryCatch(
            clusterProfiler::enrichKEGG(
                gene = unique(gene_entrez$ENTREZID),
                universe = unique(universe_entrez$ENTREZID),
                organism = kegg_organism,
                pAdjustMethod = "BH",
                pvalueCutoff = pvalue_cutoff,
                qvalueCutoff = pvalue_cutoff
            ),
            error = function(e) {
                warning("KEGG ORA failed: ", conditionMessage(e))
                NULL
            }
        )
    }

    list(
        result = result,
        n_genes_selected = length(selected_genes),
        n_genes_mapped_to_entrez = nrow(gene_entrez),
        n_universe_genes = length(universe_genes)
    )

}

run_kegg_ora_dense <- function(results_table, gene_col, padj_col, from_type, org_db, kegg_organism, alpha) {
    run_kegg_ora_padj_only(results_table, gene_col, padj_col, from_type, org_db, kegg_organism, alpha, pvalue_cutoff = 1)
}

#===============================================================================
# Map a Source's Tested Genes to ENTREZID, Keeping log2FoldChange and the
# padj < alpha "Selected" Subset (the Common ID Space Needed to Compare
# the External SYMBOL-Based File Against This Project's Own Ensembl IDs)
#===============================================================================

get_gene_lfc_by_entrez <- function(results_table, gene_col, lfc_col, padj_col, alpha, from_type, org_db) {

    tested <- data.frame(
        gene = results_table[[gene_col]],
        log2FoldChange = suppressWarnings(as.numeric(results_table[[lfc_col]])),
        padj = suppressWarnings(as.numeric(results_table[[padj_col]])),
        stringsAsFactors = FALSE
    )

    tested <- tested[!is.na(tested$gene) & !is.na(tested$log2FoldChange), ]

    mapped <- suppressMessages(
        clusterProfiler::bitr(unique(tested$gene), fromType = from_type, toType = "ENTREZID", OrgDb = org_db)
    )

    merged <- merge(tested, mapped, by.x = "gene", by.y = from_type)
    merged <- merged[!duplicated(merged$ENTREZID), ]

    sig <- !is.na(merged$padj) & merged$padj < alpha

    list(
        lfc_by_entrez = stats::setNames(merged$log2FoldChange, merged$ENTREZID),
        selected_entrez = unique(merged$ENTREZID[sig]),
        selected_up_entrez = unique(merged$ENTREZID[sig & merged$log2FoldChange > 0]),
        selected_down_entrez = unique(merged$ENTREZID[sig & merged$log2FoldChange < 0])
    )

}

#===============================================================================
# Save a Compare-Exp KEGG Result (Table + Dotplot) -- Mirrors
# save_enrichment_result() in R/enrichment.R
#===============================================================================

save_compare_exp_result <- function(result, filename_base, table_dir, figure_dir, figure_width, figure_height, figure_dpi, save_png, save_pdf) {

    if (is.null(result)) {
        return(invisible(FALSE))
    }

    result_df <- as.data.frame(result)

    save_csv(result_df, file.path(table_dir, paste0(filename_base, ".csv")), row.names = FALSE)

    if (nrow(result_df) == 0) {
        return(invisible(FALSE))
    }

    n_show <- min(20, nrow(result_df))

    dot_plot <- enrichplot::dotplot(result, showCategory = n_show) +
        ggplot2::ggtitle(gsub("_", " ", filename_base))

    plot_height <- max(figure_height, n_show * 0.5)

    if (isTRUE(save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_base, "_dotplot.png")),
            plot = dot_plot, width = figure_width, height = plot_height, dpi = figure_dpi
        )

    }

    if (isTRUE(save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_base, "_dotplot.pdf")),
            plot = dot_plot, width = figure_width, height = plot_height
        )

    }

    invisible(TRUE)

}

#===============================================================================
# Build the List of (Source, Label, Results Table, Column Names, ID Type)
# Tuples Shared by run_compare_exp_rosi_ora() and
# run_compare_exp_pathway_correlation()
#===============================================================================
#
# `compound` selects the external file's sheet ("<compound>_vs_DMSO") and
# the DESeq2 group-name prefix in both internal datasets
# ("DESeq2_group_<compound>_<dose>_vs_<reference_group>_raw.csv"). The
# dose lists are compound-specific (and dataset-specific within a
# compound) because AMP0786 and AMP0961 each test a different dose
# series for a given compound -- there's no way to discover them
# automatically here without also picking up combo-treatment contrasts
# (e.g. "Rosi_1000nM_T007_0.1uM"), which are deliberately excluded by
# only ever constructing pure "<compound>_<dose>" group names.
#===============================================================================

get_compare_exp_sources <- function(
    compound = "Rosi",
    amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"),
    amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    include_external = TRUE
) {

    sources <- list()

    if (include_external) {

        sources[[length(sources) + 1]] <- list(
            source_label = "external", label = paste0(compound, "_vs_DMSO_24h"),
            results_table = read_external_compound_results(external_path, paste0(compound, "_vs_DMSO")),
            gene_col = "gene", padj_col = "padj", from_type = "SYMBOL"
        )

    }

    for (dose in amp0786_doses) {

        f <- file.path(amp0786_dir, paste0("DESeq2_group_", compound, "_", dose, "_vs_DMSO_raw.csv"))

        if (!file.exists(f)) {
            warning("Missing AMP0786 file, skipping: ", f)
            next
        }

        sources[[length(sources) + 1]] <- list(
            source_label = "AMP0786", label = paste0(compound, "_", dose),
            results_table = utils::read.csv(f),
            gene_col = "gene", padj_col = "padj", from_type = "ENSEMBL"
        )

    }

    for (dose in amp0961_doses) {

        f <- file.path(amp0961_dir, paste0("DESeq2_group_", compound, "_", dose, "_vs_Control_raw.csv"))

        if (!file.exists(f)) {
            warning("Missing AMP0961 file, skipping: ", f)
            next
        }

        sources[[length(sources) + 1]] <- list(
            source_label = "AMP0961", label = paste0(compound, "_", dose),
            results_table = utils::read.csv(f),
            gene_col = "gene", padj_col = "padj", from_type = "ENSEMBL"
        )

    }

    sources

}

#===============================================================================
# Run the Full Cross-Dataset Rosi KEGG ORA Comparison
#===============================================================================

run_compare_exp_rosi_ora <- function(
    compound = "Rosi",
    amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"),
    amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    include_external = TRUE,
    output_dir = file.path("results/compare_exp", compound),
    alpha = 0.05,
    figure_width = 8,
    figure_height = 6,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-dataset %s KEGG ORA comparison...", compound))

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required.")
    }

    org_pkg <- get_organism_db("human")
    org_db <- get(org_pkg, envir = asNamespace(org_pkg))
    kegg_organism <- get_kegg_organism_code("human")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    sources <- get_compare_exp_sources(compound, amp0786_doses, amp0961_doses, external_path, amp0786_dir, amp0961_dir, include_external)

    summary_rows <- list()

    for (src in sources) {

        message(sprintf("Running KEGG ORA: %s (%s)...", src$source_label, src$label))

        ora <- run_kegg_ora_padj_only(
            src$results_table, src$gene_col, src$padj_col, src$from_type, org_db, kegg_organism, alpha
        )

        filename_base <- paste0("KEGG_", src$source_label, "_", src$label)

        save_compare_exp_result(
            ora$result, filename_base, table_dir, figure_dir,
            figure_width, figure_height, figure_dpi, save_png, save_pdf
        )

        n_sig_terms <- if (is.null(ora$result)) 0 else sum(as.data.frame(ora$result)$p.adjust < alpha, na.rm = TRUE)

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
            source = src$source_label,
            label = src$label,
            n_genes_selected = ora$n_genes_selected,
            n_genes_mapped_to_entrez = ora$n_genes_mapped_to_entrez,
            n_universe_genes = ora$n_universe_genes,
            n_significant_kegg_terms = n_sig_terms,
            stringsAsFactors = FALSE
        )

    }

    summary_table <- do.call(rbind, summary_rows)

    save_csv(
        summary_table,
        file.path(table_dir, sprintf("compare_exp_%s_ora_summary.csv", tolower(compound))),
        row.names = FALSE
    )

    message(sprintf("Cross-dataset %s KEGG ORA comparison complete.", compound))

    summary_table

}

#===============================================================================
# Run Source-vs-Source KEGG Pathway Correlation
#===============================================================================
#
# Pearson correlation of -log10(pvalue) across EVERY KEGG pathway each
# source's ORA actually tested (run_kegg_ora_dense(), pvalueCutoff = 1),
# not just the pathways significant in run_compare_exp_rosi_ora()'s
# saved tables. Mirrors R/gsea_compound_correlation.R's "Dense"
# convention and its documented reason: correlating only
# already-significant subsets is unstable when sources have very
# different numbers of significant genes (here: 162 for AMP0961
# Rosi_0.125nM vs. 4445 for AMP0961 Rosi_1000nM), and breaks outright
# for a source with zero significant pathways.
#===============================================================================

run_compare_exp_pathway_correlation <- function(
    compound = "Rosi",
    amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"),
    amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    include_external = TRUE,
    output_dir = file.path("results/compare_exp", compound),
    alpha = 0.05,
    figure_width = 8,
    figure_height = 6,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-dataset %s KEGG pathway correlation (dense)...", compound))

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required.")
    }

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    org_pkg <- get_organism_db("human")
    org_db <- get(org_pkg, envir = asNamespace(org_pkg))
    kegg_organism <- get_kegg_organism_code("human")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    sources <- get_compare_exp_sources(compound, amp0786_doses, amp0961_doses, external_path, amp0786_dir, amp0961_dir, include_external)

    source_names <- vapply(sources, function(s) paste0(s$source_label, "_", s$label), character(1))

    per_source <- lapply(sources, function(src) {

        message(sprintf("Running dense KEGG ORA: %s (%s)...", src$source_label, src$label))

        ora <- run_kegg_ora_dense(
            src$results_table, src$gene_col, src$padj_col, src$from_type, org_db, kegg_organism, alpha
        )

        if (is.null(ora$result)) {
            return(data.frame(ID = character(0), Description = character(0), neg_log10_p = numeric(0)))
        }

        result_df <- as.data.frame(ora$result)

        data.frame(
            ID = result_df$ID,
            Description = result_df$Description,
            neg_log10_p = -log10(result_df$pvalue),
            stringsAsFactors = FALSE
        )

    })

    names(per_source) <- source_names

    #===========================================================================
    # Build the Pathway x Source neg_log10_p Matrix
    #===========================================================================

    score_matrix <- Reduce(
        function(a, b) merge(a, b, by = "ID", all = TRUE),
        lapply(source_names, function(nm) {
            df <- per_source[[nm]][, c("ID", "neg_log10_p")]
            names(df)[2] <- nm
            df
        })
    )

    rownames(score_matrix) <- score_matrix$ID
    score_matrix$ID <- NULL
    score_matrix <- score_matrix[, source_names, drop = FALSE]

    correlation_matrix <- stats::cor(
        as.matrix(score_matrix),
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    save_csv(correlation_matrix, file.path(table_dir, sprintf("compare_exp_%s_pathway_correlation.csv", tolower(compound))))

    #===========================================================================
    # Heatmap
    #===========================================================================

    heatmap_height <- max(figure_height, length(source_names) * 0.5)
    heatmap_width <- heatmap_height + 2

    can_cluster <- !anyNA(correlation_matrix)

    if (!can_cluster) {
        message("Some source pairs share no tested KEGG pathways in common (NA correlation); disabling heatmap clustering.")
    }

    for (ext in c("png", "pdf")) {

        save_this_ext <- if (ext == "png") isTRUE(save_png) else isTRUE(save_pdf)

        if (!save_this_ext) {
            next
        }

        args <- list(
            correlation_matrix,
            cluster_rows = can_cluster,
            cluster_cols = can_cluster,
            display_numbers = TRUE,
            number_format = "%.2f",
            main = sprintf("Cross-Dataset %s KEGG Pathway Correlation (-log10 p, All Tested Pathways)", compound),
            filename = file.path(figure_dir, sprintf("Compare_exp_%s_pathway_correlation_heatmap.%s", tolower(compound), ext)),
            width = heatmap_width,
            height = heatmap_height
        )

        if (ext == "png") {
            args$res <- figure_dpi
        }

        do.call(pheatmap::pheatmap, args)

    }

    message(sprintf("Cross-dataset %s KEGG pathway correlation complete.", compound))

    correlation_matrix

}

#===============================================================================
# Run Source-vs-Source Correlation of the Genes Selected for the Pathway
# Analysis
#===============================================================================
#
# Pearson correlation of log2FoldChange across the UNION of genes that
# were "selected" (padj < alpha) for KEGG ORA in AT LEAST ONE of the 9
# sources -- not the full tested-gene universe (that's a broader
# question R/gene_rank_overlap.R already covers), and not restricted to
# genes significant in BOTH sources of a given pair (same instability
# this file's pathway correlation avoids -- see
# run_compare_exp_pathway_correlation()'s header comment). A source
# contributes its own log2FoldChange for a gene in this union even when
# that gene wasn't significant for that particular source, as long as it
# was tested, so the matrix stays dense rather than mostly NA.
#
# Genes are mapped to ENTREZID (via get_gene_lfc_by_entrez()) as the
# common ID space needed to compare the external file's gene SYMBOLs
# against this project's own Ensembl IDs.
#===============================================================================

build_and_save_gene_lfc_correlation <- function(
    per_source, source_names, gene_union, org_db,
    table_dir, figure_dir, compound, file_tag, title_suffix,
    figure_width, figure_height, figure_dpi, save_png, save_pdf
) {

    if (length(gene_union) == 0) {
        warning("No genes in this set (", file_tag, "); skipping correlation.")
        return(NULL)
    }

    lfc_matrix <- vapply(
        source_names,
        function(nm) unname(per_source[[nm]]$lfc_by_entrez[gene_union]),
        numeric(length(gene_union))
    )

    rownames(lfc_matrix) <- gene_union

    symbol_lookup <- suppressMessages(
        clusterProfiler::bitr(gene_union, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org_db)
    )

    gene_symbol <- stats::setNames(symbol_lookup$SYMBOL, symbol_lookup$ENTREZID)[gene_union]

    lfc_matrix_out <- as.data.frame(lfc_matrix)
    lfc_matrix_out$ENTREZID <- rownames(lfc_matrix)
    lfc_matrix_out$gene_symbol <- unname(gene_symbol)
    lfc_matrix_out <- lfc_matrix_out[, c("ENTREZID", "gene_symbol", source_names)]

    save_csv(
        lfc_matrix_out,
        file.path(table_dir, sprintf("compare_exp_%s_selected_gene_lfc_matrix_%s.csv", tolower(compound), file_tag)),
        row.names = FALSE
    )

    correlation_matrix <- stats::cor(
        lfc_matrix,
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    save_csv(
        correlation_matrix,
        file.path(table_dir, sprintf("compare_exp_%s_selected_gene_lfc_correlation_%s.csv", tolower(compound), file_tag))
    )

    heatmap_height <- max(figure_height, length(source_names) * 0.5)
    heatmap_width <- heatmap_height + 2

    can_cluster <- !anyNA(correlation_matrix)

    if (!can_cluster) {
        message(
            "Some source pairs share no commonly-tested genes in this set (",
            file_tag, "); disabling heatmap clustering."
        )
    }

    for (ext in c("png", "pdf")) {

        save_this_ext <- if (ext == "png") isTRUE(save_png) else isTRUE(save_pdf)

        if (!save_this_ext) {
            next
        }

        args <- list(
            correlation_matrix,
            cluster_rows = can_cluster,
            cluster_cols = can_cluster,
            display_numbers = TRUE,
            number_format = "%.2f",
            main = sprintf("Cross-Dataset %s Selected-Gene log2FoldChange Correlation (%s)", compound, title_suffix),
            filename = file.path(
                figure_dir,
                sprintf("Compare_exp_%s_selected_gene_lfc_correlation_heatmap_%s.%s", tolower(compound), file_tag, ext)
            ),
            width = heatmap_width,
            height = heatmap_height
        )

        if (ext == "png") {
            args$res <- figure_dpi
        }

        do.call(pheatmap::pheatmap, args)

    }

    correlation_matrix

}

run_compare_exp_gene_correlation <- function(
    compound = "Rosi",
    amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"),
    amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    include_external = TRUE,
    output_dir = file.path("results/compare_exp", compound),
    alpha = 0.05,
    figure_width = 8,
    figure_height = 6,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-dataset %s selected-gene log2FoldChange correlation...", compound))

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required.")
    }

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    org_pkg <- get_organism_db("human")
    org_db <- get(org_pkg, envir = asNamespace(org_pkg))

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    sources <- get_compare_exp_sources(compound, amp0786_doses, amp0961_doses, external_path, amp0786_dir, amp0961_dir, include_external)

    source_names <- vapply(sources, function(s) paste0(s$source_label, "_", s$label), character(1))

    per_source <- lapply(sources, function(src) {
        get_gene_lfc_by_entrez(src$results_table, src$gene_col, "log2FoldChange", src$padj_col, alpha, src$from_type, org_db)
    })

    names(per_source) <- source_names

    #===========================================================================
    # Combined (All Selected Genes, Either Direction), Up-Only, and
    # Down-Only Selected-Gene Unions
    #===========================================================================
    #
    # A gene's direction is source-specific (up in one source, down in
    # another) -- the up/down split is on a PER-SOURCE basis, same as
    # classify_genes() elsewhere in this project: "selected_up_entrez" for
    # a source means padj < alpha AND log2FoldChange > 0 IN THAT SOURCE.
    # The combined matrix still uses every source's own log2FoldChange
    # for a gene regardless of that source's own direction call, so a
    # gene that's genuinely discordant (up here, down there) still shows
    # up as a point dragging the correlation down, rather than being
    # silently dropped.
    #===========================================================================

    selected_union <- unique(unlist(lapply(per_source, `[[`, "selected_entrez")))
    selected_up_union <- unique(unlist(lapply(per_source, `[[`, "selected_up_entrez")))
    selected_down_union <- unique(unlist(lapply(per_source, `[[`, "selected_down_entrez")))

    message(sprintf(
        "%d genes selected (padj < %.2f) in at least one of %d sources (%d up, %d down).",
        length(selected_union), alpha, length(source_names), length(selected_up_union), length(selected_down_union)
    ))

    combined <- build_and_save_gene_lfc_correlation(
        per_source, source_names, selected_union, org_db, table_dir, figure_dir,
        compound, "combined", "All Selected Genes",
        figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    up <- build_and_save_gene_lfc_correlation(
        per_source, source_names, selected_up_union, org_db, table_dir, figure_dir,
        compound, "up", "Up-Regulated",
        figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    down <- build_and_save_gene_lfc_correlation(
        per_source, source_names, selected_down_union, org_db, table_dir, figure_dir,
        compound, "down", "Down-Regulated",
        figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    message(sprintf("Cross-dataset %s selected-gene log2FoldChange correlation complete.", compound))

    list(combined = combined, up = up, down = down)

}

#===============================================================================
# Scatter Plot: Up- vs. Down-Regulated Gene log2FoldChange Correlation,
# One Point Per Source Pair
#===============================================================================
#
# Reads the up/down correlation matrices run_compare_exp_gene_correlation()
# already saved (rather than recomputing) and plots every unique pair of
# the 9 sources as one point, up_r on x, down_r on y, against a dashed
# y = x reference line -- a direct visual of the up > down asymmetry
# already seen numerically (e.g. AMP0786_Rosi_1000nM vs. the external
# reference: up r = 0.46, down r = 0.26).
#===============================================================================

run_compare_exp_updown_scatter <- function(
    compound = "Rosi",
    output_dir = file.path("results/compare_exp", compound),
    figure_width = 9,
    figure_height = 7,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-dataset %s up- vs. down-correlation scatter plot...", compound))

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    up_matrix <- as.matrix(utils::read.csv(
        file.path(table_dir, sprintf("compare_exp_%s_selected_gene_lfc_correlation_up.csv", tolower(compound))),
        row.names = 1, check.names = FALSE
    ))

    down_matrix <- as.matrix(utils::read.csv(
        file.path(table_dir, sprintf("compare_exp_%s_selected_gene_lfc_correlation_down.csv", tolower(compound))),
        row.names = 1, check.names = FALSE
    ))

    source_names <- rownames(up_matrix)

    classify_pair <- function(a, b) {

        prefix <- function(x) sub("_.*$", "", x)

        pa <- prefix(a)
        pb <- prefix(b)

        if (pa == "external" || pb == "external") {
            "External vs. Internal"
        } else if (pa == pb) {
            paste0("Within ", pa)
        } else {
            "AMP0786 vs. AMP0961"
        }

    }

    pair_idx <- utils::combn(source_names, 2)

    pair_rows <- lapply(seq_len(ncol(pair_idx)), function(k) {

        a <- pair_idx[1, k]
        b <- pair_idx[2, k]

        data.frame(
            source_a = a, source_b = b,
            up_r = up_matrix[a, b], down_r = down_matrix[a, b],
            pair_type = classify_pair(a, b),
            stringsAsFactors = FALSE
        )

    })

    pair_table <- do.call(rbind, pair_rows)

    save_csv(
        pair_table,
        file.path(table_dir, sprintf("compare_exp_%s_updown_correlation_pairs.csv", tolower(compound))),
        row.names = FALSE
    )

    message(sprintf(
        "%d source pairs; up > down in %d, down > up in %d.",
        nrow(pair_table), sum(pair_table$up_r > pair_table$down_r), sum(pair_table$up_r < pair_table$down_r)
    ))

    axis_limits <- range(c(pair_table$up_r, pair_table$down_r), na.rm = TRUE)

    scatter_plot <- ggplot2::ggplot(pair_table, ggplot2::aes(x = up_r, y = down_r, color = pair_type)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::geom_point(size = 3, alpha = 0.85) +
        ggplot2::coord_equal(xlim = axis_limits, ylim = axis_limits) +
        ggplot2::labs(
            title = "Up- vs. Down-Regulated Gene log2FoldChange Correlation, by Source Pair",
            subtitle = paste0(
                sprintf("Each point = one pair of the %d %s sources (external paper + AMP0786 + AMP0961 doses)\n", length(source_names), compound),
                "Dashed line = equal up/down correlation; points below it mean up-regulated genes correlate more strongly"
            ),
            x = "Pearson correlation (log2FoldChange), up-regulated genes",
            y = "Pearson correlation (log2FoldChange), down-regulated genes",
            color = "Pair type"
        ) +
        ggplot2::theme_bw()

    if (isTRUE(save_png)) {
        ggplot2::ggsave(
            file.path(figure_dir, sprintf("Compare_exp_%s_updown_correlation_scatter.png", tolower(compound))),
            plot = scatter_plot, width = figure_width, height = figure_height, dpi = figure_dpi
        )
    }

    if (isTRUE(save_pdf)) {
        ggplot2::ggsave(
            file.path(figure_dir, sprintf("Compare_exp_%s_updown_correlation_scatter.pdf", tolower(compound))),
            plot = scatter_plot, width = figure_width, height = figure_height
        )
    }

    message(sprintf("Cross-dataset %s up- vs. down-correlation scatter plot complete.", compound))

    list(pair_table = pair_table, scatter_plot = scatter_plot)

}

#===============================================================================
# Count of Up- vs. Down-Regulated Genes Selected (padj < alpha) for ORA,
# Per Source
#===============================================================================
#
# How many of the padj < alpha genes fed into each source's KEGG ORA
# (run_compare_exp_rosi_ora()) were up- vs. down-regulated. Same
# per-source selected_up_entrez/selected_down_entrez sets
# get_gene_lfc_by_entrez() already computes for the up/down correlation
# split -- this just counts them instead of correlating their LFCs.
#===============================================================================

run_compare_exp_updown_gene_counts <- function(
    compound = "Rosi",
    amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"),
    amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    include_external = TRUE,
    output_dir = file.path("results/compare_exp", compound),
    alpha = 0.05,
    figure_width = 9,
    figure_height = 6,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-dataset %s up- vs. down-regulated selected gene counts...", compound))

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    org_pkg <- get_organism_db("human")
    org_db <- get(org_pkg, envir = asNamespace(org_pkg))

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    sources <- get_compare_exp_sources(compound, amp0786_doses, amp0961_doses, external_path, amp0786_dir, amp0961_dir, include_external)

    source_names <- vapply(sources, function(s) paste0(s$source_label, "_", s$label), character(1))

    per_source <- lapply(sources, function(src) {
        get_gene_lfc_by_entrez(src$results_table, src$gene_col, "log2FoldChange", src$padj_col, alpha, src$from_type, org_db)
    })

    names(per_source) <- source_names

    count_table <- data.frame(
        source = source_names,
        n_up = vapply(per_source, function(s) length(s$selected_up_entrez), integer(1)),
        n_down = vapply(per_source, function(s) length(s$selected_down_entrez), integer(1)),
        stringsAsFactors = FALSE
    )

    count_table$n_total <- count_table$n_up + count_table$n_down
    count_table$pct_up <- count_table$n_up / count_table$n_total * 100

    save_csv(
        count_table,
        file.path(table_dir, sprintf("compare_exp_%s_updown_gene_counts.csv", tolower(compound))),
        row.names = FALSE
    )

    #===========================================================================
    # Grouped Bar Chart
    #===========================================================================

    long_table <- rbind(
        data.frame(source = count_table$source, direction = "Up", n = count_table$n_up, stringsAsFactors = FALSE),
        data.frame(source = count_table$source, direction = "Down", n = count_table$n_down, stringsAsFactors = FALSE)
    )

    long_table$source <- factor(long_table$source, levels = source_names)
    long_table$direction <- factor(long_table$direction, levels = c("Down", "Up"))

    bar_plot <- ggplot2::ggplot(long_table, ggplot2::aes(x = source, y = n, fill = direction)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) +
        ggplot2::scale_fill_manual(values = c(Up = "darkorange", Down = "steelblue")) +
        ggplot2::labs(
            title = sprintf("Up- vs. Down-Regulated %s Genes Selected for KEGG ORA", compound),
            subtitle = sprintf("Genes with padj < %.2f, per source", alpha),
            x = "Source", y = "Number of genes selected", fill = "Direction"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 15)
        )

    if (isTRUE(save_png)) {
        ggplot2::ggsave(
            file.path(figure_dir, sprintf("Compare_exp_%s_updown_gene_counts.png", tolower(compound))),
            plot = bar_plot, width = figure_width, height = figure_height, dpi = figure_dpi
        )
    }

    if (isTRUE(save_pdf)) {
        ggplot2::ggsave(
            file.path(figure_dir, sprintf("Compare_exp_%s_updown_gene_counts.pdf", tolower(compound))),
            plot = bar_plot, width = figure_width, height = figure_height
        )
    }

    message(sprintf("Cross-dataset %s up- vs. down-regulated selected gene counts complete.", compound))

    list(count_table = count_table, bar_plot = bar_plot)

}

#===============================================================================
# Build Top-Dose Sources for Every AMP0961 Compound -- the "Different
# Drug" Null Baseline for the Concordance Proof
#===============================================================================
#
# Same top-dose discovery pattern as R/gene_rank_overlap.R /
# R/gsea_compound_correlation.R (parse_dose_group() + which.max(dose_nM)
# per compound, from R/dose_response.R), but surfaced here as its own
# standalone sources list rather than folded into one of those modules,
# since this one feeds the same (source_label, label, results_table,
# gene_col, padj_col, from_type) shape get_compare_exp_sources() returns
# -- every lower-level compare_exp helper (run_kegg_ora_dense(),
# get_gene_lfc_by_entrez(), build_and_save_gene_lfc_correlation()) is
# source-list-driven, not tied to get_compare_exp_sources() specifically,
# so they all apply here unchanged.
#===============================================================================

get_amp0961_all_compound_top_dose_sources <- function(
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    reference_group = "Control"
) {

    raw_files <- list.files(amp0961_dir, pattern = "^DESeq2_group_.*_raw\\.csv$")

    group_names <- sub(
        paste0("_vs_", reference_group, "_raw\\.csv$"), "",
        sub("^DESeq2_group_", "", raw_files)
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
        lapply(split(dose_table, dose_table$compound), function(d) d[which.max(d$dose_nM), ])
    )

    compounds <- sort(top_dose_table$compound)

    lapply(compounds, function(cmpd) {

        group <- top_dose_table$group[top_dose_table$compound == cmpd]
        f <- file.path(amp0961_dir, paste0("DESeq2_group_", group, "_vs_", reference_group, "_raw.csv"))

        list(
            source_label = "AMP0961", label = cmpd,
            results_table = utils::read.csv(f),
            gene_col = "gene", padj_col = "padj", from_type = "ENSEMBL"
        )

    })

}

#===============================================================================
# Run the "Different Drug" Null Baseline: All AMP0961 Compounds (Top
# Dose), Pairwise KEGG Pathway and Selected-Gene log2FoldChange
# Correlation
#===============================================================================
#
# Part 2 of the pasted control-anchoring methodology: show that drugs
# REPEATED across experiments (Rosi, MRL24, GW1929, T007 -- already
# correlated against each other via compare_exp/<compound>/) correlate
# more strongly than DIFFERENT drugs do within the same batch-corrected
# screen. This builds that "different drug" null distribution across
# every pair of AMP0961's ~20 compounds (top dose each), using the exact
# same KEGG-ORA/gene-LFC methodology as the matched-drug comparisons, so
# the two are directly comparable (see run_compare_exp_concordance_proof()).
#===============================================================================

run_amp0961_compound_null_baseline <- function(
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    reference_group = "Control",
    output_dir = "results/compare_exp/concordance_proof",
    alpha = 0.05,
    figure_width = 10,
    figure_height = 9,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message("Running AMP0961 all-compound top-dose null baseline...")

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required.")
    }

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    org_pkg <- get_organism_db("human")
    org_db <- get(org_pkg, envir = asNamespace(org_pkg))
    kegg_organism <- get_kegg_organism_code("human")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    sources <- get_amp0961_all_compound_top_dose_sources(amp0961_dir, reference_group)
    source_names <- vapply(sources, function(s) s$label, character(1))

    message(sprintf("Null baseline across %d AMP0961 compounds (top dose each).", length(sources)))

    #===========================================================================
    # Dense KEGG Pathway Correlation
    #===========================================================================

    per_source_pathway <- lapply(sources, function(src) {

        ora <- run_kegg_ora_dense(src$results_table, src$gene_col, src$padj_col, src$from_type, org_db, kegg_organism, alpha)

        if (is.null(ora$result)) {
            return(data.frame(ID = character(0), neg_log10_p = numeric(0)))
        }

        result_df <- as.data.frame(ora$result)

        data.frame(ID = result_df$ID, neg_log10_p = -log10(result_df$pvalue), stringsAsFactors = FALSE)

    })

    names(per_source_pathway) <- source_names

    pathway_score_matrix <- Reduce(
        function(a, b) merge(a, b, by = "ID", all = TRUE),
        lapply(source_names, function(nm) {
            df <- per_source_pathway[[nm]]
            names(df)[2] <- nm
            df
        })
    )

    rownames(pathway_score_matrix) <- pathway_score_matrix$ID
    pathway_score_matrix$ID <- NULL
    pathway_score_matrix <- pathway_score_matrix[, source_names, drop = FALSE]

    pathway_correlation_matrix <- stats::cor(
        as.matrix(pathway_score_matrix),
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    save_csv(pathway_correlation_matrix, file.path(table_dir, "amp0961_null_baseline_pathway_correlation.csv"))

    heatmap_height <- max(figure_height, length(source_names) * 0.3)
    heatmap_width <- heatmap_height + 2
    can_cluster <- !anyNA(pathway_correlation_matrix)

    if (!can_cluster) {
        message("Some compound pairs share no tested KEGG pathways in common (NA correlation); disabling heatmap clustering.")
    }

    for (ext in c("png", "pdf")) {

        save_this_ext <- if (ext == "png") isTRUE(save_png) else isTRUE(save_pdf)

        if (!save_this_ext) {
            next
        }

        args <- list(
            pathway_correlation_matrix,
            cluster_rows = can_cluster, cluster_cols = can_cluster,
            display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6,
            main = "AMP0961 Null Baseline: All-Compound KEGG Pathway Correlation (Top Dose)",
            filename = file.path(figure_dir, paste0("Amp0961_null_baseline_pathway_correlation_heatmap.", ext)),
            width = heatmap_width, height = heatmap_height
        )

        if (ext == "png") {
            args$res <- figure_dpi
        }

        do.call(pheatmap::pheatmap, args)

    }

    #===========================================================================
    # Selected-Gene log2FoldChange Correlation (Reuses
    # build_and_save_gene_lfc_correlation() Unchanged)
    #===========================================================================

    per_source_genes <- lapply(sources, function(src) {
        get_gene_lfc_by_entrez(src$results_table, src$gene_col, "log2FoldChange", src$padj_col, alpha, src$from_type, org_db)
    })

    names(per_source_genes) <- source_names

    selected_union <- unique(unlist(lapply(per_source_genes, `[[`, "selected_entrez")))

    gene_correlation_matrix <- build_and_save_gene_lfc_correlation(
        per_source_genes, source_names, selected_union, org_db, table_dir, figure_dir,
        "amp0961_null_baseline", "combined", "All AMP0961 Compounds, Top Dose",
        figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    message("AMP0961 all-compound top-dose null baseline complete.")

    list(
        pathway_correlation_matrix = pathway_correlation_matrix,
        gene_correlation_matrix = gene_correlation_matrix,
        sources = sources
    )

}

#===============================================================================
# Concordance Proof: Matched Same-Drug Correlations vs. the AMP0961
# Different-Drug Null Baseline
#===============================================================================
#
# Pulls the already-saved matched same-drug cross-experiment correlation
# (external-vs-AMP0961 and AMP0786-vs-AMP0961, both top dose, "combined"
# i.e. not up/down split) for every compound run through compare_exp/ so
# far (Rosi, MRL24, GW1929, T007), and reports each one's percentile rank
# within run_amp0961_compound_null_baseline()'s different-drug
# distribution -- the actual quantitative version of "repeated drugs
# concordant across experiments, unlike different drugs within the
# screen" the pasted methodology calls for.
#===============================================================================

run_compare_exp_concordance_proof <- function(
    output_dir = "results/compare_exp",
    null_baseline_dir = "results/compare_exp/concordance_proof",
    figure_width = 9,
    figure_height = 7,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message("Running concordance proof: matched same-drug vs. different-drug null baseline...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    table_dir <- file.path(null_baseline_dir, "Tables")
    figure_dir <- file.path(null_baseline_dir, "Figures")

    pathway_null <- as.matrix(utils::read.csv(
        file.path(table_dir, "amp0961_null_baseline_pathway_correlation.csv"),
        row.names = 1, check.names = FALSE
    ))

    gene_null <- as.matrix(utils::read.csv(
        file.path(table_dir, "compare_exp_amp0961_null_baseline_selected_gene_lfc_correlation_combined.csv"),
        row.names = 1, check.names = FALSE
    ))

    get_offdiag <- function(m) {
        diag(m) <- NA
        as.numeric(m[!is.na(m)])
    }

    pathway_null_values <- get_offdiag(pathway_null)
    gene_null_values <- get_offdiag(gene_null)

    percentile_rank <- function(x, distribution) mean(distribution <= x, na.rm = TRUE) * 100

    compound_output_dir <- function(compound) {
        file.path(output_dir, compound)
    }

    read_corr_value <- function(compound, filename, a, b) {

        f <- file.path(compound_output_dir(compound), "Tables", filename)
        m <- as.matrix(utils::read.csv(f, row.names = 1, check.names = FALSE))
        m[a, b]

    }

    #===========================================================================
    # Matched Same-Drug Pairs (Top Dose Each, From the compare_exp/ Runs
    # Already Completed This Session)
    #===========================================================================

    matched_pairs <- list(
        list(compound = "Rosi", partner_type = "External", amp0961_label = "AMP0961_Rosi_1000nM", partner_label = "external_Rosi_vs_DMSO_24h"),
        list(compound = "Rosi", partner_type = "AMP0786", amp0961_label = "AMP0961_Rosi_1000nM", partner_label = "AMP0786_Rosi_1000nM"),
        list(compound = "MRL24", partner_type = "External", amp0961_label = "AMP0961_MRL24_100nM", partner_label = "external_MRL24_vs_DMSO_24h"),
        list(compound = "MRL24", partner_type = "AMP0786", amp0961_label = "AMP0961_MRL24_100nM", partner_label = "AMP0786_MRL24_1000nM"),
        list(compound = "GW1929", partner_type = "External", amp0961_label = "AMP0961_GW1929_100nM", partner_label = "external_GW1929_vs_DMSO_24h"),
        list(compound = "T007", partner_type = "AMP0786", amp0961_label = "AMP0961_T007_100nM", partner_label = "AMP0786_T007_10uM")
    )

    rows <- lapply(matched_pairs, function(p) {

        pathway_file <- sprintf("compare_exp_%s_pathway_correlation.csv", tolower(p$compound))
        gene_file <- sprintf("compare_exp_%s_selected_gene_lfc_correlation_combined.csv", tolower(p$compound))

        pathway_r <- read_corr_value(p$compound, pathway_file, p$amp0961_label, p$partner_label)
        gene_r <- read_corr_value(p$compound, gene_file, p$amp0961_label, p$partner_label)

        data.frame(
            compound = p$compound, partner_type = p$partner_type,
            amp0961_label = p$amp0961_label, partner_label = p$partner_label,
            pathway_r = pathway_r, pathway_percentile = percentile_rank(pathway_r, pathway_null_values),
            gene_lfc_r = gene_r, gene_lfc_percentile = percentile_rank(gene_r, gene_null_values),
            stringsAsFactors = FALSE
        )

    })

    proof_table <- do.call(rbind, rows)

    save_csv(proof_table, file.path(table_dir, "concordance_proof_pairs.csv"), row.names = FALSE)

    message(sprintf(
        "Matched same-drug pairs: pathway percentile range %.0f-%.0f, gene-LFC percentile range %.0f-%.0f (out of 100 vs. the different-drug null; higher = more concordant than chance).",
        min(proof_table$pathway_percentile), max(proof_table$pathway_percentile),
        min(proof_table$gene_lfc_percentile), max(proof_table$gene_lfc_percentile)
    ))

    #===========================================================================
    # Strip Plot: Null Distribution (Grey) vs. Matched-Drug Points (Colored)
    #===========================================================================

    null_long <- rbind(
        data.frame(metric = "KEGG pathway correlation (-log10 p, all tested pathways)", value = pathway_null_values, stringsAsFactors = FALSE),
        data.frame(metric = "Selected-gene log2FoldChange correlation", value = gene_null_values, stringsAsFactors = FALSE)
    )

    proof_table$pair_label <- paste0(proof_table$compound, " (", proof_table$partner_type, ")")

    matched_long <- rbind(
        data.frame(metric = "KEGG pathway correlation (-log10 p, all tested pathways)", value = proof_table$pathway_r, pair_label = proof_table$pair_label, stringsAsFactors = FALSE),
        data.frame(metric = "Selected-gene log2FoldChange correlation", value = proof_table$gene_lfc_r, pair_label = proof_table$pair_label, stringsAsFactors = FALSE)
    )

    set.seed(1)

    strip_plot <- ggplot2::ggplot() +
        ggplot2::geom_jitter(
            data = null_long, ggplot2::aes(x = value, y = 0),
            height = 0.12, width = 0, alpha = 0.4, color = "grey50", size = 1.6
        ) +
        ggplot2::geom_point(
            data = matched_long, ggplot2::aes(x = value, y = 0.35, color = pair_label),
            size = 3.2
        )

    if (requireNamespace("ggrepel", quietly = TRUE)) {
        strip_plot <- strip_plot + ggrepel::geom_text_repel(
            data = matched_long, ggplot2::aes(x = value, y = 0.35, label = pair_label, color = pair_label),
            size = 3, show.legend = FALSE, direction = "both", min.segment.length = 0, seed = 1
        )
    } else {
        strip_plot <- strip_plot + ggplot2::geom_text(
            data = matched_long, ggplot2::aes(x = value, y = 0.5, label = pair_label, color = pair_label),
            size = 3, show.legend = FALSE
        )
    }

    strip_plot <- strip_plot +
        ggplot2::facet_wrap(~metric, ncol = 1, scales = "free_x") +
        ggplot2::scale_y_continuous(limits = c(-0.3, 0.8), breaks = NULL) +
        ggplot2::labs(
            title = "Concordance Proof: Same Drug Across Experiments vs. Different Drugs Within AMP0961",
            subtitle = paste0(
                "Grey = different-drug null baseline (all ~20 AMP0961 compounds, top dose, every pair)\n",
                "Colored = matched same-drug pairs (external paper or AMP0786 vs. AMP0961, top dose)"
            ),
            x = "Pearson correlation", y = NULL, color = "Matched pair"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(size = 9))

    if (isTRUE(save_png)) {
        ggplot2::ggsave(
            file.path(figure_dir, "Concordance_proof_scatter.png"),
            plot = strip_plot, width = figure_width, height = figure_height, dpi = figure_dpi
        )
    }

    if (isTRUE(save_pdf)) {
        ggplot2::ggsave(
            file.path(figure_dir, "Concordance_proof_scatter.pdf"),
            plot = strip_plot, width = figure_width, height = figure_height
        )
    }

    message("Concordance proof complete.")

    proof_table

}
