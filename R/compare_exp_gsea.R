################################################################################
# RNAseqPipeline
#
# File: compare_exp_gsea.R
#
# Description:
#   GSEA (rank-based, shrunken log2FoldChange) parallel to
#   R/compare_exp_rosi_ora.R's over-representation (ORA) suite, for the
#   same cross-dataset comparison: an external published reference vs.
#   AMP0786 vs. AMP0961, per compound, at every pure dose tested.
#
#   "GSEA ranked by shrunken log2FoldChange" is this project's own
#   pre-existing convention -- R/gsea.R's run_gsea(rank_by =
#   "log2FoldChange"), which writes to a "GSEA_LFC" directory. This file
#   reuses R/gsea.R's low-level building blocks directly
#   (build_ranked_gene_vector, build_ranked_entrez_vector, run_gsea_kegg,
#   split_gsea_result_by_significance, save_gsea_result) rather than
#   reimplementing gseKEGG plumbing, via a small local "config stub" list
#   in place of a real config object (compare_exp's whole module stays
#   config-independent, explicit function arguments instead).
#
#   Internal sources (AMP0786/AMP0961) read "*_shrunk.csv" -- the same
#   apeglm-shrunk files GSEA_LFC already reads elsewhere in this project,
#   including AMP0961's documented stale-but-accepted leftovers (see
#   R/gsea.R's header; not a new problem introduced here).
#
#   The external published reference has no raw counts available (only a
#   DESeq2-normalized count table, confirmed non-integer), so apeglm
#   shrinkage can't be refit there. Instead, apply_ashr_shrinkage() calls
#   ashr::ash() directly on the external sheet's own already-present
#   log2FoldChange + lfcSE columns -- the same thing
#   DESeq2::lfcShrink(type = "ashr") does internally, but usable without
#   a DESeqDataSet since ashr's method only needs a point estimate and
#   its standard error per gene.
#
#   GSEA doesn't threshold genes the way ORA does (ORA: a fixed
#   significant-gene list at padj < alpha; GSEA: the full ranked vector).
#   For the companion figures ORA built from its selected-gene list
#   (gene-level correlation, up/down split, up/down counts), this keeps
#   the SAME padj < alpha selection (same column, now read from the
#   shrunk/ashr source) but correlates the SHRUNKEN log2FoldChange
#   instead of raw -- same gene membership as the ORA figures, different
#   (shrunken) effect-size input, so the two stay directly comparable.
#   Reuses get_gene_lfc_by_entrez()/build_and_save_gene_lfc_correlation()
#   from R/compare_exp_rosi_ora.R unchanged (both are already generic
#   over lfc_col/padj_col, not ORA-specific).
#
#   Pathway correlation uses NES (not -log10(p) like the ORA version) --
#   deliberately follows this project's OWN established GSEA convention,
#   R/gsea_compound_correlation.R, which already correlates NES as GSEA's
#   primary comparable score (direction + magnitude in one number).
#
#   The KEGG dot plot's dot SIZE is leading-edge gene count (split
#   core_enrichment on "/", count -- the same idiom
#   R/gsea_compound_correlation.R already uses for its own leading-edge
#   work) in place of ORA's Count column (no direct analog in a gseaResult).
#
#   Output nested under results/compare_exp/<compound>/GSEA_LFC_<alpha>/,
#   a sibling of the ORA_<alpha>/ directory already there.
#
################################################################################

#===============================================================================
# ashr Shrinkage for the External Reference (No Raw Counts Available)
#===============================================================================

apply_ashr_shrinkage <- function(results_table, lfc_col = "log2FoldChange", lfcse_col = "lfcSE") {

    if (!requireNamespace("ashr", quietly = TRUE)) {
        stop("Package 'ashr' is required.")
    }

    lfc <- suppressWarnings(as.numeric(results_table[[lfc_col]]))
    lfcse <- suppressWarnings(as.numeric(results_table[[lfcse_col]]))

    valid <- !is.na(lfc) & !is.na(lfcse)

    fit <- ashr::ash(
        betahat = lfc[valid],
        sebetahat = lfcse[valid],
        mixcompdist = "normal",
        method = "shrink"
    )

    results_table[[lfc_col]][valid] <- fit$result$PosteriorMean

    results_table

}

#===============================================================================
# A Small Config Stub for Reusing R/gsea.R's config-Driven Functions
#===============================================================================

make_gsea_config_stub <- function(alpha, figure_width = 8, figure_height = 6, figure_dpi = 600, save_png = TRUE, save_pdf = TRUE) {

    list(
        gsea_min_gene_set_size = 10,
        gsea_max_gene_set_size = 500,
        random_seed = 1234,
        alpha = alpha,
        run_gsea_dense = TRUE,
        figure_width = figure_width,
        figure_height = figure_height,
        figure_dpi = figure_dpi,
        save_png = save_png,
        save_pdf = save_pdf
    )

}

#===============================================================================
# Map a Ranked Gene Vector to ENTREZID for gseKEGG -- Generalized Over
# `from_type`
#===============================================================================
#
# R/gsea.R's own build_ranked_entrez_vector() hardcodes fromType =
# "ENSEMBL" (correct for every OTHER caller in this project, which only
# ever ranks this project's own Ensembl-ID DESeq2 output). compare_exp
# also needs to rank the external reference file's gene SYMBOLs, so this
# is a from_type-generalized copy of the same logic (same one-to-many
# bitr collapse: keep the mapped entry with the largest |stat| per
# ENTREZID, then re-sort), not a modification of the original.
#===============================================================================

build_ranked_entrez_vector_generic <- function(ranked_gene_vector, org_pkg, from_type, padj_by_gene = NULL) {

    org_db <- get(org_pkg, envir = asNamespace(org_pkg))

    id_map <- suppressMessages(
        clusterProfiler::bitr(
            names(ranked_gene_vector), fromType = from_type, toType = "ENTREZID", OrgDb = org_db
        )
    )

    if (nrow(id_map) == 0) {
        ranked <- numeric(0)
        if (is.null(padj_by_gene)) return(ranked)
        return(list(ranked = ranked, padj = numeric(0)))
    }

    id_map$stat <- ranked_gene_vector[id_map[[from_type]]]

    id_map <- id_map[order(abs(id_map$stat), decreasing = TRUE), ]

    id_map <- id_map[!duplicated(id_map$ENTREZID), ]

    ranked <- id_map$stat

    names(ranked) <- id_map$ENTREZID

    ranked <- sort(ranked, decreasing = TRUE)

    if (is.null(padj_by_gene)) {
        return(ranked)
    }

    # Same deduplicated (one row per ENTREZID) id_map used to build
    # `ranked` above, so this padj is the one belonging to the exact
    # original gene ID that was kept for each ENTREZID, not an arbitrary
    # one from the one-to-many collapse.
    padj <- padj_by_gene[id_map[[from_type]]]
    names(padj) <- id_map$ENTREZID
    padj <- padj[names(ranked)]

    list(ranked = ranked, padj = padj)

}

#===============================================================================
# Run KEGG GSEA for One Source's Ranked (Shrunken LFC) Gene Vector
#===============================================================================
#
# Always runs dense (pvalueCutoff = 1 via gsea_config$run_gsea_dense =
# TRUE) and splits locally -- one gseKEGG() call yields both the
# significant-only view (for the per-source table/dotplot/summary) and
# the dense view (for pathway correlation), reusing R/gsea.R's own
# split_gsea_result_by_significance() rather than two separate calls.
#===============================================================================

run_gsea_kegg_for_source <- function(results_table, gene_col, lfc_col, from_type, org_pkg, kegg_organism, gsea_config, padj_col = "padj") {

    ranked_gene_vector <- build_ranked_gene_vector(results_table, gene_col = gene_col, stat_col = lfc_col)

    padj_by_gene <- NULL
    if (padj_col %in% names(results_table)) {
        padj_by_gene <- suppressWarnings(as.numeric(results_table[[padj_col]]))
        names(padj_by_gene) <- results_table[[gene_col]]
    }

    entrez_result <- build_ranked_entrez_vector_generic(ranked_gene_vector, org_pkg, from_type, padj_by_gene)

    if (is.list(entrez_result)) {
        ranked_entrez_vector <- entrez_result$ranked
        ranked_entrez_padj <- entrez_result$padj
    } else {
        ranked_entrez_vector <- entrez_result
        ranked_entrez_padj <- NULL
    }

    dense_result <- run_gsea_kegg(ranked_entrez_vector, kegg_organism, gsea_config)

    split_result <- split_gsea_result_by_significance(dense_result, gsea_config$alpha)

    list(
        significant = split_result$significant,
        dense = split_result$dense,
        n_ranked = length(ranked_gene_vector),
        n_mapped_to_entrez = length(ranked_entrez_vector),
        ranked_entrez_vector = ranked_entrez_vector,
        ranked_entrez_padj = ranked_entrez_padj
    )

}

#===============================================================================
# Build the List of GSEA Sources for a Compound (Shrunken LFC
# Everywhere) -- Mirrors get_compare_exp_sources() in
# R/compare_exp_rosi_ora.R
#===============================================================================

get_compare_exp_gsea_sources <- function(
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

        external_results <- read_external_compound_results(external_path, paste0(compound, "_vs_DMSO"))
        external_results <- apply_ashr_shrinkage(external_results, "log2FoldChange", "lfcSE")

        sources[[length(sources) + 1]] <- list(
            source_label = "external", label = paste0(compound, "_vs_DMSO_24h"),
            results_table = external_results,
            gene_col = "gene", padj_col = "padj", from_type = "SYMBOL"
        )

    }

    for (dose in amp0786_doses) {

        f <- file.path(amp0786_dir, paste0("DESeq2_group_", compound, "_", dose, "_vs_DMSO_shrunk.csv"))

        if (!file.exists(f)) {
            warning("Missing AMP0786 shrunk file, skipping: ", f)
            next
        }

        sources[[length(sources) + 1]] <- list(
            source_label = "AMP0786", label = paste0(compound, "_", dose),
            results_table = utils::read.csv(f),
            gene_col = "gene", padj_col = "padj", from_type = "ENSEMBL"
        )

    }

    for (dose in amp0961_doses) {

        f <- file.path(amp0961_dir, paste0("DESeq2_group_", compound, "_", dose, "_vs_Control_shrunk.csv"))

        if (!file.exists(f)) {
            warning("Missing AMP0961 shrunk file, skipping: ", f)
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
# Run the Full Cross-Dataset GSEA KEGG Comparison
#===============================================================================

run_compare_exp_gsea <- function(
    compound = "Rosi",
    amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"),
    amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    include_external = TRUE,
    alpha = 0.05,
    output_dir = file.path("results/compare_exp", compound, sprintf("GSEA_LFC_%s", alpha)),
    figure_width = 8,
    figure_height = 6,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-dataset %s GSEA KEGG comparison...", compound))

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required.")
    }

    org_pkg <- get_organism_db("human")
    kegg_organism <- get_kegg_organism_code("human")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    gsea_config <- make_gsea_config_stub(alpha, figure_width, figure_height, figure_dpi, save_png, save_pdf)

    sources <- get_compare_exp_gsea_sources(compound, amp0786_doses, amp0961_doses, external_path, amp0786_dir, amp0961_dir, include_external)

    summary_rows <- list()

    for (src in sources) {

        message(sprintf("Running GSEA KEGG: %s (%s)...", src$source_label, src$label))

        gsea_result <- run_gsea_kegg_for_source(src$results_table, src$gene_col, "log2FoldChange", src$from_type, org_pkg, kegg_organism, gsea_config)

        filename_base <- paste0("GSEA_KEGG_", src$source_label, "_", src$label)

        save_gsea_result(gsea_result$significant, filename_base, table_dir, figure_dir, gsea_config)

        # Every tested pathway (not just the significant ones) -- no
        # dotplot, a ~300-pathway dotplot isn't readable, and
        # save_gsea_result()'s dotplot is only ever built from the
        # significant-only result above anyway.
        save_gsea_result(gsea_result$dense, paste0(filename_base, "_dense"), table_dir, figure_dir, gsea_config, make_dotplot = FALSE)

        # The ranked (shrunken log2FoldChange, ENTREZID) vector actually
        # fed to gseKEGG() for this source -- same ranking
        # build_ranked_gene_vector()/build_ranked_entrez_vector_generic()
        # produced inside run_gsea_kegg_for_source(), already sorted
        # descending. padj is the original (pre-shrinkage) DESeq2/ashr
        # gene-level adjusted p-value for that same gene, for reference
        # only -- it plays no role in the ranking itself, which is by
        # shrunken LFC.
        ranked_list <- data.frame(
            rank = seq_along(gsea_result$ranked_entrez_vector),
            ENTREZID = names(gsea_result$ranked_entrez_vector),
            shrunken_log2FoldChange = as.numeric(gsea_result$ranked_entrez_vector),
            padj = as.numeric(gsea_result$ranked_entrez_padj[names(gsea_result$ranked_entrez_vector)]),
            stringsAsFactors = FALSE
        )

        save_csv(
            ranked_list,
            file.path(table_dir, paste0(filename_base, "_ranked_gene_list.csv")),
            row.names = FALSE
        )

        n_sig_terms <- if (is.null(gsea_result$significant)) 0 else nrow(as.data.frame(gsea_result$significant))

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
            source = src$source_label,
            label = src$label,
            n_ranked_genes = gsea_result$n_ranked,
            n_mapped_to_entrez = gsea_result$n_mapped_to_entrez,
            n_significant_kegg_terms = n_sig_terms,
            stringsAsFactors = FALSE
        )

    }

    summary_table <- do.call(rbind, summary_rows)

    save_csv(
        summary_table,
        file.path(table_dir, sprintf("compare_exp_%s_gsea_summary.csv", tolower(compound))),
        row.names = FALSE
    )

    message(sprintf("Cross-dataset %s GSEA KEGG comparison complete.", compound))

    summary_table

}

#===============================================================================
# Run Source-vs-Source KEGG Pathway Correlation (Dense NES)
#===============================================================================

run_compare_exp_gsea_pathway_correlation <- function(
    compound = "Rosi",
    amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"),
    amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    include_external = TRUE,
    alpha = 0.05,
    output_dir = file.path("results/compare_exp", compound, sprintf("GSEA_LFC_%s", alpha)),
    figure_width = 8,
    figure_height = 6,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-dataset %s GSEA KEGG pathway correlation (dense NES)...", compound))

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required.")
    }

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    org_pkg <- get_organism_db("human")
    kegg_organism <- get_kegg_organism_code("human")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    gsea_config <- make_gsea_config_stub(alpha, figure_width, figure_height, figure_dpi, save_png, save_pdf)

    sources <- get_compare_exp_gsea_sources(compound, amp0786_doses, amp0961_doses, external_path, amp0786_dir, amp0961_dir, include_external)

    source_names <- vapply(sources, function(s) paste0(s$source_label, "_", s$label), character(1))

    per_source <- lapply(sources, function(src) {

        message(sprintf("Running dense GSEA KEGG: %s (%s)...", src$source_label, src$label))

        gsea_result <- run_gsea_kegg_for_source(src$results_table, src$gene_col, "log2FoldChange", src$from_type, org_pkg, kegg_organism, gsea_config)

        if (is.null(gsea_result$dense)) {
            return(data.frame(ID = character(0), Description = character(0), NES = numeric(0)))
        }

        result_df <- as.data.frame(gsea_result$dense)

        data.frame(ID = result_df$ID, Description = result_df$Description, NES = result_df$NES, stringsAsFactors = FALSE)

    })

    names(per_source) <- source_names

    score_matrix <- Reduce(
        function(a, b) merge(a, b, by = "ID", all = TRUE),
        lapply(source_names, function(nm) {
            df <- per_source[[nm]][, c("ID", "NES")]
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

    save_csv(correlation_matrix, file.path(table_dir, sprintf("compare_exp_%s_gsea_pathway_correlation.csv", tolower(compound))))

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
            main = sprintf("Cross-Dataset %s GSEA KEGG Pathway Correlation (NES, All Tested Pathways)", compound),
            filename = file.path(figure_dir, sprintf("Compare_exp_%s_gsea_pathway_correlation_heatmap.%s", tolower(compound), ext)),
            width = heatmap_width,
            height = heatmap_height
        )

        if (ext == "png") {
            args$res <- figure_dpi
        }

        do.call(pheatmap::pheatmap, args)

    }

    message(sprintf("Cross-dataset %s GSEA KEGG pathway correlation complete.", compound))

    correlation_matrix

}

#===============================================================================
# Run Source-vs-Source Correlation of the Shrunken log2FoldChange for
# Genes Selected (padj < alpha) for GSEA -- Combined/Up/Down
#===============================================================================

run_compare_exp_gsea_gene_correlation <- function(
    compound = "Rosi",
    amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"),
    amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    include_external = TRUE,
    alpha = 0.05,
    output_dir = file.path("results/compare_exp", compound, sprintf("GSEA_LFC_%s", alpha)),
    figure_width = 8,
    figure_height = 6,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-dataset %s GSEA selected-gene shrunken log2FoldChange correlation...", compound))

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

    sources <- get_compare_exp_gsea_sources(compound, amp0786_doses, amp0961_doses, external_path, amp0786_dir, amp0961_dir, include_external)

    source_names <- vapply(sources, function(s) paste0(s$source_label, "_", s$label), character(1))

    per_source <- lapply(sources, function(src) {
        get_gene_lfc_by_entrez(src$results_table, src$gene_col, "log2FoldChange", src$padj_col, alpha, src$from_type, org_db)
    })

    names(per_source) <- source_names

    selected_union <- unique(unlist(lapply(per_source, `[[`, "selected_entrez")))
    selected_up_union <- unique(unlist(lapply(per_source, `[[`, "selected_up_entrez")))
    selected_down_union <- unique(unlist(lapply(per_source, `[[`, "selected_down_entrez")))

    message(sprintf(
        "%d genes selected (padj < %.2f) in at least one of %d sources (%d up, %d down).",
        length(selected_union), alpha, length(source_names), length(selected_up_union), length(selected_down_union)
    ))

    combined <- build_and_save_gene_lfc_correlation(
        per_source, source_names, selected_union, org_db, table_dir, figure_dir,
        compound, "combined", "All Selected Genes (Shrunken LFC)",
        figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    up <- build_and_save_gene_lfc_correlation(
        per_source, source_names, selected_up_union, org_db, table_dir, figure_dir,
        compound, "up", "Up-Regulated (Shrunken LFC)",
        figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    down <- build_and_save_gene_lfc_correlation(
        per_source, source_names, selected_down_union, org_db, table_dir, figure_dir,
        compound, "down", "Down-Regulated (Shrunken LFC)",
        figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    message(sprintf("Cross-dataset %s GSEA selected-gene shrunken log2FoldChange correlation complete.", compound))

    list(combined = combined, up = up, down = down)

}

#===============================================================================
# Scatter Plot: Up- vs. Down-Regulated Shrunken LFC Correlation, One
# Point Per Source Pair
#===============================================================================

run_compare_exp_gsea_updown_scatter <- function(
    compound = "Rosi",
    alpha = 0.05,
    output_dir = file.path("results/compare_exp", compound, sprintf("GSEA_LFC_%s", alpha)),
    figure_width = 9,
    figure_height = 7,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-dataset %s GSEA up- vs. down-correlation scatter plot...", compound))

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
            title = "Up- vs. Down-Regulated Shrunken log2FoldChange Correlation (GSEA), by Source Pair",
            subtitle = paste0(
                sprintf("Each point = one pair of the %d %s sources (external paper + AMP0786 + AMP0961 doses)\n", length(source_names), compound),
                "Dashed line = equal up/down correlation; points below it mean up-regulated genes correlate more strongly"
            ),
            x = "Pearson correlation (shrunken log2FoldChange), up-regulated genes",
            y = "Pearson correlation (shrunken log2FoldChange), down-regulated genes",
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

    message(sprintf("Cross-dataset %s GSEA up- vs. down-correlation scatter plot complete.", compound))

    list(pair_table = pair_table, scatter_plot = scatter_plot)

}

#===============================================================================
# Count of Up- vs. Down-Regulated Genes Selected (padj < alpha), Per
# Source
#===============================================================================

run_compare_exp_gsea_updown_gene_counts <- function(
    compound = "Rosi",
    amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"),
    amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    include_external = TRUE,
    alpha = 0.05,
    output_dir = file.path("results/compare_exp", compound, sprintf("GSEA_LFC_%s", alpha)),
    figure_width = 9,
    figure_height = 6,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-dataset %s GSEA up- vs. down-regulated selected gene counts...", compound))

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

    sources <- get_compare_exp_gsea_sources(compound, amp0786_doses, amp0961_doses, external_path, amp0786_dir, amp0961_dir, include_external)

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
            title = sprintf("Up- vs. Down-Regulated %s Genes Selected for GSEA (Shrunken LFC)", compound),
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

    message(sprintf("Cross-dataset %s GSEA up- vs. down-regulated selected gene counts complete.", compound))

    list(count_table = count_table, bar_plot = bar_plot)

}

#===============================================================================
# Build, Select, Plot, and Save One Direction's (Up or Down) KEGG Dot
# Plot -- Shared by run_compare_exp_gsea_kegg_dotplot()
#===============================================================================

build_and_save_gsea_direction_dotplot <- function(
    dot_table, x_levels, compound, alpha, direction_label, file_tag,
    table_dir, figure_dir, max_pathways, figure_width, figure_height, figure_dpi, save_png, save_pdf
) {

    if (nrow(dot_table) == 0) {
        warning("No ", tolower(direction_label), " pathways found for ", compound, "; skipping that plot.")
        return(NULL)
    }

    n_sources_by_pathway <- tapply(dot_table$x_label, dot_table$Description, function(x) length(unique(x)))
    best_p_by_pathway <- tapply(dot_table$p.adjust, dot_table$Description, min)

    pathway_summary <- data.frame(
        Description = names(n_sources_by_pathway),
        n_sources = as.integer(n_sources_by_pathway),
        best_p = best_p_by_pathway[names(n_sources_by_pathway)],
        stringsAsFactors = FALSE
    )

    pathway_summary <- pathway_summary[order(-pathway_summary$n_sources, pathway_summary$best_p), ]

    selected_pathways <- utils::head(pathway_summary$Description, max_pathways)

    dot_table <- dot_table[dot_table$Description %in% selected_pathways, ]

    dot_table$x_label <- factor(dot_table$x_label, levels = x_levels)
    dot_table$Description <- factor(dot_table$Description, levels = rev(selected_pathways))

    n_pathways_shown <- length(selected_pathways)

    dot_plot <- ggplot2::ggplot(dot_table, ggplot2::aes(x = x_label, y = Description, size = Count, color = p.adjust)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient(low = "red", high = "blue") +
        ggplot2::scale_size_continuous(range = c(2, 9)) +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::labs(
            title = sprintf("%s GSEA KEGG Pathway Enrichment Across Doses and Experiments\n%s", compound, direction_label),
            subtitle = paste0(
                sprintf("Significant, %s pathways (padj < %s, NES %s 0) in at least one source\n", tolower(direction_label), alpha, if (file_tag == "up") ">" else "<"),
                sprintf("Top %d shown, ranked by breadth then significance; size = leading-edge gene count", n_pathways_shown)
            ),
            x = NULL, y = NULL, size = "Leading-edge\ngene count", color = "Adjusted\np-value"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 8),
            axis.text.y = ggplot2::element_text(size = 7),
            panel.grid.major = ggplot2::element_line(color = "grey92")
        )

    save_csv(
        dot_table[, c("x_label", "Description", "Count", "p.adjust", "NES")],
        file.path(table_dir, sprintf("compare_exp_%s_gsea_kegg_dotplot_%s.csv", tolower(compound), file_tag)),
        row.names = FALSE
    )

    plot_height <- max(figure_height, 2 + n_pathways_shown * 0.22)

    if (isTRUE(save_png)) {
        ggplot2::ggsave(
            file.path(figure_dir, sprintf("Compare_exp_%s_gsea_kegg_dotplot_%s.png", tolower(compound), file_tag)),
            plot = dot_plot, width = figure_width, height = plot_height, dpi = figure_dpi, limitsize = FALSE
        )
    }

    if (isTRUE(save_pdf)) {
        ggplot2::ggsave(
            file.path(figure_dir, sprintf("Compare_exp_%s_gsea_kegg_dotplot_%s.pdf", tolower(compound), file_tag)),
            plot = dot_plot, width = figure_width, height = plot_height, limitsize = FALSE
        )
    }

    list(dot_table = dot_table, plot = dot_plot)

}

#===============================================================================
# KEGG Dot Plots: Pathway (y) x Dose/Experiment (x), Size = Leading-Edge
# Gene Count, Color = Adjusted p-Value -- Split Into Up- and
# Down-Regulated Pathways
#===============================================================================
#
# Reads the already-saved, per-source significant-only GSEA KEGG tables
# run_compare_exp_gsea() wrote (GSEA_KEGG_<source_label>_<label>.csv) --
# does not rerun gseKEGG. Mirrors run_compare_exp_kegg_dotplot() in
# R/compare_exp_rosi_ora.R, except dot SIZE is the leading-edge gene
# count (core_enrichment split on "/") in place of ORA's Count column
# (a gseaResult doesn't have one), and -- since GSEA's NES is inherently
# directional, unlike ORA's hypergeometric result -- pathways are split
# into two separate plots by the sign of NES in each (pathway, source)
# cell: NES > 0 ("up", enriched toward the up-regulated end of the
# ranked list) and NES < 0 ("down"). A pathway with a significant hit in
# BOTH directions across different sources/doses (the same kind of
# direction split explored earlier for gene-set membership in
# R/geneset_direction_mix.R) legitimately appears in both plots -- that
# split is itself the point, not an artifact to suppress.
#===============================================================================

run_compare_exp_gsea_kegg_dotplot <- function(
    compound = "Rosi",
    amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"),
    amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"),
    include_external = TRUE,
    alpha = 0.05,
    output_dir = file.path("results/compare_exp", compound, sprintf("GSEA_LFC_%s", alpha)),
    max_pathways = 40,
    figure_width = 10,
    figure_height = 9,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running %s GSEA KEGG dot plots (up/down)...", compound))

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    source_labels <- get_compare_exp_source_labels(compound, amp0786_doses, amp0961_doses, include_external)

    x_levels <- vapply(source_labels, function(s) s$x_label, character(1))

    rows <- list()

    for (src in source_labels) {

        f <- file.path(table_dir, sprintf("GSEA_KEGG_%s_%s.csv", src$source_label, src$label))

        if (!file.exists(f)) {
            warning("Missing GSEA KEGG result file (run run_compare_exp_gsea() first?), skipping: ", f)
            next
        }

        result_df <- utils::read.csv(f)

        if (nrow(result_df) == 0) {
            next
        }

        leading_edge_count <- vapply(
            strsplit(as.character(result_df$core_enrichment), "/", fixed = TRUE),
            function(g) length(unique(g)),
            integer(1)
        )

        rows[[length(rows) + 1]] <- data.frame(
            x_label = src$x_label,
            Description = result_df$Description,
            Count = leading_edge_count,
            p.adjust = result_df$p.adjust,
            NES = result_df$NES,
            stringsAsFactors = FALSE
        )

    }

    if (length(rows) == 0) {
        stop("No significant GSEA KEGG results found for any source of ", compound, " -- nothing to plot.")
    }

    dot_table <- do.call(rbind, rows)

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    up_result <- build_and_save_gsea_direction_dotplot(
        dot_table[dot_table$NES > 0, ], x_levels, compound, alpha, "Up-Regulated", "up",
        table_dir, figure_dir, max_pathways, figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    down_result <- build_and_save_gsea_direction_dotplot(
        dot_table[dot_table$NES < 0, ], x_levels, compound, alpha, "Down-Regulated", "down",
        table_dir, figure_dir, max_pathways, figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    message(sprintf("%s GSEA KEGG dot plots (up/down) complete.", compound))

    list(up = up_result, down = down_result)

}

#===============================================================================
# Build Top-Dose GSEA Sources for Every AMP0961 Compound -- the
# "Different Drug" Null Baseline
#===============================================================================

get_amp0961_all_compound_top_dose_gsea_sources <- function(
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    reference_group = "Control"
) {

    shrunk_files <- list.files(amp0961_dir, pattern = "^DESeq2_group_.*_shrunk\\.csv$")

    group_names <- sub(
        paste0("_vs_", reference_group, "_shrunk\\.csv$"), "",
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
        lapply(split(dose_table, dose_table$compound), function(d) d[which.max(d$dose_nM), ])
    )

    compounds <- sort(top_dose_table$compound)

    lapply(compounds, function(cmpd) {

        group <- top_dose_table$group[top_dose_table$compound == cmpd]
        f <- file.path(amp0961_dir, paste0("DESeq2_group_", group, "_vs_", reference_group, "_shrunk.csv"))

        list(
            source_label = "AMP0961", label = cmpd,
            results_table = utils::read.csv(f),
            gene_col = "gene", padj_col = "padj", from_type = "ENSEMBL"
        )

    })

}

#===============================================================================
# Run the GSEA "Different Drug" Null Baseline: All AMP0961 Compounds
# (Top Dose), Pairwise NES and Shrunken-LFC Correlation
#===============================================================================

run_amp0961_compound_gsea_null_baseline <- function(
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    reference_group = "Control",
    alpha = 0.05,
    output_dir = file.path("results/compare_exp/concordance_proof", sprintf("GSEA_LFC_%s", alpha)),
    figure_width = 10,
    figure_height = 9,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message("Running AMP0961 all-compound top-dose GSEA null baseline...")

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

    gsea_config <- make_gsea_config_stub(alpha, figure_width, figure_height, figure_dpi, save_png, save_pdf)

    sources <- get_amp0961_all_compound_top_dose_gsea_sources(amp0961_dir, reference_group)
    source_names <- vapply(sources, function(s) s$label, character(1))

    message(sprintf("GSEA null baseline across %d AMP0961 compounds (top dose each).", length(sources)))

    per_source_pathway <- lapply(sources, function(src) {

        message(sprintf("Running dense GSEA KEGG: %s...", src$label))

        gsea_result <- run_gsea_kegg_for_source(src$results_table, src$gene_col, "log2FoldChange", src$from_type, org_pkg, kegg_organism, gsea_config)

        if (is.null(gsea_result$dense)) {
            return(data.frame(ID = character(0), NES = numeric(0)))
        }

        result_df <- as.data.frame(gsea_result$dense)

        data.frame(ID = result_df$ID, NES = result_df$NES, stringsAsFactors = FALSE)

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

    save_csv(pathway_correlation_matrix, file.path(table_dir, "amp0961_gsea_null_baseline_pathway_correlation.csv"))

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
            main = "AMP0961 GSEA Null Baseline: All-Compound KEGG Pathway Correlation (NES, Top Dose)",
            filename = file.path(figure_dir, paste0("Amp0961_gsea_null_baseline_pathway_correlation_heatmap.", ext)),
            width = heatmap_width, height = heatmap_height
        )

        if (ext == "png") {
            args$res <- figure_dpi
        }

        do.call(pheatmap::pheatmap, args)

    }

    per_source_genes <- lapply(sources, function(src) {
        get_gene_lfc_by_entrez(src$results_table, src$gene_col, "log2FoldChange", src$padj_col, alpha, src$from_type, org_db)
    })

    names(per_source_genes) <- source_names

    selected_union <- unique(unlist(lapply(per_source_genes, `[[`, "selected_entrez")))

    gene_correlation_matrix <- build_and_save_gene_lfc_correlation(
        per_source_genes, source_names, selected_union, org_db, table_dir, figure_dir,
        "amp0961_gsea_null_baseline", "combined", "All AMP0961 Compounds, Top Dose (GSEA Shrunken LFC)",
        figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    message("AMP0961 all-compound top-dose GSEA null baseline complete.")

    list(
        pathway_correlation_matrix = pathway_correlation_matrix,
        gene_correlation_matrix = gene_correlation_matrix,
        sources = sources
    )

}

#===============================================================================
# GSEA Concordance Proof: Matched Same-Drug Correlations vs. the
# AMP0961 Different-Drug Null Baseline
#===============================================================================

run_compare_exp_gsea_concordance_proof <- function(
    alpha = 0.05,
    output_dir = "results/compare_exp",
    null_baseline_dir = file.path("results/compare_exp/concordance_proof", sprintf("GSEA_LFC_%s", alpha)),
    figure_width = 9,
    figure_height = 7,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message("Running GSEA concordance proof: matched same-drug vs. different-drug null baseline...")

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    table_dir <- file.path(null_baseline_dir, "Tables")
    figure_dir <- file.path(null_baseline_dir, "Figures")

    pathway_null <- as.matrix(utils::read.csv(
        file.path(table_dir, "amp0961_gsea_null_baseline_pathway_correlation.csv"),
        row.names = 1, check.names = FALSE
    ))

    gene_null <- as.matrix(utils::read.csv(
        file.path(table_dir, "compare_exp_amp0961_gsea_null_baseline_selected_gene_lfc_correlation_combined.csv"),
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
        file.path(output_dir, compound, sprintf("GSEA_LFC_%s", alpha))
    }

    read_corr_value <- function(compound, filename, a, b) {
        f <- file.path(compound_output_dir(compound), "Tables", filename)
        m <- as.matrix(utils::read.csv(f, row.names = 1, check.names = FALSE))
        m[a, b]
    }

    matched_pairs <- list(
        list(compound = "Rosi", partner_type = "External", amp0961_label = "AMP0961_Rosi_1000nM", partner_label = "external_Rosi_vs_DMSO_24h"),
        list(compound = "Rosi", partner_type = "AMP0786", amp0961_label = "AMP0961_Rosi_1000nM", partner_label = "AMP0786_Rosi_1000nM"),
        list(compound = "MRL24", partner_type = "External", amp0961_label = "AMP0961_MRL24_100nM", partner_label = "external_MRL24_vs_DMSO_24h"),
        list(compound = "MRL24", partner_type = "AMP0786", amp0961_label = "AMP0961_MRL24_100nM", partner_label = "AMP0786_MRL24_1000nM"),
        list(compound = "GW1929", partner_type = "External", amp0961_label = "AMP0961_GW1929_100nM", partner_label = "external_GW1929_vs_DMSO_24h"),
        list(compound = "T007", partner_type = "AMP0786", amp0961_label = "AMP0961_T007_100nM", partner_label = "AMP0786_T007_10uM")
    )

    rows <- lapply(matched_pairs, function(p) {

        pathway_file <- sprintf("compare_exp_%s_gsea_pathway_correlation.csv", tolower(p$compound))
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
        "Matched same-drug pairs (GSEA): pathway percentile range %.0f-%.0f, gene-LFC percentile range %.0f-%.0f (out of 100 vs. the different-drug null).",
        min(proof_table$pathway_percentile), max(proof_table$pathway_percentile),
        min(proof_table$gene_lfc_percentile), max(proof_table$gene_lfc_percentile)
    ))

    null_long <- rbind(
        data.frame(metric = "KEGG pathway correlation (NES, all tested pathways)", value = pathway_null_values, stringsAsFactors = FALSE),
        data.frame(metric = "Selected-gene shrunken log2FoldChange correlation", value = gene_null_values, stringsAsFactors = FALSE)
    )

    proof_table$pair_label <- paste0(proof_table$compound, " (", proof_table$partner_type, ")")

    matched_long <- rbind(
        data.frame(metric = "KEGG pathway correlation (NES, all tested pathways)", value = proof_table$pathway_r, pair_label = proof_table$pair_label, stringsAsFactors = FALSE),
        data.frame(metric = "Selected-gene shrunken log2FoldChange correlation", value = proof_table$gene_lfc_r, pair_label = proof_table$pair_label, stringsAsFactors = FALSE)
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
            title = "GSEA Concordance Proof: Same Drug Across Experiments vs. Different Drugs Within AMP0961",
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

    message("GSEA concordance proof complete.")

    proof_table

}

#===============================================================================
# Run the Full GSEA compare_exp Suite for Every Compound Already Run
# Through the ORA Suite
#===============================================================================

run_compare_exp_gsea_all <- function(max_pathways = 40) {

    compound_args <- list(
        list(compound = "Rosi", amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"), amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"), include_external = TRUE),
        list(compound = "MRL24", amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"), amp0961_doses = c("0.0125nM", "0.25nM", "5nM", "100nM"), include_external = TRUE),
        list(compound = "GW1929", amp0786_doses = character(0), amp0961_doses = c("0.0125nM", "0.25nM", "5nM", "100nM"), include_external = TRUE),
        list(compound = "T007", amp0786_doses = c("10uM"), amp0961_doses = c("1nM", "10nM", "100nM"), include_external = FALSE)
    )

    for (args in compound_args) {

        do.call(run_compare_exp_gsea, args)
        do.call(run_compare_exp_gsea_pathway_correlation, args)
        do.call(run_compare_exp_gsea_gene_correlation, args)
        run_compare_exp_gsea_updown_scatter(compound = args$compound)
        do.call(run_compare_exp_gsea_updown_gene_counts, args)
        do.call(run_compare_exp_gsea_kegg_dotplot, c(args, list(max_pathways = max_pathways)))

    }

    run_amp0961_compound_gsea_null_baseline()
    run_compare_exp_gsea_concordance_proof()

    message("All compare_exp GSEA analyses complete.")

    invisible(NULL)

}

#===============================================================================
# Combined Cross-Experiment Pathway NES Clustering for the Repeated
# Drugs (Rosi, MRL24, GW1929)
#===============================================================================
#
# Reads the NES column directly from the already-saved, significant-only
# GSEA_KEGG_<source>_<label>.csv files in
# results/compare_exp/<drug>/GSEA_LFC_<alpha>/Tables/ (written by
# run_compare_exp_gsea()) -- does not recompute anything, and does not
# build any combined expression matrix the way the gene/pathway UMAP in
# R/compare_exp_umap.R does. Every one of those files is already a
# single drug-dose-vs-that-experiment's-own-local-control contrast (NES
# is computed by ranking on log2FoldChange relative to each experiment's
# own reference group), so there is no separate "control" row to handle
# here -- a given source simply has no file/row for a group that doesn't
# exist (e.g. AMP0786 has no GW1929 data at all).
#
# Only the significant-only files are read (not a dense recomputation
# like run_compare_exp_gsea_pathway_correlation()'s per-compound
# version), so a pathway/source cell with no significant hit there is
# NA, not 0 -- correlation uses pairwise.complete.obs accordingly. This
# mirrors the same "significant-only correlation is less stable than a
# dense one" caveat already documented for compare_exp's other
# correlation work; flagged here for transparency, not silently hidden.
#===============================================================================

read_saved_gsea_kegg_long <- function(drug, alpha = 0.05, output_dir = "results/compare_exp") {

    table_dir <- file.path(output_dir, drug, sprintf("GSEA_LFC_%s", alpha), "Tables")

    files <- list.files(table_dir, pattern = "^GSEA_KEGG_.*\\.csv$")

    if (length(files) == 0) {
        warning("No saved GSEA_KEGG files found for ", drug, " in ", table_dir)
        return(data.frame(source_group = character(0), ID = character(0), Description = character(0), NES = numeric(0)))
    }

    rows <- lapply(files, function(f) {

        df <- utils::read.csv(file.path(table_dir, f))

        if (nrow(df) == 0) {
            return(NULL)
        }

        label <- sub("^GSEA_KEGG_", "", sub("\\.csv$", "", f))

        data.frame(
            source_group = label,
            ID = df$ID, Description = df$Description, NES = df$NES,
            stringsAsFactors = FALSE
        )

    })

    do.call(rbind, rows[!vapply(rows, is.null, logical(1))])

}

#===============================================================================
# Build the Combined, Dense (All Tested Pathways) NES Matrix for a Set
# of Repeated Drugs -- Shared by the Correlation Heatmap and the UMAP
#===============================================================================
#
# One column per (source, dose-group) -- e.g. "AMP0786_Rosi_1000nM",
# "external_Rosi_vs_DMSO_24h" -- across every drug in `drugs`. Dense =
# every tested pathway, not just the significant ones, needed for a
# non-saturated correlation/embedding (see
# run_compare_exp_repeated_drug_nes_clustering()'s header).
#===============================================================================

build_combined_dense_nes_matrix <- function(
    drugs = c("Rosi", "MRL24", "GW1929"),
    drug_args = list(
        Rosi = list(amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"), amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"), include_external = TRUE),
        MRL24 = list(amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"), amp0961_doses = c("0.0125nM", "0.25nM", "5nM", "100nM"), include_external = TRUE),
        GW1929 = list(amp0786_doses = character(0), amp0961_doses = c("0.0125nM", "0.25nM", "5nM", "100nM"), include_external = TRUE)
    ),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    alpha = 0.05
) {

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required.")
    }

    org_pkg <- get_organism_db("human")
    kegg_organism <- get_kegg_organism_code("human")
    gsea_config <- make_gsea_config_stub(alpha)

    per_source_pathway <- list()

    for (drug in drugs) {

        args <- drug_args[[drug]]

        sources <- get_compare_exp_gsea_sources(
            drug, args$amp0786_doses, args$amp0961_doses, external_path, amp0786_dir, amp0961_dir, args$include_external
        )

        for (src in sources) {

            label <- paste0(src$source_label, "_", src$label)

            message(sprintf("Running dense GSEA KEGG: %s...", label))

            gsea_result <- run_gsea_kegg_for_source(
                src$results_table, src$gene_col, "log2FoldChange", src$from_type, org_pkg, kegg_organism, gsea_config
            )

            if (is.null(gsea_result$dense)) {
                next
            }

            result_df <- as.data.frame(gsea_result$dense)

            per_source_pathway[[label]] <- data.frame(ID = result_df$ID, NES = result_df$NES, stringsAsFactors = FALSE)

        }

    }

    source_names <- names(per_source_pathway)

    message(sprintf("Dense NES computed for %d source_groups.", length(source_names)))

    nes_matrix <- Reduce(
        function(a, b) merge(a, b, by = "ID", all = TRUE),
        lapply(source_names, function(nm) {
            df <- per_source_pathway[[nm]]
            names(df)[2] <- nm
            df
        })
    )

    rownames(nes_matrix) <- nes_matrix$ID
    nes_matrix$ID <- NULL

    as.matrix(nes_matrix[, source_names, drop = FALSE])

}

#===============================================================================
# Parse a source_group Label (e.g. "AMP0786_Rosi_1000nM",
# "external_GW1929_vs_DMSO_24h") Into (source, compound) for Plot
# Color/Shape
#===============================================================================

parse_source_group_label <- function(source_group, drugs) {

    source_label <- ifelse(startsWith(source_group, "external_"), "external", sub("_.*$", "", source_group))

    compound <- vapply(source_group, function(g) {
        hit <- drugs[vapply(drugs, function(d) grepl(paste0("(^|_)", d, "(_|$)"), g), logical(1))]
        if (length(hit) == 0) NA_character_ else hit[1]
    }, character(1))

    data.frame(sample = source_group, source = source_label, compound = unname(compound), stringsAsFactors = FALSE)

}

run_compare_exp_repeated_drug_nes_clustering <- function(
    drugs = c("Rosi", "MRL24", "GW1929"),
    drug_args = list(
        Rosi = list(amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"), amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"), include_external = TRUE),
        MRL24 = list(amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"), amp0961_doses = c("0.0125nM", "0.25nM", "5nM", "100nM"), include_external = TRUE),
        GW1929 = list(amp0786_doses = character(0), amp0961_doses = c("0.0125nM", "0.25nM", "5nM", "100nM"), include_external = TRUE)
    ),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    alpha = 0.05,
    nes_output_dir = "results/compare_exp/NES_comparison",
    figure_width = 12,
    figure_height = 11,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running combined cross-experiment pathway NES clustering (%s, dense)...", paste(drugs, collapse = ", ")))

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    nes_matrix <- build_combined_dense_nes_matrix(drugs, drug_args, external_path, amp0786_dir, amp0961_dir, alpha)

    correlation_matrix <- stats::cor(nes_matrix, method = "pearson", use = "pairwise.complete.obs")

    table_dir <- file.path(nes_output_dir, "Tables")
    figure_dir <- file.path(nes_output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    file_tag <- paste(tolower(drugs), collapse = "_")

    save_csv(nes_matrix, file.path(table_dir, paste0("repeated_drug_nes_matrix_", file_tag, ".csv")))
    save_csv(correlation_matrix, file.path(table_dir, paste0("repeated_drug_nes_correlation_", file_tag, ".csv")))

    can_cluster <- !anyNA(correlation_matrix)

    if (!can_cluster) {
        message("Some source-group pairs share no tested pathways in common (NA correlation); disabling heatmap clustering.")
    }

    heatmap_height <- max(figure_height, ncol(correlation_matrix) * 0.35)
    heatmap_width <- heatmap_height + 2

    for (ext in c("png", "pdf")) {

        save_this_ext <- if (ext == "png") isTRUE(save_png) else isTRUE(save_pdf)

        if (!save_this_ext) {
            next
        }

        args <- list(
            correlation_matrix,
            cluster_rows = can_cluster, cluster_cols = can_cluster,
            display_numbers = TRUE, number_format = "%.2f", fontsize_number = 7,
            main = sprintf(
                "Cross-Experiment Pathway NES Correlation: %s (Dense, All Tested Pathways)",
                paste(drugs, collapse = ", ")
            ),
            filename = file.path(figure_dir, paste0("repeated_drug_nes_correlation_heatmap_", file_tag, ".", ext)),
            width = heatmap_width, height = heatmap_height
        )

        if (ext == "png") {
            args$res <- figure_dpi
        }

        do.call(pheatmap::pheatmap, args)

    }

    message("Combined cross-experiment pathway NES clustering complete.")

    list(nes_matrix = nes_matrix, correlation_matrix = correlation_matrix)

}

#===============================================================================
# Combined Cross-Experiment Pathway UMAP -- NES-Based, Replaces the
# GSVA-Score-Based Pathway UMAP in R/compare_exp_umap.R
#===============================================================================
#
# One point per (source, dose-group) -- NOT per sample, unlike the
# gene-level UMAP in R/compare_exp_umap.R, since NES only exists at the
# contrast level (it's the output of one gseKEGG run on one ranked
# drug-vs-local-control vector, not a per-sample score). Only ~20-25
# points total for 3 drugs, so this is a much smaller embedding than the
# 209-sample gene-level UMAP -- shown for direct comparability with
# run_compare_exp_repeated_drug_nes_clustering()'s heatmap/dendrogram,
# not because a UMAP is the most natural tool for this few points.
#
# Reuses build_combined_dense_nes_matrix() (same dense NES computation
# as the correlation heatmap) and run_umap_with_shape() from
# R/compare_exp_umap.R (source that file first) for the actual
# embedding/plot, so both UMAPs look and save identically.
#===============================================================================

run_compare_exp_repeated_drug_nes_umap <- function(
    drugs = c("Rosi", "MRL24", "GW1929"),
    drug_args = list(
        Rosi = list(amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"), amp0961_doses = c("0.125nM", "2.5nM", "50nM", "1000nM"), include_external = TRUE),
        MRL24 = list(amp0786_doses = c("1nM", "10nM", "100nM", "1000nM"), amp0961_doses = c("0.0125nM", "0.25nM", "5nM", "100nM"), include_external = TRUE),
        GW1929 = list(amp0786_doses = character(0), amp0961_doses = c("0.0125nM", "0.25nM", "5nM", "100nM"), include_external = TRUE)
    ),
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    amp0786_dir = "results/AMP0786/DESeq2/Tables",
    amp0961_dir = "results/AMP0961/DESeq2/Tables",
    alpha = 0.05,
    nes_output_dir = "results/compare_exp/NES_comparison",
    n_neighbors = 5,
    min_dist = 0.1,
    random_seed = 1234,
    figure_width = 9,
    figure_height = 7,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running combined cross-experiment pathway NES UMAP (%s, dense)...", paste(drugs, collapse = ", ")))

    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    nes_matrix <- build_combined_dense_nes_matrix(drugs, drug_args, external_path, amp0786_dir, amp0961_dir, alpha)

    group_metadata <- parse_source_group_label(colnames(nes_matrix), drugs)

    result <- run_umap_with_shape(
        nes_matrix, group_metadata,
        sprintf("Cross-Experiment Pathway UMAP: KEGG NES (%s, Dense)", paste(drugs, collapse = ", ")),
        paste0("repeated_drug_nes_", paste(tolower(drugs), collapse = "_")),
        nes_output_dir, n_neighbors, min_dist, random_seed,
        figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    message("Combined cross-experiment pathway NES UMAP complete.")

    result

}
