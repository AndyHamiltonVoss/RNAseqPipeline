################################################################################
# RNAseqPipeline
#
# File: compare_exp_umap.R
#
# Description:
#   Cross-experiment sample-level UMAPs (gene expression and KEGG
#   pathway GSVA score), spanning all three compare_exp sources -- the
#   external published reference, AMP0786, and AMP0961 -- restricted to
#   Rosi/MRL24/GW1929/T007 and each source's own control/DMSO, at every
#   dose each source actually tested (not just one matched dose). A
#   compound simply has no points from a source that never tested it
#   (e.g. AMP0786 never tested GW1929, the external paper never tested
#   T007).
#
#   Uses each source's OWN already-computed normalized counts
#   (DESeq2_normalized_counts.csv for AMP0786/AMP0961, the external
#   file's "Normalized count table" sheet) rather than a freshly
#   joint-renormalized raw-count matrix -- there's no statistically
#   meaningful way to compute one shared size-factor normalization
#   across three different pipelines/platforms anyway.
#
#   Since the three sources are independently normalized on different
#   scales, the gene-level UMAP z-scores every gene WITHIN its own
#   source before combining (removes each pipeline's own
#   magnitude/dynamic-range differences, the standard trick for joint
#   embedding of independently-normalized datasets without a formal
#   batch-correction model). The pathway step deliberately does NOT use
#   the z-scored matrix -- GSVA's kcdf step already estimates each
#   gene's own cross-sample distribution internally, so pre-z-scoring
#   would be redundant there, unlike the gene-level UMAP which has no
#   such step of its own.
#
#   Mirrors R/umap.R's single-dataset UMAP shape (top-variable-gene
#   subset -> uwot::umap() -> ggplot2), and reuses
#   build_kegg_gene_sets()/GSVA plumbing from R/gsva.R's run_gsva() for
#   the pathway step. Output in results/compare_exp/UMAP/ -- a
#   cross-compound, cross-source view, not scoped to one compound's
#   subfolder like the rest of compare_exp/.
#
################################################################################

#===============================================================================
# Build One Source's (sample, group) Rows From Its Normalized-Count
# Column Names
#===============================================================================
#
# A group is kept only if it's the source's own control label, or
# parse_dose_group() (R/dose_response.R) successfully parses it AND the
# parsed compound is in `compounds` -- the same exact-match logic used
# throughout this project, which naturally excludes combo arms (e.g.
# "Rosi_1000nM_T007_0.1uM" doesn't match parse_dose_group()'s
# single-compound regex) and every compound not asked for.
#===============================================================================

build_source_sample_rows <- function(cols, groups, source_label, control_group_name, compounds) {

    rows <- list()

    for (i in seq_along(cols)) {

        g <- groups[i]

        if (g == control_group_name) {
            rows[[length(rows) + 1]] <- data.frame(
                original_sample = cols[i], source = source_label, group = "Control", compound = "Control",
                stringsAsFactors = FALSE
            )
            next
        }

        parsed <- parse_dose_group(g)

        if (!is.null(parsed) && parsed$compound %in% compounds) {
            rows[[length(rows) + 1]] <- data.frame(
                original_sample = cols[i], source = source_label, group = g, compound = parsed$compound,
                stringsAsFactors = FALSE
            )
        }

    }

    do.call(rbind, rows)

}

#===============================================================================
# Combined Sample Metadata Across All 3 compare_exp Sources
#===============================================================================

get_combined_umap_sample_metadata <- function(
    amp0786_file = "results/AMP0786/DESeq2/Tables/DESeq2_normalized_counts.csv",
    amp0961_file = "results/AMP0961/DESeq2/Tables/DESeq2_normalized_counts.csv",
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    compounds = c("Rosi", "MRL24", "GW1929", "T007"),
    top_dose_only = FALSE
) {

    if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required.")
    }

    amp0786_cols <- names(utils::read.csv(amp0786_file, nrows = 0, check.names = FALSE))[-1]
    amp0786_groups <- sub("_rep[0-9]+$", "", amp0786_cols)
    amp0786_meta <- build_source_sample_rows(amp0786_cols, amp0786_groups, "AMP0786", "DMSO", compounds)

    amp0961_cols <- names(utils::read.csv(amp0961_file, nrows = 0, check.names = FALSE))[-1]
    amp0961_groups <- sub("_rep-[0-9]+$", "", amp0961_cols)
    amp0961_meta <- build_source_sample_rows(amp0961_cols, amp0961_groups, "AMP0961", "Control", compounds)

    # n_max = 0 would silently drop the blank first (gene-ID) column
    # name that n_max > 0 keeps as "...1" -- read.header with n_max = 1
    # instead so column position [-1] reliably means "drop the gene-ID
    # column", not "drop the first real sample".
    external_header <- suppressMessages(readxl::read_excel(external_path, sheet = "Normalized count table", skip = 2, n_max = 1))
    external_cols <- names(external_header)[-1]
    external_groups_raw <- sub("_rep[0-9]+$", "", external_cols)

    external_compound_map <- c(dmso = "Control", gw1929 = "GW1929", rosi = "Rosi", mrl24 = "MRL24")
    external_compound <- unname(external_compound_map[external_groups_raw])

    keep_ext <- !is.na(external_compound) & (external_compound == "Control" | external_compound %in% compounds)

    external_meta <- data.frame(
        original_sample = external_cols[keep_ext],
        source = "external",
        group = external_compound[keep_ext],
        compound = external_compound[keep_ext],
        stringsAsFactors = FALSE
    )

    combined <- rbind(amp0786_meta, amp0961_meta, external_meta)
    combined$sample <- paste0(combined$source, "__", combined$original_sample)

    # The external source has only one (unlabeled) dose per compound to
    # begin with, so parse_dose_group() returns NULL for it and it's
    # left untouched here regardless of top_dose_only -- only AMP0786/
    # AMP0961 (dose encoded in `group`, e.g. "Rosi_1000nM") get filtered.
    if (isTRUE(top_dose_only)) {

        non_control <- combined$compound != "Control"

        parsed <- lapply(combined$group[non_control], parse_dose_group)

        combined$dose_nM <- NA_real_
        combined$dose_nM[non_control] <- vapply(parsed, function(p) if (is.null(p)) NA_real_ else p$dose_nM, numeric(1))

        max_dose <- stats::ave(
            combined$dose_nM, combined$source, combined$compound,
            FUN = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
        )

        is_top_or_exempt <- combined$compound == "Control" | is.na(combined$dose_nM) | combined$dose_nM == max_dose

        combined <- combined[is_top_or_exempt, ]
        combined$dose_nM <- NULL

    }

    message(sprintf(
        "Combined UMAP metadata: %d samples (%d AMP0786, %d AMP0961, %d external)%s.",
        nrow(combined),
        sum(combined$source == "AMP0786"), sum(combined$source == "AMP0961"), sum(combined$source == "external"),
        if (isTRUE(top_dose_only)) ", top dose only" else ""
    ))

    combined

}

#===============================================================================
# Combined, Gene-Harmonized Expression Matrix (log2 and Within-Source
# z-Scored Versions)
#===============================================================================

get_combined_umap_expression_matrix <- function(
    combined_metadata,
    amp0786_file = "results/AMP0786/DESeq2/Tables/DESeq2_normalized_counts.csv",
    amp0961_file = "results/AMP0961/DESeq2/Tables/DESeq2_normalized_counts.csv",
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx"
) {

    if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required.")
    }

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required.")
    }

    org_pkg <- get_organism_db("human")
    org_db <- get(org_pkg, envir = asNamespace(org_pkg))

    read_normalized <- function(f) {
        df <- utils::read.csv(f, check.names = FALSE, row.names = 1)
        as.matrix(df)
    }

    amp0786_mat <- read_normalized(amp0786_file)
    amp0961_mat <- read_normalized(amp0961_file)

    external_raw <- suppressMessages(readxl::read_excel(external_path, sheet = "Normalized count table", skip = 2))
    names(external_raw)[1] <- "gene_symbol"
    external_raw <- as.data.frame(external_raw)

    symbol_map <- suppressMessages(
        clusterProfiler::bitr(external_raw$gene_symbol, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org_db)
    )
    symbol_map <- symbol_map[!duplicated(symbol_map$SYMBOL), ]
    symbol_map <- symbol_map[!duplicated(symbol_map$ENSEMBL), ]

    external_mat <- as.matrix(external_raw[match(symbol_map$SYMBOL, external_raw$gene_symbol), -1, drop = FALSE])
    rownames(external_mat) <- symbol_map$ENSEMBL

    shared_genes <- Reduce(intersect, list(rownames(amp0786_mat), rownames(amp0961_mat), rownames(external_mat)))

    message(sprintf(
        "Shared genes across all 3 sources: %d (%d AMP0786 / %d AMP0961 / %d external-mapped).",
        length(shared_genes), nrow(amp0786_mat), nrow(amp0961_mat), nrow(external_mat)
    ))

    amp0786_meta <- combined_metadata[combined_metadata$source == "AMP0786", ]
    amp0961_meta <- combined_metadata[combined_metadata$source == "AMP0961", ]
    external_meta <- combined_metadata[combined_metadata$source == "external", ]

    amp0786_sub <- amp0786_mat[shared_genes, amp0786_meta$original_sample, drop = FALSE]
    amp0961_sub <- amp0961_mat[shared_genes, amp0961_meta$original_sample, drop = FALSE]
    external_sub <- external_mat[shared_genes, external_meta$original_sample, drop = FALSE]

    colnames(amp0786_sub) <- amp0786_meta$sample
    colnames(amp0961_sub) <- amp0961_meta$sample
    colnames(external_sub) <- external_meta$sample

    log2_combined <- cbind(log2(amp0786_sub + 1), log2(amp0961_sub + 1), log2(external_sub + 1))
    log2_combined <- log2_combined[, combined_metadata$sample, drop = FALSE]

    z_score_within_source <- function(mat, samples) {
        sub <- mat[, samples, drop = FALSE]
        means <- rowMeans(sub)
        sds <- apply(sub, 1, stats::sd)
        sds[sds == 0] <- NA
        (sub - means) / sds
    }

    z_combined <- cbind(
        z_score_within_source(log2_combined, amp0786_meta$sample),
        z_score_within_source(log2_combined, amp0961_meta$sample),
        z_score_within_source(log2_combined, external_meta$sample)
    )
    z_combined <- z_combined[, combined_metadata$sample, drop = FALSE]

    list(log2_matrix = log2_combined, z_matrix = z_combined, shared_genes = shared_genes)

}

#===============================================================================
# Run UMAP on a (Features x Samples) Matrix and Plot, Colored by
# Compound, Shaped by Source
#===============================================================================

run_umap_with_shape <- function(
    feature_by_sample_matrix, combined_metadata, title, filename_tag, output_dir,
    n_neighbors = 15, min_dist = 0.1, random_seed = 1234,
    figure_width = 9, figure_height = 7, figure_dpi = 600, save_png = TRUE, save_pdf = TRUE
) {

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    mat <- feature_by_sample_matrix[stats::complete.cases(feature_by_sample_matrix), , drop = FALSE]

    n_samples <- ncol(mat)
    n_neighbors_use <- min(n_neighbors, n_samples - 1)

    set.seed(random_seed)

    embedding <- uwot::umap(
        t(mat),
        n_neighbors = n_neighbors_use,
        min_dist = min_dist,
        n_components = 2
    )

    scores <- data.frame(
        sample = colnames(mat),
        UMAP1 = embedding[, 1],
        UMAP2 = embedding[, 2],
        stringsAsFactors = FALSE
    )

    scores <- merge(scores, combined_metadata, by = "sample")

    save_csv(
        scores,
        file.path(table_dir, paste0(filename_tag, "_umap_coordinates.csv")),
        row.names = FALSE
    )

    shape_values <- c(AMP0786 = 16, AMP0961 = 17, external = 15)

    umap_plot <- ggplot2::ggplot(scores, ggplot2::aes(x = UMAP1, y = UMAP2, color = compound, shape = source)) +
        ggplot2::geom_point(size = 3, alpha = 0.85) +
        ggplot2::scale_shape_manual(values = shape_values) +
        ggplot2::labs(title = title, x = "UMAP1", y = "UMAP2", color = "Compound", shape = "Experiment") +
        ggplot2::theme_bw()

    if (isTRUE(save_png)) {
        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_tag, "_umap.png")),
            plot = umap_plot, width = figure_width, height = figure_height, dpi = figure_dpi
        )
    }

    if (isTRUE(save_pdf)) {
        ggplot2::ggsave(
            file.path(figure_dir, paste0(filename_tag, "_umap.pdf")),
            plot = umap_plot, width = figure_width, height = figure_height
        )
    }

    list(scores = scores, plot = umap_plot)

}

#===============================================================================
# Cross-Experiment Sample-Level GENE UMAP
#===============================================================================

run_compare_exp_sample_umap <- function(
    amp0786_file = "results/AMP0786/DESeq2/Tables/DESeq2_normalized_counts.csv",
    amp0961_file = "results/AMP0961/DESeq2/Tables/DESeq2_normalized_counts.csv",
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    compounds = c("Rosi", "MRL24", "GW1929", "T007"),
    include_control = TRUE,
    top_dose_only = FALSE,
    top_n_genes = NULL,
    output_dir = "results/compare_exp/UMAP",
    n_neighbors = 15,
    min_dist = 0.1,
    random_seed = 1234,
    figure_width = 9,
    figure_height = 7,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running cross-experiment sample-level gene UMAP%s...", if (isTRUE(top_dose_only)) " (top dose only)" else ""))

    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    combined_metadata <- get_combined_umap_sample_metadata(amp0786_file, amp0961_file, external_path, compounds, top_dose_only)

    if (!isTRUE(include_control)) {
        combined_metadata <- combined_metadata[combined_metadata$compound != "Control", ]
    }

    expr <- get_combined_umap_expression_matrix(combined_metadata, amp0786_file, amp0961_file, external_path)

    z_mat <- expr$z_matrix[stats::complete.cases(expr$z_matrix), , drop = FALSE]

    # top_n_genes = NULL (default) -> use every shared, complete-case
    # gene, matching the pathway UMAP/the correlation heatmaps (neither
    # of those restrict to a "top variable" subset either). Pass a
    # number to subset to the most variable genes instead, same as
    # before.
    if (is.null(top_n_genes)) {

        embedding_genes <- rownames(z_mat)
        message(sprintf("Embedding %d samples on all %d shared genes (within-source z-scored).", ncol(z_mat), length(embedding_genes)))

    } else {

        gene_vars <- apply(z_mat, 1, stats::var)
        n_top <- min(top_n_genes, length(gene_vars))
        embedding_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_top)]
        message(sprintf("Embedding %d samples on top %d variable genes (within-source z-scored).", ncol(z_mat), n_top))

    }

    default_compounds <- c("Rosi", "MRL24", "GW1929", "T007")
    compound_tag <- if (setequal(compounds, default_compounds)) "" else paste0("_", paste(tolower(sort(compounds)), collapse = "_"))
    control_tag <- if (!isTRUE(include_control)) "_no_control" else ""

    title_suffix <- paste0(
        if (isTRUE(top_dose_only)) " (Top Dose Only)" else "",
        if (nzchar(compound_tag)) sprintf(" [%s]", paste(compounds, collapse = ", ")) else "",
        if (!isTRUE(include_control)) " (No Control)" else ""
    )
    gene_set_label <- if (is.null(top_n_genes)) "All Genes" else "Top Variable Genes"
    file_tag <- paste0(
        "gene", if (is.null(top_n_genes)) "_all" else "", if (isTRUE(top_dose_only)) "_top_dose" else "",
        compound_tag, control_tag
    )

    # The compound-scope suffix can make the title long enough to clip
    # against a fixed-width canvas -- widen a little whenever it's
    # present instead of letting it run off the edge.
    plot_width <- if (nzchar(compound_tag) || nzchar(control_tag)) figure_width + 2 else figure_width

    result <- run_umap_with_shape(
        z_mat[embedding_genes, , drop = FALSE], combined_metadata,
        sprintf("Cross-Experiment Sample UMAP: Gene Expression (%s)%s", gene_set_label, title_suffix),
        file_tag, output_dir, n_neighbors, min_dist, random_seed,
        plot_width, figure_height, figure_dpi, save_png, save_pdf
    )

    message("Cross-experiment sample-level gene UMAP complete.")

    result

}

#===============================================================================
# Cross-Experiment Sample-Level PATHWAY (KEGG GSVA) UMAP
#===============================================================================

#===============================================================================
# Combined, Cross-Experiment KEGG GSVA Score Matrix -- Shared by the
# Pathway UMAP and the Pathway Correlation Heatmap
#===============================================================================

get_combined_gsva_scores <- function(
    combined_metadata, expr,
    gsva_kcdf = "Gaussian",
    gsva_min_gene_set_size = 10,
    gsva_max_gene_set_size = 500,
    table_dir = NULL,
    file_tag = "all"
) {

    if (!requireNamespace("GSVA", quietly = TRUE)) {
        stop("Package 'GSVA' is required.")
    }

    org_pkg <- get_organism_db("human")
    kegg_organism <- get_kegg_organism_code("human")

    gene_sets <- build_kegg_gene_sets(
        kegg_organism = kegg_organism,
        org_pkg = org_pkg,
        universe_gene_ids = rownames(expr$log2_matrix)
    )

    message(sprintf("Scoring %d candidate KEGG pathways across %d samples...", length(gene_sets), ncol(expr$log2_matrix)))

    gsva_param <- GSVA::gsvaParam(
        exprData = expr$log2_matrix,
        geneSets = gene_sets,
        kcdf = gsva_kcdf,
        minSize = gsva_min_gene_set_size,
        maxSize = gsva_max_gene_set_size
    )

    gsva_scores <- GSVA::gsva(gsva_param, verbose = FALSE)

    message(sprintf("%d pathways passed the gene-set size filter.", nrow(gsva_scores)))

    if (!is.null(table_dir)) {
        dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
        save_csv(as.data.frame(gsva_scores), file.path(table_dir, paste0("combined_gsva_scores_", file_tag, ".csv")))
    }

    as.matrix(gsva_scores)

}

run_compare_exp_pathway_umap <- function(
    amp0786_file = "results/AMP0786/DESeq2/Tables/DESeq2_normalized_counts.csv",
    amp0961_file = "results/AMP0961/DESeq2/Tables/DESeq2_normalized_counts.csv",
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    compounds = c("Rosi", "MRL24", "GW1929", "T007"),
    gsva_kcdf = "Gaussian",
    gsva_min_gene_set_size = 10,
    gsva_max_gene_set_size = 500,
    output_dir = "results/compare_exp/UMAP",
    n_neighbors = 15,
    min_dist = 0.1,
    random_seed = 1234,
    figure_width = 9,
    figure_height = 7,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message("Running cross-experiment sample-level pathway (GSVA) UMAP...")

    if (!requireNamespace("GSVA", quietly = TRUE)) {
        stop("Package 'GSVA' is required.")
    }

    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    combined_metadata <- get_combined_umap_sample_metadata(amp0786_file, amp0961_file, external_path, compounds)

    expr <- get_combined_umap_expression_matrix(combined_metadata, amp0786_file, amp0961_file, external_path)

    gsva_scores <- get_combined_gsva_scores(
        combined_metadata, expr, gsva_kcdf, gsva_min_gene_set_size, gsva_max_gene_set_size,
        table_dir = file.path(output_dir, "Tables"), file_tag = paste(tolower(compounds), collapse = "_")
    )

    result <- run_umap_with_shape(
        gsva_scores, combined_metadata,
        "Cross-Experiment Sample UMAP: KEGG Pathway GSVA Score",
        "pathway", output_dir, n_neighbors, min_dist, random_seed,
        figure_width, figure_height, figure_dpi, save_png, save_pdf
    )

    message("Cross-experiment sample-level pathway (GSVA) UMAP complete.")

    result

}

#===============================================================================
# Run Both
#===============================================================================

run_compare_exp_umap_all <- function() {

    gene_result <- run_compare_exp_sample_umap()
    pathway_result <- run_compare_exp_pathway_umap()

    message("Cross-experiment gene and pathway UMAPs complete.")

    list(gene = gene_result, pathway = pathway_result)

}

#===============================================================================
# Combined Cross-Experiment Gene Expression Correlation Heatmap
#===============================================================================
#
# One heatmap covering every (source, dose group) combination for the
# given compounds (default Rosi + MRL24, the two compounds tested in
# all 3 sources) plus each source's own Control/DMSO -- reuses the same
# combined, within-source z-scored expression matrix as
# run_compare_exp_sample_umap() (get_combined_umap_sample_metadata() /
# get_combined_umap_expression_matrix()), averaged to one profile per
# (source, group) over its replicate samples, then Pearson-correlated.
# Unlike the per-compound correlation heatmaps in
# R/compare_exp_rosi_ora.R / R/compare_exp_gsea.R (one compound's doses
# across sources at a time), this puts multiple compounds' dose series
# side by side in ONE matrix, so cross-compound as well as
# cross-experiment structure is visible together.
#===============================================================================

run_compare_exp_combined_correlation_heatmap <- function(
    amp0786_file = "results/AMP0786/DESeq2/Tables/DESeq2_normalized_counts.csv",
    amp0961_file = "results/AMP0961/DESeq2/Tables/DESeq2_normalized_counts.csv",
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    compounds = c("Rosi", "MRL24"),
    output_dir = "results/compare_exp/UMAP",
    figure_width = 10,
    figure_height = 9,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running combined cross-experiment correlation heatmap (%s)...", paste(compounds, collapse = ", ")))

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    combined_metadata <- get_combined_umap_sample_metadata(amp0786_file, amp0961_file, external_path, compounds)

    expr <- get_combined_umap_expression_matrix(combined_metadata, amp0786_file, amp0961_file, external_path)

    z_mat <- expr$z_matrix[stats::complete.cases(expr$z_matrix), , drop = FALSE]

    combined_metadata$source_group <- paste0(combined_metadata$source, "_", combined_metadata$group)

    group_levels <- unique(combined_metadata$source_group)

    group_profile <- vapply(group_levels, function(g) {
        samples <- combined_metadata$sample[combined_metadata$source_group == g]
        rowMeans(z_mat[, samples, drop = FALSE])
    }, numeric(nrow(z_mat)))

    colnames(group_profile) <- group_levels

    correlation_matrix <- stats::cor(group_profile, method = "pearson", use = "pairwise.complete.obs")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    file_tag <- paste(tolower(compounds), collapse = "_")

    save_csv(correlation_matrix, file.path(table_dir, paste0("combined_", file_tag, "_correlation.csv")))

    can_cluster <- !anyNA(correlation_matrix)

    if (!can_cluster) {
        message("Some group pairs share no commonly-tested genes (NA correlation); disabling heatmap clustering.")
    }

    heatmap_height <- max(figure_height, length(group_levels) * 0.35)
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
            main = sprintf("Cross-Experiment Gene Expression Correlation: %s", paste(compounds, collapse = " & ")),
            filename = file.path(figure_dir, paste0("combined_", file_tag, "_correlation_heatmap.", ext)),
            width = heatmap_width, height = heatmap_height
        )

        if (ext == "png") {
            args$res <- figure_dpi
        }

        do.call(pheatmap::pheatmap, args)

    }

    message("Combined cross-experiment correlation heatmap complete.")

    correlation_matrix

}

#===============================================================================
# Combined Cross-Experiment KEGG Pathway (GSVA Score) Correlation
# Heatmap
#===============================================================================
#
# Pathway-level analog of run_compare_exp_combined_correlation_heatmap():
# same (source, dose group) averaging and Pearson correlation, but on
# per-sample KEGG GSVA scores (get_combined_gsva_scores(), shared with
# run_compare_exp_pathway_umap()) instead of z-scored gene expression.
#===============================================================================

run_compare_exp_combined_pathway_correlation_heatmap <- function(
    amp0786_file = "results/AMP0786/DESeq2/Tables/DESeq2_normalized_counts.csv",
    amp0961_file = "results/AMP0961/DESeq2/Tables/DESeq2_normalized_counts.csv",
    external_path = "/mnt/d/DRUGseq/ref/1-s2.0-S0026895X24018340-mmc4-992_supp_datafile2_24h.xlsx",
    compounds = c("Rosi", "MRL24"),
    gsva_kcdf = "Gaussian",
    gsva_min_gene_set_size = 10,
    gsva_max_gene_set_size = 500,
    output_dir = "results/compare_exp/UMAP",
    figure_width = 10,
    figure_height = 9,
    figure_dpi = 600,
    save_png = TRUE,
    save_pdf = TRUE
) {

    message(sprintf("Running combined cross-experiment pathway correlation heatmap (%s)...", paste(compounds, collapse = ", ")))

    if (!requireNamespace("GSVA", quietly = TRUE)) {
        stop("Package 'GSVA' is required.")
    }

    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        stop("Package 'pheatmap' is required.")
    }

    combined_metadata <- get_combined_umap_sample_metadata(amp0786_file, amp0961_file, external_path, compounds)

    expr <- get_combined_umap_expression_matrix(combined_metadata, amp0786_file, amp0961_file, external_path)

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    gsva_scores <- get_combined_gsva_scores(
        combined_metadata, expr, gsva_kcdf, gsva_min_gene_set_size, gsva_max_gene_set_size,
        table_dir = table_dir, file_tag = paste(tolower(compounds), collapse = "_")
    )

    combined_metadata$source_group <- paste0(combined_metadata$source, "_", combined_metadata$group)

    group_levels <- unique(combined_metadata$source_group)

    group_profile <- vapply(group_levels, function(g) {
        samples <- combined_metadata$sample[combined_metadata$source_group == g]
        rowMeans(gsva_scores[, samples, drop = FALSE])
    }, numeric(nrow(gsva_scores)))

    colnames(group_profile) <- group_levels

    correlation_matrix <- stats::cor(group_profile, method = "pearson", use = "pairwise.complete.obs")

    file_tag <- paste(tolower(compounds), collapse = "_")

    save_csv(correlation_matrix, file.path(table_dir, paste0("combined_", file_tag, "_pathway_correlation.csv")))

    can_cluster <- !anyNA(correlation_matrix)

    if (!can_cluster) {
        message("Some group pairs share no commonly-tested pathways (NA correlation); disabling heatmap clustering.")
    }

    heatmap_height <- max(figure_height, length(group_levels) * 0.35)
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
            main = sprintf("Cross-Experiment KEGG Pathway (GSVA) Correlation: %s", paste(compounds, collapse = " & ")),
            filename = file.path(figure_dir, paste0("combined_", file_tag, "_pathway_correlation_heatmap.", ext)),
            width = heatmap_width, height = heatmap_height
        )

        if (ext == "png") {
            args$res <- figure_dpi
        }

        do.call(pheatmap::pheatmap, args)

    }

    message("Combined cross-experiment pathway correlation heatmap complete.")

    correlation_matrix

}
