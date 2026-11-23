################################################################################
# RNAseqPipeline
#
# File: umap.R
#
# Description:
#   Sample-level UMAP: nonlinear 2D embedding on variance-stabilized counts,
#   same top-variable-gene input as the PCA in qc.R.
#
################################################################################

#===============================================================================
# Run UMAP
#===============================================================================

run_umap <- function(counts, metadata, config) {

    message("Running UMAP...")

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.")
    }

    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    counts <- counts[, rownames(metadata), drop = FALSE]

    vst_data <- tryCatch(
        DESeq2::vst(round(as.matrix(counts)), blind = TRUE),
        error = function(e) {
            warning(
                "DESeq2::vst() failed (", conditionMessage(e), "); ",
                "falling back to log2(counts + 1)."
            )
            log2(as.matrix(counts) + 1)
        }
    )

    gene_vars <- apply(vst_data, 1, stats::var)

    n_top <- min(config$pca_top_variable_genes, nrow(vst_data))

    top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_top)]

    n_samples <- ncol(vst_data)

    n_neighbors <- min(config$umap_n_neighbors, n_samples - 1)

    set.seed(config$random_seed)

    embedding <- uwot::umap(
        t(vst_data[top_genes, , drop = FALSE]),
        n_neighbors = n_neighbors,
        min_dist = config$umap_min_dist,
        n_components = 2
    )

    scores <- data.frame(
        sample = colnames(vst_data),
        UMAP1 = embedding[, 1],
        UMAP2 = embedding[, 2]
    )

    scores <- merge(
        scores,
        metadata,
        by = "sample"
    )

    output_dir <- file.path(config$output_directory, "QC", "UMAP")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        scores,
        file.path(table_dir, "umap_coordinates.csv"),
        row.names = FALSE
    )

    n_groups <- length(unique(scores$group))

    umap_plot <- ggplot2::ggplot(
        scores,
        ggplot2::aes(x = UMAP1, y = UMAP2, color = group)
    ) +
        ggplot2::geom_point(size = 2.5, alpha = 0.85) +
        ggplot2::labs(
            title = "UMAP of Samples",
            x = "UMAP1",
            y = "UMAP2",
            color = "Group"
        ) +
        ggplot2::theme_bw()

    # A legend with too many entries pushes the plot panel itself off the
    # canvas, so it is dropped above a group-count threshold rather than
    # rendering an unreadable (or invisible) plot.
    if (n_groups > config$umap_max_legend_groups) {

        message(
            sprintf(
                "%d groups exceeds umap_max_legend_groups (%d); legend suppressed.",
                n_groups,
                config$umap_max_legend_groups
            )
        )

        umap_plot <- umap_plot + ggplot2::theme(legend.position = "none")

    }

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "UMAP_plot.png"),
            plot = umap_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "UMAP_plot.pdf"),
            plot = umap_plot,
            width = config$figure_width,
            height = config$figure_height
        )

    }

    message("UMAP complete.")

    return(
        list(
            scores = scores,
            plot = umap_plot
        )
    )

}

#===============================================================================
# Run UMAP on limma::removeBatchEffect-Corrected Expression
#===============================================================================
#
# Same VST + top-variable-gene + UMAP pipeline as run_umap(), but with
# limma::removeBatchEffect() applied to the VST matrix first -- batch =
# metadata$library, design = model.matrix(~ group) so the group effect is
# protected rather than regressed out along with the batch effect. Lets a
# multi-experiment/multi-library dataset (e.g. results/compare_both) be visually
# checked for whether the batch structure seen in the uncorrected UMAP
# (run_umap()) is a library/experiment artifact that correction resolves,
# or whether it persists. Requires a metadata "library" column with more
# than one level.
#===============================================================================

run_umap_batch_corrected <- function(counts, metadata, config) {

    message("Running batch-corrected UMAP...")

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.")
    }

    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required.")
    }

    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    if (!"library" %in% names(metadata) || length(unique(metadata$library)) < 2) {
        stop("run_umap_batch_corrected() requires a metadata 'library' column with more than one level.")
    }

    counts <- counts[, rownames(metadata), drop = FALSE]

    vst_data <- tryCatch(
        DESeq2::vst(round(as.matrix(counts)), blind = TRUE),
        error = function(e) {
            warning(
                "DESeq2::vst() failed (", conditionMessage(e), "); ",
                "falling back to log2(counts + 1)."
            )
            log2(as.matrix(counts) + 1)
        }
    )

    design_matrix <- stats::model.matrix(~ group, data = metadata[colnames(vst_data), , drop = FALSE])

    corrected_data <- limma::removeBatchEffect(
        vst_data,
        batch = metadata[colnames(vst_data), "library"],
        design = design_matrix
    )

    gene_vars <- apply(corrected_data, 1, stats::var)

    n_top <- min(config$pca_top_variable_genes, nrow(corrected_data))

    top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_top)]

    n_samples <- ncol(corrected_data)

    n_neighbors <- min(config$umap_n_neighbors, n_samples - 1)

    set.seed(config$random_seed)

    embedding <- uwot::umap(
        t(corrected_data[top_genes, , drop = FALSE]),
        n_neighbors = n_neighbors,
        min_dist = config$umap_min_dist,
        n_components = 2
    )

    scores <- data.frame(
        sample = colnames(corrected_data),
        UMAP1 = embedding[, 1],
        UMAP2 = embedding[, 2]
    )

    scores <- merge(
        scores,
        metadata,
        by = "sample"
    )

    output_dir <- file.path(config$output_directory, "QC", "UMAP_batch_corrected")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        scores,
        file.path(table_dir, "umap_coordinates_batch_corrected.csv"),
        row.names = FALSE
    )

    n_groups <- length(unique(scores$group))

    umap_plot <- ggplot2::ggplot(
        scores,
        ggplot2::aes(x = UMAP1, y = UMAP2, color = group)
    ) +
        ggplot2::geom_point(size = 2.5, alpha = 0.85) +
        ggplot2::labs(
            title = "UMAP of Samples (limma Batch-Corrected)",
            x = "UMAP1",
            y = "UMAP2",
            color = "Group"
        ) +
        ggplot2::theme_bw()

    # A legend with too many entries pushes the plot panel itself off the
    # canvas, so it is dropped above a group-count threshold rather than
    # rendering an unreadable (or invisible) plot.
    if (n_groups > config$umap_max_legend_groups) {

        message(
            sprintf(
                "%d groups exceeds umap_max_legend_groups (%d); legend suppressed.",
                n_groups,
                config$umap_max_legend_groups
            )
        )

        umap_plot <- umap_plot + ggplot2::theme(legend.position = "none")

    }

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "UMAP_plot_batch_corrected.png"),
            plot = umap_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "UMAP_plot_batch_corrected.pdf"),
            plot = umap_plot,
            width = config$figure_width,
            height = config$figure_height
        )

    }

    message("Batch-corrected UMAP complete.")

    return(
        list(
            scores = scores,
            plot = umap_plot,
            corrected_data = corrected_data
        )
    )

}

#===============================================================================
# Run UMAP on Each Compound's Top Tested Dose + Control
#===============================================================================
#
# Same VST + top-variable-gene + UMAP pipeline as run_umap(), but restricted
# first to config$reference_group plus each compound's single highest
# tested dose (parse_dose_group(), same selection as
# save_gsva_top_dose_heatmap() in R/gsva.R) -- drops lower doses and T007
# co-treatment combos (which don't match parse_dose_group()'s pattern) so
# the embedding isn't dominated by near-duplicate low-dose/vehicle-like
# samples.
#===============================================================================

run_umap_top_dose <- function(counts, metadata, config) {

    message("Running top-dose UMAP...")

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.")
    }

    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    groups <- unique(as.character(metadata$group))

    parsed <- lapply(groups, parse_dose_group)

    keep <- vapply(parsed, function(x) !is.null(x), logical(1))

    dose_table <- data.frame(
        group = groups[keep],
        compound = vapply(parsed[keep], `[[`, character(1), "compound"),
        dose_nM = vapply(parsed[keep], `[[`, numeric(1), "dose_nM")
    )

    top_dose_table <- do.call(
        rbind,
        lapply(split(dose_table, dose_table$compound), function(d) {
            d[which.max(d$dose_nM), ]
        })
    )

    selected_groups <- c(config$reference_group, top_dose_table$group)

    metadata <- metadata[metadata$group %in% selected_groups, , drop = FALSE]

    message(sprintf(
        "Top-dose UMAP: %d compounds at their highest tested dose + %s, %d samples total.",
        nrow(top_dose_table), config$reference_group, nrow(metadata)
    ))

    counts <- counts[, rownames(metadata), drop = FALSE]

    vst_data <- tryCatch(
        DESeq2::vst(round(as.matrix(counts)), blind = TRUE),
        error = function(e) {
            warning(
                "DESeq2::vst() failed (", conditionMessage(e), "); ",
                "falling back to log2(counts + 1)."
            )
            log2(as.matrix(counts) + 1)
        }
    )

    gene_vars <- apply(vst_data, 1, stats::var)

    n_top <- min(config$pca_top_variable_genes, nrow(vst_data))

    top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_top)]

    n_samples <- ncol(vst_data)

    n_neighbors <- min(config$umap_n_neighbors, n_samples - 1)

    set.seed(config$random_seed)

    embedding <- uwot::umap(
        t(vst_data[top_genes, , drop = FALSE]),
        n_neighbors = n_neighbors,
        min_dist = config$umap_min_dist,
        n_components = 2
    )

    scores <- data.frame(
        sample = colnames(vst_data),
        UMAP1 = embedding[, 1],
        UMAP2 = embedding[, 2]
    )

    scores <- merge(
        scores,
        metadata,
        by = "sample"
    )

    output_dir <- file.path(config$output_directory, "QC", "UMAP_top_dose")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        scores,
        file.path(table_dir, "umap_coordinates_top_dose.csv"),
        row.names = FALSE
    )

    n_groups <- length(unique(scores$group))

    umap_plot <- ggplot2::ggplot(
        scores,
        ggplot2::aes(x = UMAP1, y = UMAP2, color = group)
    ) +
        ggplot2::geom_point(size = 2.5, alpha = 0.85) +
        ggplot2::labs(
            title = "UMAP of Samples (Top Tested Dose per Compound + Control)",
            x = "UMAP1",
            y = "UMAP2",
            color = "Group"
        ) +
        ggplot2::theme_bw()

    # A legend with too many entries pushes the plot panel itself off the
    # canvas, so it is dropped above a group-count threshold rather than
    # rendering an unreadable (or invisible) plot.
    if (n_groups > config$umap_max_legend_groups) {

        message(
            sprintf(
                "%d groups exceeds umap_max_legend_groups (%d); legend suppressed.",
                n_groups,
                config$umap_max_legend_groups
            )
        )

        umap_plot <- umap_plot + ggplot2::theme(legend.position = "none")

    }

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "UMAP_plot_top_dose.png"),
            plot = umap_plot,
            width = config$figure_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "UMAP_plot_top_dose.pdf"),
            plot = umap_plot,
            width = config$figure_width,
            height = config$figure_height
        )

    }

    message("Top-dose UMAP complete.")

    return(
        list(
            scores = scores,
            plot = umap_plot,
            top_dose_table = top_dose_table
        )
    )

}

#===============================================================================
# Run UMAP on Each Compound's Top Tested Dose + Control, Colored by
# Compound and Shaped by Library
#===============================================================================
#
# Same subset -> VST -> top-variable-gene -> embed pipeline as
# run_umap_top_dose() (duplicated rather than refactored out, matching
# this file's existing one-function-per-view convention), but plotted
# with color = compound (not the full "<compound>_<dose>" group label)
# and shape = library, so a compound's samples split across its two
# plates are visually distinguishable from each other -- a per-sample
# check for whether library/plate identity creates structure alongside
# or instead of compound identity. Output lands in the SAME directory as
# run_umap_top_dose() (same sample set, different aesthetic), under
# distinct filenames.
#===============================================================================

run_umap_top_dose_by_library <- function(counts, metadata, config) {

    message("Running top-dose UMAP by compound/library...")

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.")
    }

    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    groups <- unique(as.character(metadata$group))

    parsed <- lapply(groups, parse_dose_group)

    keep <- vapply(parsed, function(x) !is.null(x), logical(1))

    dose_table <- data.frame(
        group = groups[keep],
        compound = vapply(parsed[keep], `[[`, character(1), "compound"),
        dose_nM = vapply(parsed[keep], `[[`, numeric(1), "dose_nM")
    )

    top_dose_table <- do.call(
        rbind,
        lapply(split(dose_table, dose_table$compound), function(d) {
            d[which.max(d$dose_nM), ]
        })
    )

    selected_groups <- c(config$reference_group, top_dose_table$group)

    metadata <- metadata[metadata$group %in% selected_groups, , drop = FALSE]

    message(sprintf(
        "Top-dose UMAP by compound/library: %d compounds at their highest tested dose + %s, %d samples total.",
        nrow(top_dose_table), config$reference_group, nrow(metadata)
    ))

    counts <- counts[, rownames(metadata), drop = FALSE]

    vst_data <- tryCatch(
        DESeq2::vst(round(as.matrix(counts)), blind = TRUE),
        error = function(e) {
            warning(
                "DESeq2::vst() failed (", conditionMessage(e), "); ",
                "falling back to log2(counts + 1)."
            )
            log2(as.matrix(counts) + 1)
        }
    )

    gene_vars <- apply(vst_data, 1, stats::var)

    n_top <- min(config$pca_top_variable_genes, nrow(vst_data))

    top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_top)]

    n_samples <- ncol(vst_data)

    n_neighbors <- min(config$umap_n_neighbors, n_samples - 1)

    set.seed(config$random_seed)

    embedding <- uwot::umap(
        t(vst_data[top_genes, , drop = FALSE]),
        n_neighbors = n_neighbors,
        min_dist = config$umap_min_dist,
        n_components = 2
    )

    scores <- data.frame(
        sample = colnames(vst_data),
        UMAP1 = embedding[, 1],
        UMAP2 = embedding[, 2]
    )

    scores <- merge(
        scores,
        metadata,
        by = "sample"
    )

    # "Control" doesn't match parse_dose_group()'s "<compound>_<dose>"
    # regex, so it comes back NA here rather than a real compound name.
    scores$compound <- vapply(
        scores$group,
        function(g) {
            parsed_group <- parse_dose_group(as.character(g))
            if (is.null(parsed_group)) as.character(g) else parsed_group$compound
        },
        character(1)
    )

    output_dir <- file.path(config$output_directory, "QC", "UMAP_top_dose")

    table_dir <- file.path(output_dir, "Tables")
    figure_dir <- file.path(output_dir, "Figures")

    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    save_csv(
        scores,
        file.path(table_dir, "umap_coordinates_top_dose_by_library.csv"),
        row.names = FALSE
    )

    compounds <- sort(unique(scores$compound))

    base_colors <- c(Rosi = "#1b9e77", MRL24 = "#d95f02")

    unmatched_compounds <- setdiff(compounds, names(base_colors))

    if (length(unmatched_compounds) > 0) {

        extra_palette <- grDevices::colorRampPalette(
            c("#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#1f78b4", "#33a02c", "#e31a1c")
        )(length(unmatched_compounds))

        compound_colors <- c(base_colors, stats::setNames(extra_palette, unmatched_compounds))

    } else {

        compound_colors <- base_colors

    }

    compound_colors <- compound_colors[compounds]

    libraries <- sort(unique(scores$library))

    # ggplot2's default shape scale only supports 6 distinct values
    # before recycling (with a warning); up to 10 libraries can appear
    # here since Control alone spans all 10.
    shape_palette <- c(16, 17, 15, 18, 3, 4, 8, 0, 1, 2, 5, 6)
    library_shapes <- stats::setNames(shape_palette[seq_along(libraries)], libraries)

    umap_plot <- ggplot2::ggplot(
        scores,
        ggplot2::aes(x = UMAP1, y = UMAP2, color = compound, shape = library)
    ) +
        ggplot2::geom_point(size = 2.5, alpha = 0.85) +
        ggplot2::scale_color_manual(values = compound_colors) +
        ggplot2::scale_shape_manual(values = library_shapes) +
        ggplot2::labs(
            title = "UMAP of Samples (Top Tested Dose per Compound + Control)",
            x = "UMAP1",
            y = "UMAP2",
            color = "Compound",
            shape = "Library"
        ) +
        ggplot2::theme_bw()

    legend_width <- config$figure_width + 3

    if (isTRUE(config$save_png)) {

        ggplot2::ggsave(
            file.path(figure_dir, "UMAP_plot_top_dose_by_library.png"),
            plot = umap_plot,
            width = legend_width,
            height = config$figure_height,
            dpi = config$figure_dpi
        )

    }

    if (isTRUE(config$save_pdf)) {

        ggplot2::ggsave(
            file.path(figure_dir, "UMAP_plot_top_dose_by_library.pdf"),
            plot = umap_plot,
            width = legend_width,
            height = config$figure_height
        )

    }

    message("Top-dose UMAP by compound/library complete.")

    return(
        list(
            scores = scores,
            plot = umap_plot,
            top_dose_table = top_dose_table
        )
    )

}
