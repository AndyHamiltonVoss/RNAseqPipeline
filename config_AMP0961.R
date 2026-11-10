################################################################################
# RNAseqPipeline
#
# File: config_AMP0961.R
#
# Description:
#   Configuration for the AMP0961 BRB-seq compound screen (600 samples, 10
#   libraries/plates, ~20 compounds + Control). Run via:
#
#       Rscript DGE_pipeline.R config_AMP0961.R
#
#   run_qc, run_saturation_curve, run_umap, and run_deseq2 (library-adjusted,
#   ~ library + group) are enabled. edgeR/limma/compare_methods/GO/KEGG/GSVA
#   assume the AMP0786-style two-compound dose-response design
#   (dose_response_compounds, core_gene_set_subsets, etc.) and have not been
#   adapted to AMP0961's much larger compound panel yet, so they stay
#   disabled here pending that work.
#
################################################################################

#===============================================================================
# PROJECT SETTINGS
#===============================================================================

config <- list(

    #---------------------------------------------------------------------------
    # General Information
    #---------------------------------------------------------------------------

    project_name = "AMP0961",

    experiment_description =
        "BRB-seq compound screen (~20 compounds + Control, 10 libraries x 60 samples)",

    organism = "human",      # Options: "human", "mouse"

    random_seed = 1234,

    threads = 8,


    #---------------------------------------------------------------------------
    # Input Files
    #---------------------------------------------------------------------------

    counts_file = file.path(
        "/mnt/d/DRUGseq/experiment/AMP0961_results/AMP0961/counts.umi/counts",
        c(
            "SP_AV1_01_LIB00006625.umi.counts.sampleIDs.txt",
            "SP_AV2_01_LIB00006626.umi.counts.sampleIDs.txt",
            "SP_AV3_01_LIB00006627.umi.counts.sampleIDs.txt",
            "SP_AV4_01_LIB00006628.umi.counts.sampleIDs.txt",
            "SP_AV5_01_LIB00006629.umi.counts.sampleIDs.txt",
            "SP_AV6_01_LIB00006630.umi.counts.sampleIDs.txt",
            "SP_AV7_01_LIB00006631.umi.counts.sampleIDs.txt",
            "SP_AV8_01_LIB00006632.umi.counts.sampleIDs.txt",
            "SP_AV9_01_LIB00006633.umi.counts.sampleIDs.txt",
            "SP_AV10_01_LIB00006634.umi.counts.sampleIDs.txt"
        )
    ),

    metadata_file = NULL,

    annotation_file = NULL,


    #---------------------------------------------------------------------------
    # Output
    #---------------------------------------------------------------------------

    output_directory = "Output_AMP0961",

    overwrite_output = FALSE,

    save_workspace = TRUE,

    save_r_objects = TRUE,


    #---------------------------------------------------------------------------
    # Experimental Design
    #---------------------------------------------------------------------------

    reference_group = "Control",

    design_formula = "~ group",

    replicate_pattern = "_rep-?[0-9]+$",


    #---------------------------------------------------------------------------
    # Gene Filtering
    #---------------------------------------------------------------------------

    minimum_counts = 10,

    minimum_samples = 2,

    remove_low_count_genes = TRUE,


    #---------------------------------------------------------------------------
    # Differential Expression
    #---------------------------------------------------------------------------

    alpha = 0.05,

    log2fc_cutoff = 1,

    # apeglm shrinkage is disabled here: it has to account for the full
    # ~104-column design matrix (10 library + 95 group coefficients) per
    # gene per contrast, which projected to >23 hours across all 94
    # contrasts on this machine. Raw (unshrunken) Wald-test log2FC/p-values
    # are used instead -- significance calls are unaffected, only LFC
    # magnitude stability for low-count genes is lost.
    lfc_shrinkage = FALSE,

    shrinkage_method = "apeglm",


    #---------------------------------------------------------------------------
    # Pipeline Modules
    #---------------------------------------------------------------------------

    run_qc = TRUE,

    run_saturation_curve = TRUE,

    run_ercc_cv = TRUE,

    run_ercc_dose_heatmap = TRUE,

    run_umap = TRUE,

    run_deseq2 = TRUE,

    run_limma = FALSE,

    run_edger = FALSE,

    compare_methods = FALSE,

    run_go = FALSE,

    run_kegg = FALSE,

    run_gsea = FALSE,

    run_gsva = FALSE,


    #---------------------------------------------------------------------------
    # GO / KEGG Enrichment
    #---------------------------------------------------------------------------

    go_ontology = "BP",

    go_kegg_log2fc_cutoff = 0.5,

    exclude_groups_pattern = "T007",


    #---------------------------------------------------------------------------
    # GSVA (Gene Set Variation Analysis)
    #---------------------------------------------------------------------------

    gsva_min_gene_set_size = 10,

    gsva_max_gene_set_size = 500,

    gsva_kcdf = "Gaussian",

    gsva_heatmap_top_pathways = 30,

    gsva_pairwise_comparisons = list(),


    #---------------------------------------------------------------------------
    # Method Comparison / Dose-Response
    #---------------------------------------------------------------------------

    dose_response_compounds = character(0),

    core_gene_set_subsets = list(),


    #---------------------------------------------------------------------------
    # Figure Options
    #---------------------------------------------------------------------------

    save_pdf = TRUE,

    save_png = TRUE,

    figure_width = 8,

    figure_height = 6,

    figure_dpi = 600,


    #---------------------------------------------------------------------------
    # PCA Options
    #---------------------------------------------------------------------------

    pca_top_variable_genes = 500,

    pca_center = TRUE,

    pca_scale = TRUE,


    #---------------------------------------------------------------------------
    # UMAP Options
    #---------------------------------------------------------------------------

    umap_n_neighbors = 15,

    umap_min_dist = 0.1,

    umap_max_legend_groups = 30,


    #---------------------------------------------------------------------------
    # Heatmap Options
    #---------------------------------------------------------------------------

    heatmap_top_genes = 50,

    cluster_rows = TRUE,

    cluster_columns = TRUE,


    #---------------------------------------------------------------------------
    # Gene Detection Saturation Curve
    #---------------------------------------------------------------------------
    #
    # Uses non-deduplicated READ counts (counts.read/, not counts.umi/) against
    # the true total sequencing depth from the provider's per-sample QC table,
    # not the main counts_file above. One path per library, same order as
    # counts_file.
    #---------------------------------------------------------------------------

    saturation_reads_file = file.path(
        "/mnt/d/DRUGseq/experiment/AMP0961_results/AMP0961/counts.read/counts",
        c(
            "SP_AV1_01_LIB00006625.read.counts.sampleIDs.txt",
            "SP_AV2_01_LIB00006626.read.counts.sampleIDs.txt",
            "SP_AV3_01_LIB00006627.read.counts.sampleIDs.txt",
            "SP_AV4_01_LIB00006628.read.counts.sampleIDs.txt",
            "SP_AV5_01_LIB00006629.read.counts.sampleIDs.txt",
            "SP_AV6_01_LIB00006630.read.counts.sampleIDs.txt",
            "SP_AV7_01_LIB00006631.read.counts.sampleIDs.txt",
            "SP_AV8_01_LIB00006632.read.counts.sampleIDs.txt",
            "SP_AV9_01_LIB00006633.read.counts.sampleIDs.txt",
            "SP_AV10_01_LIB00006634.read.counts.sampleIDs.txt"
        )
    ),

    saturation_qc_file = file.path(
        "/mnt/d/DRUGseq/experiment/AMP0961_results/AMP0961/qc",
        c(
            "SP_AV1_01_LIB00006625_matrix.alignment.by.sample.txt",
            "SP_AV2_01_LIB00006626_matrix.alignment.by.sample.txt",
            "SP_AV3_01_LIB00006627_matrix.alignment.by.sample.txt",
            "SP_AV4_01_LIB00006628_matrix.alignment.by.sample.txt",
            "SP_AV5_01_LIB00006629_matrix.alignment.by.sample.txt",
            "SP_AV6_01_LIB00006630_matrix.alignment.by.sample.txt",
            "SP_AV7_01_LIB00006631_matrix.alignment.by.sample.txt",
            "SP_AV8_01_LIB00006632_matrix.alignment.by.sample.txt",
            "SP_AV9_01_LIB00006633_matrix.alignment.by.sample.txt",
            "SP_AV10_01_LIB00006634_matrix.alignment.by.sample.txt"
        )
    ),

    saturation_gene_prefix = "ENSG",

    saturation_step = 100000,


    #---------------------------------------------------------------------------
    # ERCC Spike-In QC
    #---------------------------------------------------------------------------

    # Spike-in row ID prefix in counts_file (the main, UMI-deduplicated matrix).
    ercc_prefix = "ERCC",


    #---------------------------------------------------------------------------
    # Logging
    #---------------------------------------------------------------------------

    verbose = TRUE,

    log_to_file = TRUE

)

################################################################################
# End of configuration
################################################################################
