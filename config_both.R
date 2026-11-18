################################################################################
# RNAseqPipeline
#
# File: config_both.R
#
# Description:
#   Combined AMP0786/AMP0961 reproducibility check: DESeq2 on the merged
#   raw dataset built by R/prepare_combined_experiment_data.R
#   (Output_both/combined_raw_counts.txt, Output_both/combined_metadata.csv)
#   -- the two drugs shared between both experiments (MRL24, Rosi) at
#   their one exactly-matching dose each (MRL24_100nM, Rosi_1000nM), plus
#   library-matched controls only. design_formula = "~ library + group"
#   adjusts for both the between-experiment batch effect and AMP0961's
#   within-experiment library-to-library variation in one joint fit.
#   metadata$group has 5 levels: "Control" is pooled across BOTH
#   experiments (18 samples), but the two drugs are experiment-specific
#   (MRL24_100nM_AMP0786/MRL24_100nM_AMP0961/Rosi_1000nM_AMP0786/
#   Rosi_1000nM_AMP0961), so DESeq2 produces separate per-experiment DEG
#   results -- each vs. the one shared Control baseline -- from a single
#   joint fit (shared dispersion estimation across all 42 samples),
#   instead of pooling both experiments' drug effect into one combined
#   contrast.
#
#   Pooling Control (but not the drugs) across experiments is what makes
#   "~ library + group" full rank here: making Control experiment-specific
#   too made AMP0786's single library a perfect linear function of its own
#   group levels (confirmed empirically, 9 columns / rank 8 -- see
#   R/prepare_combined_experiment_data.R for the full derivation); this
#   version is rank 8 of 8.
#
#   run_deseq2 and run_umap are enabled; every other module (QC, saturation
#   curve, compound similarity, GO/KEGG, GSVA, etc.) assumes a full
#   single-experiment dataset and isn't meaningful for this 42-sample,
#   four-contrast reproducibility check. UMAP is the exception: it only
#   needs counts/metadata, computes its own VST, and colors by the
#   5-level "group" column, which already encodes both drug identity and
#   experiment -- a useful view of batch structure/reproducibility.
#
#   Run via:
#       Rscript DGE_pipeline.R config_both.R
#
################################################################################

#===============================================================================
# PROJECT SETTINGS
#===============================================================================

config <- list(

    #---------------------------------------------------------------------------
    # General Information
    #---------------------------------------------------------------------------

    project_name = "AMP0786_AMP0961_combined",

    experiment_description =
        "Combined MRL24_100nM / Rosi_1000nM reproducibility check across AMP0786 and AMP0961, batch-corrected for experiment/library",

    organism = "human",

    random_seed = 1234,

    threads = 8,


    #---------------------------------------------------------------------------
    # Input Files
    #---------------------------------------------------------------------------

    counts_file = "Output_both/combined_raw_counts.txt",

    metadata_file = "Output_both/combined_metadata.csv",

    annotation_file = NULL,


    #---------------------------------------------------------------------------
    # Output
    #---------------------------------------------------------------------------

    output_directory = "Output_both",

    overwrite_output = FALSE,

    save_workspace = TRUE,

    save_r_objects = TRUE,


    #---------------------------------------------------------------------------
    # Experimental Design
    #---------------------------------------------------------------------------

    reference_group = "Control",

    design_formula = "~ library + group",

    # Unused when metadata_file is set (load_metadata() reads "replicate"
    # directly instead of parsing it from sample names), kept for
    # completeness only.
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

    # Small design matrix here (8 coefficients: intercept + 3 library +
    # 4 group), unlike AMP0961's full run -- apeglm shrinkage is cheap and
    # enabled.
    lfc_shrinkage = TRUE,

    shrinkage_method = "apeglm",


    #---------------------------------------------------------------------------
    # Pipeline Modules
    #---------------------------------------------------------------------------

    run_qc = FALSE,

    run_saturation_curve = FALSE,

    run_ercc_cv = FALSE,

    run_ercc_dose_heatmap = FALSE,

    run_de_gene_heatmap = FALSE,

    run_compound_similarity = FALSE,

    run_ligand_group_venn = FALSE,

    run_t007_curated_comparison = FALSE,

    run_control_plate_correlation = FALSE,

    run_umap = TRUE,

    # limma::removeBatchEffect(batch = library, design = ~ group) on the
    # same VST matrix, to check whether the batch structure seen in the
    # uncorrected UMAP is a library/experiment artifact.
    run_umap_batch_corrected = TRUE,

    run_deseq2 = TRUE,

    run_limma = FALSE,

    run_edger = FALSE,

    compare_methods = FALSE,

    run_go = FALSE,

    run_kegg = FALSE,

    run_gsea = FALSE,

    run_gsea_compound_correlation = FALSE,

    run_gsva = FALSE,

    run_maya = FALSE,

    run_maya_differential = FALSE,

    run_gsva_compound_correlation = FALSE,

    run_pathway_downregulation_skew = FALSE,

    run_gene_expression_skew = FALSE,

    run_pathway_specificity_comparison = FALSE,

    run_gene_vs_pathway_comparison = FALSE,


    #---------------------------------------------------------------------------
    # Method Comparison / Dose-Response
    #---------------------------------------------------------------------------

    dose_response_compounds = character(0),

    venn_compound_groups = list(),

    venn_cross_group_exclude = character(0),

    t007_co_treatment_suffix = "",

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
    # ERCC Spike-In QC
    #---------------------------------------------------------------------------

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
