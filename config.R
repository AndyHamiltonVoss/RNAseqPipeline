################################################################################
# RNAseqPipeline
#
# File: config.R
#
# Description:
#   User configuration file for the RNAseqPipeline.
#
#   This is the ONLY file users should edit before running the pipeline.
#
# Author:
#   Andrew Voss
#
################################################################################

#===============================================================================
# PROJECT SETTINGS
#===============================================================================

config <- list(

    #---------------------------------------------------------------------------
    # General Information
    #---------------------------------------------------------------------------

    project_name = "AMP0786",

    experiment_description =
        "DRUG-seq analysis of Rosiglitazone and MRL24 treated samples",

    organism = "human",      # Options: "human", "mouse"

    random_seed = 1234,

    threads = 8,


    #---------------------------------------------------------------------------
    # Input Files
    #---------------------------------------------------------------------------

    counts_file =
        "/mnt/d/DRUGseq/AMP0786_results/AMP0786/counts.umi/counts/SP_L270526AV1_01.umi.counts.sampleIDs.txt",

    metadata_file = NULL,

    annotation_file = NULL,


    #---------------------------------------------------------------------------
    # Output
    #---------------------------------------------------------------------------

    output_directory = "Output",

    overwrite_output = FALSE,

    save_workspace = TRUE,

    save_r_objects = TRUE,


    #---------------------------------------------------------------------------
    # Experimental Design
    #---------------------------------------------------------------------------

    reference_group = "DMSO",

    design_formula = "~ group",

    replicate_pattern = "_rep[0-9]+$",


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

    lfc_shrinkage = TRUE,

    shrinkage_method = "apeglm",


    #---------------------------------------------------------------------------
    # Pipeline Modules
    #---------------------------------------------------------------------------

    run_qc = TRUE,

    run_deseq2 = TRUE,

    run_limma = TRUE,

    run_edger = TRUE,

    compare_methods = TRUE,

    run_go = FALSE,

    run_kegg = FALSE,

    run_gsea = FALSE,


    #---------------------------------------------------------------------------
    # Method Comparison / Dose-Response
    #---------------------------------------------------------------------------

    # Compounds to plot dose-response curves for. Each must appear in sample
    # group names as "<compound>_<dose><unit>" (e.g. "Rosi_100nM"), with no
    # other groups sharing that exact pattern (co-treatments are excluded
    # automatically since they don't match this pattern).
    dose_response_compounds = c("Rosi", "MRL24"),


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
    # Heatmap Options
    #---------------------------------------------------------------------------

    heatmap_top_genes = 50,

    cluster_rows = TRUE,

    cluster_columns = TRUE,


    #---------------------------------------------------------------------------
    # Logging
    #---------------------------------------------------------------------------

    verbose = TRUE,

    log_to_file = TRUE

)

################################################################################
# End of configuration
################################################################################################################################################
