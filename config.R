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
        "/mnt/d/DRUGseq/preliminary/AMP0786_results/AMP0786/counts.umi/counts/SP_L270526AV1_01.umi.counts.sampleIDs.txt",

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

    run_saturation_curve = TRUE,

    run_ercc_cv = TRUE,

    run_ercc_dose_heatmap = TRUE,

    run_de_gene_heatmap = FALSE,

    run_control_plate_correlation = TRUE,

    run_umap = TRUE,

    run_deseq2 = TRUE,

    run_limma = TRUE,

    run_edger = TRUE,

    compare_methods = TRUE,

    run_go = TRUE,

    run_kegg = TRUE,

    run_gsea = FALSE,

    run_gsva = TRUE,


    #---------------------------------------------------------------------------
    # GO / KEGG Enrichment
    #---------------------------------------------------------------------------

    # Gene Ontology sub-ontology for enrichGO: "BP" (biological process),
    # "MF" (molecular function), or "CC" (cellular component).
    go_ontology = "BP",

    # log2FC magnitude cutoff for selecting genes to feed into GO/KEGG
    # enrichment (independent of log2fc_cutoff used for dose-response/Venn
    # comparisons elsewhere).
    go_kegg_log2fc_cutoff = 0.5,

    # Regex; sample groups matching this pattern are skipped for GO/KEGG
    # enrichment (e.g. co-treatment/antagonist arms). Set to NULL or "" to
    # disable.
    exclude_groups_pattern = "T007",


    #---------------------------------------------------------------------------
    # GSVA (Gene Set Variation Analysis)
    #---------------------------------------------------------------------------

    # KEGG pathways smaller than this (after restricting to genes tested by
    # DESeq2) are dropped as unreliable to score.
    gsva_min_gene_set_size = 10,

    # KEGG pathways larger than this are dropped (too generic to be
    # informative as a single per-sample activity score).
    gsva_max_gene_set_size = 500,

    # Kernel used to estimate the per-gene expression CDF: "Gaussian" for
    # continuous, log-scale data (VST-transformed counts here).
    gsva_kcdf = "Gaussian",

    # Number of most-variable pathways to show in the GSVA heatmap.
    gsva_heatmap_top_pathways = 30,

    # Named list of ad hoc pairwise GSVA comparisons (direct group-vs-group
    # limma + volcano plot, not vs the DMSO reference). Each entry is
    # c(group1, group2); the plotted difference is group1 minus group2.
    # Groups must match sample group names exactly.
    gsva_pairwise_comparisons = list(
        # MRL24 saturates PPARg at 10nM; Rosi (lower affinity, ~20nM Kd) needs
        # ~100nM to reach a comparable occupancy level, so this compares
        # pathway activity between the two compounds at their respective
        # roughly-saturating doses.
        MRL24_10nM_vs_Rosi_100nM = c("MRL24_10nM", "Rosi_100nM")
    ),


    #---------------------------------------------------------------------------
    # Method Comparison / Dose-Response
    #---------------------------------------------------------------------------

    # Compounds to plot dose-response curves for. Each must appear in sample
    # group names as "<compound>_<dose><unit>" (e.g. "Rosi_100nM"), with no
    # other groups sharing that exact pattern (co-treatments are excluded
    # automatically since they don't match this pattern).
    dose_response_compounds = c("Rosi", "MRL24"),

    # Named list of custom dose-subset "core" gene sets to also run GO/KEGG
    # enrichment on (genes up/down consistently across just the listed doses),
    # in addition to the full-dose-range core set computed automatically for
    # each compound above. Doses must match group name dose labels (e.g. "100nM").
    core_gene_set_subsets = list(
        # Rosi's binding affinity for PPARg is ~20nM, so 100nM/1000nM are the
        # doses that saturate receptor occupancy; this isolates the pathway
        # signature specific to that saturating-dose regime.
        Rosi_100nM_1000nM = list(compound = "Rosi", doses = c("100nM", "1000nM")),

        # Genes consistently up/down at 100nM/1000nM but NOT already up/down at
        # 1nM/10nM: isolates genes newly recruited once receptor occupancy
        # clears the ~20nM affinity, rather than genes already active below it.
        Rosi_100nM_1000nM_newly_recruited = list(
            compound = "Rosi",
            doses = c("100nM", "1000nM"),
            exclude_doses = c("1nM", "10nM")
        )
    ),


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
    # not the main counts_file/metadata_file above.
    #---------------------------------------------------------------------------

    saturation_reads_file =
        "/mnt/d/DRUGseq/preliminary/AMP0786_results/AMP0786/counts.read/counts/SP_L270526AV1_01.read.counts.sampleIDs.txt",

    saturation_qc_file =
        "/mnt/d/DRUGseq/preliminary/AMP0786_results/AMP0786/qc/SP_L270526AV1_01_matrix.alignment.by.sample.txt",

    # Endogenous gene ID prefix; excludes ERCC spike-in rows from detection counts.
    saturation_gene_prefix = "ENSG",

    # Step size (in reads) of the downsampling depth grid.
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
################################################################################################################################################
