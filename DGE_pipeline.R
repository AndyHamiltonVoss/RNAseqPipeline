################################################################################
# RNAseqPipeline
#
# File: DGE_pipeline.R
#
# Description:
#   Main driver script for the RNAseqPipeline.
################################################################################

rm(list = ls())
graphics.off()

#===============================================================================
# Source configuration and functions
#===============================================================================

config_path <- commandArgs(trailingOnly = TRUE)[1]

if (is.na(config_path)) {
    config_path <- "config.R"
}

source(config_path)

source("R/helpers.R")
source("R/io.R")
source("R/metadata.R")
source("R/qc.R")
source("R/saturation.R")
source("R/ercc_qc.R")
source("R/umap.R")
source("R/deseq2.R")
source("R/control_plate_correlation.R")
source("R/edger.R")
source("R/limma.R")
source("R/dose_response.R")
source("R/ercc_dose_heatmap.R")
source("R/de_gene_heatmap.R")
source("R/compound_similarity.R")
source("R/venn_diagrams.R")
source("R/compound_venn.R")
source("R/ligand_group_venn.R")
source("R/dose_concordance.R")
source("R/enrichment.R")
source("R/gsva.R")
source("R/logging.R")

#===============================================================================
# Validate configuration
#===============================================================================

validate_config(config)

#===============================================================================
# Check packages
#===============================================================================

check_required_packages(c(
    "data.table",
    "DESeq2"
))

#===============================================================================
# Set up project
#===============================================================================

start_time <- Sys.time()

set.seed(config$random_seed)

create_output_directories(config$output_directory)

logfile <- create_logger(config$output_directory)

#===============================================================================
# Read and validate counts
#===============================================================================

counts <- read_counts(config$counts_file)

validate_counts(counts)

#===============================================================================
# Create and validate metadata
#===============================================================================
#
# Metadata depends only on sample (column) names, which gene filtering never
# changes, so it is created from the raw counts before filtering.
#===============================================================================

metadata <- create_metadata(
    counts = counts,
    replicate_pattern = config$replicate_pattern
)

validate_metadata(
    metadata = metadata,
    minimum_replicates = config$minimum_samples
)

save_metadata(
    metadata = metadata,
    output_directory = config$output_directory
)

#===============================================================================
# Run Gene Detection Saturation Curve
#===============================================================================
#
# Loads its own non-deduplicated read counts and QC-depth data (config$
# saturation_reads_file / saturation_qc_file) independently of the main
# counts object -- UMI-deduplicated counts discard the duplicate structure
# needed to ask "how many genes would be detected at a given raw read depth".
#===============================================================================

if (isTRUE(config$run_saturation_curve)) {

    saturation_results <- run_saturation_curve(
        metadata = metadata,
        config = config
    )

}

#===============================================================================
# Run ERCC Spike-In %CV by Plate
#===============================================================================
#
# Runs on raw, unfiltered counts so no legitimate low-abundance ERCC probe is
# dropped by the gene-level count filter before being summed.
#===============================================================================

if (isTRUE(config$run_ercc_cv)) {

    ercc_results <- run_ercc_cv(
        counts = counts,
        metadata = metadata,
        config = config
    )

}

#===============================================================================
# Run nUMI/nERCC Ratio Heatmap by Drug x Dose
#===============================================================================

if (isTRUE(config$run_ercc_dose_heatmap)) {

    ercc_dose_heatmap_results <- run_ercc_umi_ratio_heatmap(
        counts = counts,
        metadata = metadata,
        config = config
    )

}

#===============================================================================
# Filter low-count genes
#===============================================================================

counts <- filter_low_counts(
    counts = counts,
    minimum_counts = config$minimum_counts,
    minimum_samples = config$minimum_samples
)

#===============================================================================
# Run QC / PCA
#===============================================================================

if (isTRUE(config$run_qc)) {

    qc_results <- run_qc(
        counts = counts,
        metadata = metadata,
        config = config
    )

}

#===============================================================================
# Run UMAP
#===============================================================================

if (isTRUE(config$run_umap)) {

    umap_results <- run_umap(
        counts = counts,
        metadata = metadata,
        config = config
    )

}

#===============================================================================
# Run DESeq2
#===============================================================================

if (isTRUE(config$run_deseq2)) {

    deseq2_results <- run_deseq2(
        counts = counts,
        metadata = metadata,
        config = config
    )

}

#===============================================================================
# Run Inter-Plate Control-Well Correlation
#===============================================================================

if (isTRUE(config$run_deseq2) && isTRUE(config$run_control_plate_correlation)) {

    control_plate_correlation_results <- run_control_plate_correlation(
        normalized_counts = deseq2_results$normalized_counts,
        metadata = metadata,
        config = config
    )

}

#===============================================================================
# Run DE Gene Count Heatmap by Drug x Dose
#===============================================================================

if (isTRUE(config$run_deseq2) && isTRUE(config$run_de_gene_heatmap)) {

    de_gene_heatmap_results <- run_de_gene_heatmap(
        config = config
    )

}

#===============================================================================
# Run Compound-vs-Compound Gene Similarity
#===============================================================================

if (isTRUE(config$run_deseq2) && isTRUE(config$run_compound_similarity)) {

    compound_similarity_results <- run_compound_similarity(
        config = config
    )

}

#===============================================================================
# Run Ligand-Group Venn Diagrams
#===============================================================================

if (isTRUE(config$run_deseq2) && isTRUE(config$run_ligand_group_venn)) {

    ligand_group_venn_results <- run_ligand_group_venn(
        config = config
    )

}

#===============================================================================
# Run GSVA (KEGG gene sets)
#===============================================================================

if (isTRUE(config$run_gsva)) {

    gsva_results <- run_gsva(
        config = config,
        metadata = metadata
    )

}

#===============================================================================
# Run GO / KEGG enrichment
#===============================================================================

if (isTRUE(config$run_go) || isTRUE(config$run_kegg)) {

    enrichment_results <- run_enrichment(
        config = config,
        metadata = metadata
    )

}

#===============================================================================
# Run edgeR
#===============================================================================

if (isTRUE(config$run_edger)) {

    edger_results <- run_edger(
        counts = counts,
        metadata = metadata,
        config = config
    )

}

#===============================================================================
# Run limma-voom
#===============================================================================

if (isTRUE(config$run_limma)) {

    limma_results <- run_limma(
        counts = counts,
        metadata = metadata,
        config = config
    )

}

#===============================================================================
# Run dose-response comparison
#===============================================================================

if (isTRUE(config$compare_methods)) {

    dose_response_results <- run_dose_response(
        config = config,
        metadata = metadata
    )

}

#===============================================================================
# Run Venn diagram comparison
#===============================================================================

if (isTRUE(config$compare_methods)) {

    venn_results <- run_venn_diagrams(
        config = config,
        metadata = metadata
    )

}

#===============================================================================
# Save session information
#===============================================================================

write_session_info(config$output_directory)

finish_logger(
    logfile = logfile,
    start_time = start_time,
    verbose = config$verbose
)
