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

source("config.R")

source("R/helpers.R")
source("R/io.R")
source("R/metadata.R")
source("R/deseq2.R")
source("R/edger.R")
source("R/limma.R")
source("R/dose_response.R")
source("R/venn_diagrams.R")
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

counts <- filter_low_counts(
    counts = counts,
    minimum_counts = config$minimum_counts,
    minimum_samples = config$minimum_samples
)

#===============================================================================
# Create and validate metadata
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
