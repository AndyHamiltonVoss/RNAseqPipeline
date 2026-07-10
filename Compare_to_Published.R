################################################################################
# RNAseqPipeline
#
# File: Compare_to_Published.R
#
# Description:
#   Compares this study's DESeq2 results against the published dataset from
#   Rayl et al. 2024 (Molecular Pharmacology), stored in data/external/.
#   Assumes DGE_pipeline.R has already been run (this study's DESeq2 tables
#   must exist in Output/DESeq2/Tables).
################################################################################

source("config.R")

source("R/helpers.R")
source("R/io.R")
source("R/metadata.R")
source("R/deseq2.R")
source("R/dose_response.R")
source("R/venn_diagrams.R")
source("R/enrichment.R")
source("R/external_comparison.R")

check_required_packages(c(
    "readxl",
    "clusterProfiler",
    "VennDiagram",
    "ggplot2"
))

external_data_file <- file.path("data", "external", "Rayl2024_MolPharmacol_supp_datafile2_24h.xlsx")

results_dir <- file.path("analysis", "external_validation")

set.seed(config$random_seed)

reprocessed_results <- run_published_reprocessing(external_data_file, config)

reproducibility_results <- run_reproducibility_check(external_data_file, config, results_dir)

dose_comparison_results <- run_dose_comparison(external_data_file, config, results_dir)

message("External comparison complete.")
