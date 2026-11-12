################################################################################
# RNAseqPipeline
#
# File: generate_ma_plots.R
#
# Description:
#   Regenerates raw-vs-shrunk MA plots from already-saved DESeq2 result
#   tables (_raw.csv / _shrunk.csv), without reloading the dds object or
#   recomputing lfcShrink. Useful when shrinkage was computed separately
#   (e.g. as an offline batch) after the main DESeq2 run already wrote
#   raw-only MA plots.
#
################################################################################

config_path <- commandArgs(trailingOnly = TRUE)[1]

if (is.na(config_path)) {
    config_path <- "config.R"
}

source(config_path)

source("R/deseq2.R")

#===============================================================================
# Build a DESeqResults object from a saved result CSV
#===============================================================================

load_deseq2_results_csv <- function(path) {

    if (!file.exists(path)) {
        return(NULL)
    }

    df <- read.csv(path, stringsAsFactors = FALSE)

    rownames(df) <- df$gene
    df$gene <- NULL

    DESeq2::DESeqResults(S4Vectors::DataFrame(df))
}

#===============================================================================
# Regenerate MA plots for every contrast with a saved _raw.csv
#===============================================================================

table_dir <- file.path(config$output_directory, "DESeq2", "Tables")
figure_dir <- file.path(config$output_directory, "DESeq2", "Figures")

raw_files <- list.files(
    table_dir,
    pattern = "^DESeq2_.*_raw\\.csv$"
)

coef_names <- sub("^DESeq2_", "", sub("_raw\\.csv$", "", raw_files))

message("Found ", length(coef_names), " contrasts. Regenerating MA plots...")

for (coef_name in coef_names) {

    res_raw <- load_deseq2_results_csv(
        file.path(table_dir, paste0("DESeq2_", coef_name, "_raw.csv"))
    )

    res_shrunk <- load_deseq2_results_csv(
        file.path(table_dir, paste0("DESeq2_", coef_name, "_shrunk.csv"))
    )

    save_deseq2_ma_plot(
        res_raw,
        res_shrunk,
        coef_name,
        figure_dir,
        config
    )

    message("Regenerated MA plot: ", coef_name)
}

message("Done.")
