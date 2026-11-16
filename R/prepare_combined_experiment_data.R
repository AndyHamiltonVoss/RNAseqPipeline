################################################################################
# RNAseqPipeline
#
# File: prepare_combined_experiment_data.R
#
# Description:
#   Data-preparation step (not a DESeq2 run): builds one combined raw
#   count matrix and metadata table spanning both AMP0786 (pilot, single
#   library, reference group "DMSO") and AMP0961 (10-library screen,
#   reference group "Control"), restricted to the two drugs shared
#   between the experiments (MRL24, Rosi) at their one exactly-shared
#   dose each (MRL24_100nM, Rosi_1000nM), plus library-matched controls
#   only. The output is meant to feed a future combined DESeq2 run with a
#   "~ library + group" design (the same batch-adjustment pattern already
#   used for AMP0961's own 10-library data, R/deseq2.R:33-43) to separate
#   the drug effect from the between-experiment batch effect directly,
#   rather than comparing two independently-fit models post-hoc.
#
################################################################################

#===============================================================================
# Build Combined Raw Count Matrix and Metadata
#===============================================================================

build_combined_experiment_dataset <- function() {

    message("Building combined AMP0786/AMP0961 raw dataset...")

    output_dir <- "Output_both"
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    shared_drug_groups <- c("MRL24_100nM", "Rosi_1000nM")

    #===========================================================================
    # AMP0786: All 18 Samples (One Library, So Controls Need No Filtering)
    #===========================================================================

    amp0786_metadata <- read.csv("Output/Metadata/sample_metadata.csv", stringsAsFactors = FALSE)
    amp0786_metadata <- amp0786_metadata[amp0786_metadata$group %in% c("DMSO", shared_drug_groups), ]

    amp0786_metadata$library <- "SP_L270526AV1_01"
    amp0786_metadata$experiment <- "AMP0786"
    amp0786_metadata$original_group <- amp0786_metadata$group
    amp0786_metadata$group[amp0786_metadata$group == "DMSO"] <- "Control"
    amp0786_metadata$original_sample_id <- amp0786_metadata$sample
    amp0786_metadata$sample <- paste0("AMP0786__", amp0786_metadata$original_sample_id)

    #===========================================================================
    # AMP0961: Drug Samples + Only Library-Matched Controls
    #===========================================================================

    amp0961_metadata <- read.csv("Output_AMP0961/Metadata/sample_metadata.csv", stringsAsFactors = FALSE)

    drug_libraries <- unique(amp0961_metadata$library[amp0961_metadata$group %in% shared_drug_groups])

    message(sprintf(
        "AMP0961 libraries containing %s: %s",
        paste(shared_drug_groups, collapse = "/"), paste(drug_libraries, collapse = ", ")
    ))

    amp0961_keep <- amp0961_metadata$group %in% shared_drug_groups |
        (amp0961_metadata$group == "Control" & amp0961_metadata$library %in% drug_libraries)

    amp0961_metadata <- amp0961_metadata[amp0961_keep, ]

    amp0961_metadata$experiment <- "AMP0961"
    amp0961_metadata$original_group <- amp0961_metadata$group
    amp0961_metadata$original_sample_id <- amp0961_metadata$sample
    amp0961_metadata$sample <- paste0("AMP0961__", amp0961_metadata$original_sample_id)

    #===========================================================================
    # Combine Metadata
    #===========================================================================

    metadata_columns <- c("sample", "original_sample_id", "group", "original_group", "replicate", "library", "experiment")

    combined_metadata <- rbind(
        amp0786_metadata[, metadata_columns],
        amp0961_metadata[, metadata_columns]
    )

    message(sprintf(
        "Combined metadata: %d samples (%d AMP0786, %d AMP0961).",
        nrow(combined_metadata), sum(combined_metadata$experiment == "AMP0786"), sum(combined_metadata$experiment == "AMP0961")
    ))

    print(table(combined_metadata$group, combined_metadata$experiment))
    print(table(combined_metadata$library, combined_metadata$experiment))

    #===========================================================================
    # Load and Subset Raw Counts
    #===========================================================================

    amp0786_counts_file <- "/mnt/d/DRUGseq/preliminary/AMP0786_results/AMP0786/counts.umi/counts/SP_L270526AV1_01.umi.counts.sampleIDs.txt"

    amp0961_counts_file <- file.path(
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
    )

    amp0786_counts <- read_counts(amp0786_counts_file)
    amp0961_counts <- read_counts(amp0961_counts_file)

    amp0786_counts <- amp0786_counts[, amp0786_metadata$original_sample_id, drop = FALSE]
    amp0961_counts <- amp0961_counts[, amp0961_metadata$original_sample_id, drop = FALSE]

    shared_genes <- intersect(rownames(amp0786_counts), rownames(amp0961_counts))

    message(sprintf(
        "Gene ID overlap: %d shared / %d AMP0786 / %d AMP0961 (%.1f%% of the smaller set).",
        length(shared_genes), nrow(amp0786_counts), nrow(amp0961_counts),
        100 * length(shared_genes) / min(nrow(amp0786_counts), nrow(amp0961_counts))
    ))

    amp0786_counts <- amp0786_counts[shared_genes, , drop = FALSE]
    amp0961_counts <- amp0961_counts[shared_genes, , drop = FALSE]

    colnames(amp0786_counts) <- paste0("AMP0786__", colnames(amp0786_counts))
    colnames(amp0961_counts) <- paste0("AMP0961__", colnames(amp0961_counts))

    combined_counts <- cbind(amp0786_counts, amp0961_counts)
    combined_counts <- combined_counts[, combined_metadata$sample, drop = FALSE]

    #===========================================================================
    # Save Outputs
    #===========================================================================

    write.csv(combined_metadata, file.path(output_dir, "combined_metadata.csv"), row.names = FALSE)

    combined_counts_out <- data.frame(Gene_id = rownames(combined_counts), combined_counts, check.names = FALSE)

    utils::write.table(
        combined_counts_out,
        file.path(output_dir, "combined_raw_counts.txt"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    message("Combined dataset written to Output_both/.")

    return(
        list(
            metadata = combined_metadata,
            counts = combined_counts
        )
    )

}
