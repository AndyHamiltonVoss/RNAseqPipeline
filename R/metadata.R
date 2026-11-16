################################################################################
# RNAseqPipeline
#
# File: metadata.R
#
# Description:
#   Functions for automatically generating and validating sample metadata.
#
# Author:
#   Andrew Voss
#
################################################################################

#===============================================================================
#' Create Sample Metadata
#'
#' Automatically generates a metadata table from sample names.
#'
#' Sample names must end with:
#'     _rep1
#'     _rep2
#'     ...
#'
#===============================================================================

create_metadata <- function(
    counts,
    replicate_pattern = "_rep[0-9]+$"
) {

    sample_names <- colnames(counts)

    group <- sub(
        replicate_pattern,
        "",
        sample_names
    )

    replicate <- sub(
        "^.*_rep-?",
        "",
        sample_names
    )

    metadata <- data.frame(

        sample = sample_names,

        group = factor(group),

        replicate = as.integer(replicate),

        stringsAsFactors = FALSE

    )

    library_map <- attr(counts, "library_map")

    if (!is.null(library_map)) {

        metadata$library <- factor(library_map[sample_names])

    }

    rownames(metadata) <- metadata$sample

    return(metadata)

}

#===============================================================================
#' Load Sample Metadata From an External File
#'
#' Reads a pre-built metadata table (e.g. one assembled outside the normal
#' sample-name-parsing flow, such as a combined multi-experiment dataset)
#' instead of auto-deriving it from counts column names. The file's
#' "sample" column must exactly match the counts matrix's column names.
#===============================================================================

load_metadata <- function(
    metadata_file,
    counts
) {

    assert_file_exists(metadata_file)

    metadata <- utils::read.csv(metadata_file, stringsAsFactors = FALSE)

    if (!"sample" %in% names(metadata)) {
        stop("Metadata file must contain a 'sample' column.")
    }

    missing_from_counts <- setdiff(metadata$sample, colnames(counts))
    missing_from_metadata <- setdiff(colnames(counts), metadata$sample)

    if (length(missing_from_counts) > 0 || length(missing_from_metadata) > 0) {
        stop(
            "Metadata 'sample' column and counts column names must match exactly.\n",
            "In metadata but not counts: ", paste(missing_from_counts, collapse = ", "), "\n",
            "In counts but not metadata: ", paste(missing_from_metadata, collapse = ", ")
        )
    }

    metadata$group <- factor(metadata$group)

    if ("replicate" %in% names(metadata)) {
        metadata$replicate <- as.integer(metadata$replicate)
    }

    if ("library" %in% names(metadata)) {
        metadata$library <- factor(metadata$library)
    }

    rownames(metadata) <- metadata$sample

    return(metadata)

}

#===============================================================================
# Validate Metadata
#===============================================================================

validate_metadata <- function(
    metadata,
    minimum_replicates = 2
) {

    if (!all(c("sample", "group") %in% names(metadata))) {

        stop(
            "Metadata must contain 'sample' and 'group' columns."
        )

    }

    if (any(is.na(metadata$group))) {

        stop(
            "Missing experimental groups detected."
        )

    }

    replicate_counts <- table(metadata$group)

    if (any(replicate_counts < minimum_replicates)) {

        warning(

            "Some groups contain fewer than ",

            minimum_replicates,

            " replicates."

        )

    }

    invisible(TRUE)

}

#===============================================================================
# Get Experimental Groups
#===============================================================================

get_groups <- function(metadata){

    levels(metadata$group)

}

#===============================================================================
# Save Sample Metadata
#===============================================================================

save_metadata <- function(metadata, output_directory) {

    save_csv(
        metadata,
        file.path(output_directory, "Metadata", "sample_metadata.csv"),
        row.names = FALSE
    )

}


