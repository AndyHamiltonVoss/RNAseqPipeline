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
        "^.*_rep",
        "",
        sample_names
    )

    metadata <- data.frame(

        sample = sample_names,

        group = factor(group),

        replicate = as.integer(replicate),

        stringsAsFactors = FALSE

    )

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


