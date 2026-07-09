################################################################################
# RNAseqPipeline
#
# File: io.R
#
# Description:
#   Functions for reading and validating count matrices.
#
################################################################################

#===============================================================================
# Read Count Matrix
#===============================================================================

read_counts <- function(counts_file) {

    assert_file_exists(counts_file)

    message("Reading count matrix...")

    counts <- data.table::fread(
        counts_file,
        data.table = FALSE
    )

    if (ncol(counts) < 2) {
        stop("Count matrix must contain at least one sample column.")
    }

    rownames(counts) <- counts[[1]]

    counts <- counts[, -1, drop = FALSE]

    return(counts)

}

#===============================================================================
# Validate Count Matrix
#===============================================================================

validate_counts <- function(counts) {

    message("Validating count matrix...")

    if (any(is.na(counts))) {

        stop("Count matrix contains missing values.")

    }

    if (any(duplicated(rownames(counts)))) {

        duplicated_genes <- unique(
            rownames(counts)[duplicated(rownames(counts))]
        )

        stop(
            paste(
                "Duplicated gene IDs detected:\n",
                paste(head(duplicated_genes, 20),
                      collapse = "\n")
            )
        )

    }

    if (any(duplicated(colnames(counts)))) {

        stop("Duplicated sample names detected.")

    }

    if (any(counts < 0)) {

        stop("Negative counts detected.")

    }

    if (!all(vapply(counts, is.numeric, logical(1)))) {

        stop("Count matrix contains non-numeric values.")

    }

    invisible(TRUE)

}

#===============================================================================
# Filter Low Count Genes
#===============================================================================

filter_low_counts <- function(
    counts,
    minimum_counts = 10,
    minimum_samples = 2
) {

    keep <- rowSums(counts >= minimum_counts) >= minimum_samples

    filtered <- counts[keep, ]

    message(
        sprintf(
            "Filtered %d low-count genes.",
            sum(!keep)
        )
    )

    return(filtered)

}

#===============================================================================
# Get Sample Names
#===============================================================================

get_sample_names <- function(counts) {

    return(colnames(counts))

}
