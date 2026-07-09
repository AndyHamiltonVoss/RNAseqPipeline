################################################################################
# RNAseqPipeline
#
# File: helpers.R
#
# Description:
#   General helper functions used throughout the RNAseqPipeline.
#
# Author:
#   Andrew Voss
#
################################################################################

#===============================================================================
#' Create Output Directory Structure
#'
#' Creates the complete directory tree required by the pipeline.
#'
#' @param output_directory Character string.
#'
#' @return Invisibly returns the created directories.
#===============================================================================

create_output_directories <- function(output_directory) {

    stopifnot(is.character(output_directory))

    directories <- c(

        output_directory,

        file.path(output_directory, "Logs"),

        file.path(output_directory, "QC"),

        file.path(output_directory, "Figures"),

        file.path(output_directory, "DESeq2"),

        file.path(output_directory, "DESeq2", "Tables"),

        file.path(output_directory, "DESeq2", "Figures"),

        file.path(output_directory, "DESeq2", "Objects"),

        file.path(output_directory, "edgeR"),

        file.path(output_directory, "edgeR", "Tables"),

        file.path(output_directory, "edgeR", "Figures"),

        file.path(output_directory, "edgeR", "Objects"),

        file.path(output_directory, "limma"),

        file.path(output_directory, "limma", "Tables"),

        file.path(output_directory, "limma", "Figures"),

        file.path(output_directory, "limma", "Objects"),

        file.path(output_directory, "Combined"),

        file.path(output_directory, "Metadata"),

        file.path(output_directory, "RDS")

    )

    invisible(lapply(directories, function(x) {

        if (!dir.exists(x))

            dir.create(x, recursive = TRUE)

    }))

}

#===============================================================================
#' Check Required Packages
#'
#' Stops execution if required packages are missing.
#===============================================================================

check_required_packages <- function(packages){

    missing_packages <- packages[

        !sapply(packages,

                requireNamespace,

                quietly = TRUE)

    ]

    if(length(missing_packages)>0){

        stop(

            paste(

                "Missing packages:\n",

                paste(missing_packages,

                      collapse="\n")

            ),

            call.=FALSE

        )

    }

}

assert_file_exists <- function(filename){

    if(!file.exists(filename))

        stop(

            paste("File not found:",filename),

            call.=FALSE

        )

}

#==========================================================================>
#' Relevel Group Factor
#'
#' Ensures a metadata group factor is releveled with the configured
#' reference group as baseline. Stops if the reference group is absent.
#==========================================================================>

relevel_group <- function(metadata, reference_group){

    metadata$group <- factor(metadata$group)

    if (!reference_group %in% levels(metadata$group)) {
        stop(
            paste0(
                "Reference group '",
                reference_group,
                "' was not found in metadata."
            )
        )
    }

    metadata$group <- stats::relevel(
        metadata$group,
        ref = reference_group
    )

    return(metadata)

}

#==========================================================================>
#' Directory Checker
#'
#'
#==========================================================================>

assert_directory_exists <- function(directory){

    if(!dir.exists(directory))

        stop(

            paste("Directory not found:",directory),

            call.=FALSE

        )

}

#==========================================================================>
#' Save CSV
#'
#'
#==========================================================================>

save_csv <- function(data,

                     filename,

                     row.names=TRUE){

    write.csv(

        data,

        filename,

        quote=FALSE,

        row.names=row.names

    )

}

#==========================================================================>
#' Save RDS
#'
#'
#==========================================================================>

save_rds <- function(object,

                     filename){

    saveRDS(

        object,

        filename

    )

}

#==========================================================================>
#' Save RDS
#'
#'
#==========================================================================>

save_plot <- function(plot,

                      filename,

                      width=8,

                      height=6,

                      dpi=600){

    ggplot2::ggsave(

        filename,

        plot=plot,

        width=width,

        height=height,

        dpi=dpi

    )

}

#==========================================================================>
#' Save RDS
#'
#'
#==========================================================================>

timestamp <- function(){

    format(

        Sys.time(),

        "%Y-%m-%d %H:%M:%S"

    )

}

#==========================================================================>
#' Pipeline Banner
#'
#'
#==========================================================================>

pipeline_banner <- function(){

cat(
"
=========================================================
RNAseqPipeline
Differential Expression Analysis Pipeline

Author:
Andrew Voss

Version:
0.1.0

=========================================================
")
}

#==========================================================================>
#' Timer
#'
#'
#==========================================================================>

elapsed_time <- function(start_time){

    difftime(

        Sys.time(),

        start_time,

        units="mins"

    )

}

#==========================================================================>
#' Validate Config
#'
#'
#==========================================================================>

validate_config <- function(config){

    required_fields <- c(

        "counts_file",

        "output_directory",

        "reference_group",

        "threads",

        "alpha"

    )

    missing <- setdiff(

        required_fields,

        names(config)

    )

    if(length(missing)>0){

        stop(

            paste(

                "Missing config fields:\n",

                paste(missing,

                      collapse="\n")

            )

        )

    }

}


