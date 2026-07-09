################################################################################
# RNAseqPipeline
#
# File: logging.R
#
# Description:
#   Logging utilities for the RNAseqPipeline.
#
# Author:
#   Andrew Voss
#
################################################################################

#===============================================================================
# Logging Levels
#===============================================================================

LOG_LEVELS <- c(
    "INFO",
    "WARNING",
    "ERROR",
    "SUCCESS",
    "DEBUG"
)

#===============================================================================
#' Create Logger
#'
#' Initialize the log file.
#'
#' @param output_directory Output directory.
#'
#' @return Path to log file.
#===============================================================================

create_logger <- function(output_directory){

    logfile <- file.path(
        output_directory,
        "Logs",
        "analysis.log"
    )

    writeLines(
        c(
            "============================================================",
            "RNAseqPipeline Log",
            paste("Started :", Sys.time()),
            "============================================================",
            ""
        ),
        logfile
    )

    return(logfile)

}

#==========================================================================

.log_message <- function(

    logfile,

    level,

    message,

    verbose = TRUE

){

    stopifnot(level %in% LOG_LEVELS)

    entry <- sprintf(

        "[%s] %-8s %s",

        format(Sys.time(),

               "%Y-%m-%d %H:%M:%S"),

        level,

        message

    )

    if(verbose)

        cat(entry,"\n")

    cat(

        entry,

        "\n",

        file = logfile,

        append = TRUE

    )

}

#============================================================================================

log_info <- function(logfile,

                     message,

                     verbose = TRUE){

    .log_message(

        logfile,

        "INFO",

        message,

        verbose

    )

}

#==========================================================================

log_success <- function(logfile,

                        message,

                        verbose = TRUE){

    .log_message(

        logfile,

        "SUCCESS",

        message,

        verbose

    )

}

#==========================================================================

log_warning <- function(logfile,

                        message,

                        verbose = TRUE){

    .log_message(

        logfile,

        "WARNING",

        message,

        verbose

    )

}

#==========================================================================

log_debug <- function(logfile,

                      message,

                      verbose = TRUE){

    .log_message(

        logfile,

        "DEBUG",

        message,

        verbose

    )

}

#======================================================================================

log_error <- function(

    logfile,

    message,

    stop_pipeline = TRUE,

    verbose = TRUE

){

    .log_message(

        logfile,

        "ERROR",

        message,

        verbose

    )

    if(stop_pipeline)

        stop(message,

             call.=FALSE)

}

#============================================================================

write_session_info <- function(

    output_directory

){

    writeLines(

        capture.output(

            sessionInfo()

        ),

        file.path(

            output_directory,

            "Logs",

            "sessionInfo.txt"

        )

    )

}

#===========================================================================

finish_logger <- function(

    logfile,

    start_time,

    verbose = TRUE

){

    runtime <- round(

        as.numeric(

            difftime(

                Sys.time(),

                start_time,

                units = "mins"

            )

        ),

        2

    )

    log_success(

        logfile,

        paste(

            "Pipeline completed successfully."

        ),

        verbose

    )

    log_info(

        logfile,

        paste(

            "Elapsed time:",

            runtime,

            "minutes"

        ),

        verbose

    )

}

#=================================================================================


