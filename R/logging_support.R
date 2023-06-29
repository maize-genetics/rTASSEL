#' @title Start TASSEL logging information
#'
#' @description This function will create a file for storing logging output
#'    from TASSEL.
#'
#' @param fullPath full working path of log file location. If \code{NULL},
#'    logging file will be added to current working directory.
#' @param fileName name of logging file. If \code{NULL}, filename will resort
#'    to "rTASSEL_log".
#' @param ... additional parameters to be added
#'
#' @export
startLogger <- function(fullPath = NULL, fileName = NULL, ...) {
    if (is.null(fileName)) {
        fileName <- "rTASSEL_log"
    }

    ## Remove slash from end of path (for purely aesthetics purposes only)
    #if (grepl(pattern = "/$", x = fullPath)) {
    #    fullPath <- gsub(pattern = "/$", replacement = "", x = fullPath)
    #}

    if (is.null(fullPath)) {
        file.create(fileName)
        rtlog <- file.path(getwd(), fileName)
    } else {
        if (grepl(pattern = "~", x = fullPath)) {
            stop(
                paste0(
                    "It seems that you are using a '~' instead of your full",
                    " home directory path.\n",
                    "  Consider using: ", Sys.getenv("HOME")
                )
            )
        }
        file.create(file.path(fullPath, fileName))
        rtlog <- file.path(fullPath, fileName)
    }

    rJava::.jcall(
        "net.maizegenetics/util/LoggingUtils",
        "V",
        "setupLogfile",
        rtlog
    )

    message("TASSEL logging file created at: ", rtlog)
}
