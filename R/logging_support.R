#' @title Start TASSEL logging information
#'
#' @description This function will create a file for storing logging output
#'    from TASSEL.
#'
#' @param filePath
#' File path and name of log file location. If \code{NULL}, a logging file
#' (\code{rtasel_log.txt} will be added to current working directory.
#'
#' @export
startLogger <- function(filePath = NULL, verbose = TRUE) {
    if (is.null(filePath)) {
        filePath <- "rtassel_log.txt"
    } else {
        if (!dir.exists(dirname(filePath))) {
            rlang::stop("No such directory exists provided for logging file")
        }

        filePath <- normalizePath(filePath, mustWork = FALSE)
    }

    jLogUtils <- rJava::.jnew(TASSEL_JVM$LOGGING_UTILS)
    jLogUtils$closeLogfile()
    jLogUtils$setupLogfile(filePath)

    if (verbose) message("TASSEL logging file created at: ", filePath)
}
