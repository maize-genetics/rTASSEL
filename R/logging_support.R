#' @title Start TASSEL logging information
#'
#' @description This function will create a file for storing logging output
#'    from TASSEL.
#'
#' @param filePath
#' File path and name of log file location. If \code{NULL}, a logging file (\code{rtasel_log}will
#' be added to current working directory.
#'
#' @export
startLogger <- function(filePath = NULL) {
    if (is.null(fileName)) {
        filePath <- "rtassel_log.txt"
    } else {
        # TODO
    }


    rJava::.jcall(
        "net.maizegenetics/util/LoggingUtils",
        "V",
        "setupLogfile",
        rtlog
    )

    message("TASSEL logging file created at: ", rtlog)
}
