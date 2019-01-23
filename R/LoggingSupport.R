## Logger

startLogger <- function() {
    file.create("rTASSEL_log")
    rtlog <- paste0(getwd(), "/rTASSEL_log")
    rJava::.jcall(
        "net.maizegenetics/util/LoggingUtils",
        "V",
        "setupLogfile",
        rtlog
    )

    message("Logging file created at: ", rtlog)
}
