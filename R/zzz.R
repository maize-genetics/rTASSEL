.onLoad <- function(libname, pkgname) {
    ## Initialize Jar
    rJava::.jpackage(pkgname, lib.loc = libname)
    rJava::.jaddClassPath(dir(file.path(getwd(), "inst/java"), full.names = TRUE))
    # rJava::.jcall(
    #     "net.maizegenetics/util/LoggingUtils",
    #     "V",
    #     "setupLoggingOff"
    # )
}

.onAttach <- function(libname, pkgname){
    msg <- paste0(
        "Welcome to rTASSEL (version ", utils::packageVersion("rTASSEL"), ") - PREVIEW", "\n",
        " \u2022 Consider starting a TASSEL log file (see ?startLogger())", "\n"
    )
    packageStartupMessage(msg)
}
