.onLoad <- function(libname, pkgname) {
    ## Initialize Jar
    rJava::.jpackage(pkgname, lib.loc = libname)
    rJava::.jaddClassPath(dir(file.path(getwd(), "inst/java"), full.names = TRUE))
    rJava::.jcall(
        "net.maizegenetics/util/LoggingUtils",
        "V",
        "setupLoggingOff"
    )
}

.onAttach <- function(libname, pkgname){
    msg <- paste0(
        "Loading rTASSEL - v",
        utils::packageVersion("rTASSEL")
    )
    packageStartupMessage(msg)
}
