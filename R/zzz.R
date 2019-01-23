.onLoad <- function(libname, pkgname) {
    ## Initialize Jar
    rJava::.jpackage(pkgname, lib.loc = libname)

    ## Turn Logging off
    rJava::.jpackage(
        "net.maizegenetics/util/LoggingUtils",
        "V",
        "setupLoggingOff"
    )
}

.onAttach <- function(libname, pkgname){
  msg <- paste0(
    "rTASSEL (experimental)\n",
    "Use with caution!\n"
  )
  packageStartupMessage(msg)
}
