.onLoad <- function(libname, pkgname) {
    ## Initialize Jar
    rJava::.jpackage(pkgname, lib.loc = libname)
    rJava::.jaddClassPath(dir(file.path(getwd(), "inst/java"), full.names = TRUE))
}

.onAttach <- function(libname, pkgname){
    msg <- paste0(
        "Welcome to rTASSEL (version ", utils::packageVersion("rTASSEL"), ")", "\n",
        " \u2022 Consider starting a TASSEL log file (see ?startLogger())", "\n"
    )
    packageStartupMessage(msg)
}
