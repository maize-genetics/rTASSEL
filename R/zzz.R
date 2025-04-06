.onLoad <- function(libname, pkgname) {
    # Initialize Jar
    rJava::.jpackage(pkgname, lib.loc = libname)
    rJava::.jaddClassPath(dir(file.path(getwd(), "inst/java"), full.names = TRUE))

    # Start up a temp logging file regardless of user input
    # NOTE: users can define a logging file whenever they want - this is to
    #       prevent initial TASSEL stdout from contaminating the console space
    #       if a user has not defined one on initial use
    startLogger(tempfile(), verbose = FALSE)
}

.onAttach <- function(libname, pkgname){
    msg <- paste0(
        "Welcome to rTASSEL (version ", utils::packageVersion("rTASSEL"), ")", "\n",
        " \u2022 Consider starting a TASSEL log file (see ?startLogger())", "\n"
    )
    packageStartupMessage(msg)
}
