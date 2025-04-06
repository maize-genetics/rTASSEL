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
    tasselVersion <- "5.2.96"
    info          <- intToUtf8(0x2139)       # ℹ
    bold_on       <- "\033[1m"           # Start bold
    bold_off      <- "\033[22m"         # End bold

    msg <- sprintf(
        paste0(
            "Welcome to rTASSEL (version ", bold_on, "%s", bold_off, ")\n",
            "%s Running TASSEL version ", bold_on, "%s", bold_off, "\n",
            "%s Consider starting a TASSEL log file (see ?startLogger())\n"
        ),
        utils::packageVersion("rTASSEL"),
        info,
        tasselVersion,
        info
    )
    packageStartupMessage(msg)
}
