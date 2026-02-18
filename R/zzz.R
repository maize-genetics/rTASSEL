.onLoad <- function(libname, pkgname) {
    resolved <- resolveJarPath(pkgname, libname)

    if (is.null(resolved$path)) {
        return(invisible())
    }

    # Initialize JVM and add JARs to classpath
    rJava::.jpackage(pkgname, lib.loc = libname)
    jars <- list.files(resolved$path, pattern = "\\.jar$", full.names = TRUE)
    rJava::.jaddClassPath(jars)

    # Start up a temp logging file regardless of user input
    # NOTE: users can define a logging file whenever they want - this is to
    #       prevent initial TASSEL stdout from contaminating the console space
    #       if a user has not defined one on initial use
    startLogger(tempfile(), verbose = FALSE)
}

.onAttach <- function(libname, pkgname){
    pkgVersion <- utils::packageVersion("rTASSEL")
    resolved   <- resolveJarPath(pkgname, libname)

    msg <- cli::cli_fmt({
        cli::cli_div(theme = list(h2 = list("margin-top" = 0, "margin-bottom" = 0)))
        cli::cli_h2("Welcome to rTASSEL (version {.val {pkgVersion}})")

        if (is.null(resolved$path)) {
            cli::cli_bullets(c(
                "i" = "TASSEL JARs not found",
                "i" = "Run {.run rTASSEL::setupTASSEL()} to download from Maven Central"
            ))
        } else {
            cli::cli_bullets(c(
                "i" = "Running TASSEL version {.val {TASSEL_MAVEN$VERSION}} ({.field {resolved$source}})",
                "i" = "Consider starting a TASSEL log file (see {.help [startLogger()](rTASSEL::startLogger)})"
            ))
        }
    })

    packageStartupMessage(paste(msg, collapse = "\n"))
}
