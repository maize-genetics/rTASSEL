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
            activeVersion <- getActiveTASSELVersion()

            # JARs from a user-supplied path need not match anything rTASSEL
            # recorded, so the recorded version says nothing about them.
            onNightly <- !identical(resolved$source, "option") &&
                isNightlyVersion(activeVersion)

            # A nightly reports only its release version to the JVM, so the
            # recorded version is the only place its build date survives.
            tasselVersion <- if (onNightly) {
                activeVersion
            } else {
                getLoadedTASSELVersion()
            }

            cli::cli_bullets(c(
                "i" = "Running TASSEL version {.val {tasselVersion}} ({.field {resolved$source}})",
                "i" = "Consider starting a TASSEL log file (see {.help [startLogger()](rTASSEL::startLogger)})"
            ))

            if (onNightly) {
                update <- checkForTASSELUpdate(channel = "nightly")

                if (!is.null(update) && isNewerNightly(update$latest, tasselVersion)) {
                    cli::cli_bullets(c(
                        "!" = "A newer nightly build {.val {update$latest}} is available on GitHub",
                        "i" = "Run {.run rTASSEL::setupTASSEL(source = \"github\")} to update"
                    ))
                }
            } else {
                update <- checkForTASSELUpdate()

                if (!is.null(update) && isNewerVersion(update$latest, tasselVersion)) {
                    cli::cli_bullets(c(
                        "!" = "TASSEL {.val {update$latest}} is available on Maven Central",
                        "i" = "Run {.run rTASSEL::setupTASSEL(version = \"{update$latest}\")} to update"
                    ))
                }
            }
        }
    })

    packageStartupMessage(paste(msg, collapse = "\n"))
}
