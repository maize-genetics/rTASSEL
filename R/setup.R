## ----
#' @title Get the TASSEL JAR filename for a given version
#'
#' @param version TASSEL version string.
#'
#' @return A character string JAR filename.
#'
#' @keywords internal
getTASSELJarName <- function(version = TASSEL_MAVEN$VERSION) {
    sprintf(
        "%s-%s-%s.jar",
        TASSEL_MAVEN$ARTIFACT_ID,
        version,
        TASSEL_MAVEN$CLASSIFIER
    )
}


## ----
#' @title Get the JAR filename for a thin (dependency-less) TASSEL release
#'
#' @param version TASSEL version string.
#'
#' @return A character string JAR filename.
#'
#' @keywords internal
getTASSELThinJarName <- function(version) {
    sprintf("%s-%s.jar", TASSEL_MAVEN$ARTIFACT_ID, version)
}


## ----
#' @title Build the Maven Central URL for a TASSEL artifact
#'
#' @param version TASSEL version string.
#' @param filename Artifact filename.
#'
#' @return A character string URL.
#'
#' @keywords internal
tasselArtifactUrl <- function(version, filename) {
    sprintf(
        "%s/%s/%s/%s/%s",
        TASSEL_MAVEN$BASE_URL,
        TASSEL_MAVEN$GROUP_PATH,
        TASSEL_MAVEN$ARTIFACT_ID,
        version,
        filename
    )
}


## ----
#' @title Get the root directory holding all cached TASSEL versions
#'
#' @return A character string path.
#'
#' @keywords internal
getTASSELJavaDir <- function() {
    file.path(tools::R_user_dir("rTASSEL", "cache"), "java")
}


## ----
#' @title Get the local cache directory for TASSEL JARs
#'
#' @param version TASSEL version string.
#'
#' @return A character string path to the versioned cache directory.
#'
#' @keywords internal
getTASSELCacheDir <- function(version = getActiveTASSELVersion()) {
    file.path(getTASSELJavaDir(), version)
}


## ----
#' @title Path to the file recording the active TASSEL version
#'
#' @return A character string path.
#'
#' @keywords internal
activeVersionFile <- function() {
    file.path(getTASSELJavaDir(), "active")
}


## ----
#' @title Get the TASSEL version rTASSEL should load
#'
#' @description
#' Resolves the active version in priority order:
#' \enumerate{
#'   \item \code{options(rTASSEL.tassel.version = ...)}
#'   \item The version recorded by the most recent \code{\link{setupTASSEL}}
#'   \item The version pinned by this release of rTASSEL
#' }
#'
#' @return A character string version.
#'
#' @keywords internal
getActiveTASSELVersion <- function() {
    fromOption <- getOption("rTASSEL.tassel.version", default = NULL)
    if (!is.null(fromOption) && nzchar(fromOption)) return(as.character(fromOption))

    versionFile <- activeVersionFile()

    if (file.exists(versionFile)) {
        recorded <- tryCatch(
            trimws(readLines(versionFile, warn = FALSE))[1],
            error = function(e) NA_character_
        )

        if (!is.na(recorded) && nzchar(recorded) &&
            dir.exists(file.path(getTASSELJavaDir(), recorded))) {
            return(recorded)
        }
    }

    TASSEL_MAVEN$VERSION
}


## ----
#' @title Get the TASSEL version actually loaded in the JVM
#'
#' @description
#' Asks TASSEL which version it is, which is the only reliable answer when
#' JARs come from a user-supplied path that need not match any version
#' rTASSEL has recorded. Falls back to the configured version if TASSEL
#' cannot be queried.
#'
#' @param fallback Version to return if the JVM cannot be queried.
#'
#' @return A character string version.
#'
#' @keywords internal
getLoadedTASSELVersion <- function(fallback = getActiveTASSELVersion()) {
    version <- tryCatch(
        rJava::J(TASSEL_JVM$VERSIONS)$INSTANCE$tasselVersion(),
        error = function(e) NULL
    )

    if (!is.character(version) || length(version) != 1L || !nzchar(version)) {
        return(fallback)
    }

    version
}


## ----
#' @title Record the TASSEL version that rTASSEL should load
#'
#' @param version TASSEL version string.
#'
#' @return \code{version}, invisibly.
#'
#' @keywords internal
setActiveTASSELVersion <- function(version) {
    versionFile <- activeVersionFile()

    tryCatch({
        dir.create(dirname(versionFile), recursive = TRUE, showWarnings = FALSE)
        writeLines(as.character(version), versionFile)
    }, error = function(e) NULL)

    invisible(version)
}


## ----
#' @title Get path to cached TASSEL JAR files
#'
#' @description
#' Returns the file path to the cached TASSEL JAR directory, or \code{NULL}
#' if no cached JARs are found. JARs are cached per-version in the
#' standard R user cache directory.
#'
#' Both packaging layouts are recognised: a single fat JAR, and a thin JAR
#' accompanied by separately resolved dependency JARs.
#'
#' @param version
#' TASSEL version string. Defaults to the currently active version.
#'
#' @return A character string path to the JAR cache directory, or
#'   \code{NULL} if no cached JARs exist.
#'
#' @keywords internal
getTASSELJarPath <- function(version = getActiveTASSELVersion()) {
    jarDir <- getTASSELCacheDir(version)

    if (!dir.exists(jarDir)) return(NULL)

    candidates <- c(
        getTASSELJarName(version),
        getTASSELThinJarName(version)
    )

    if (any(file.exists(file.path(jarDir, candidates)))) return(jarDir)

    NULL
}


## ----
#' @title Resolve the JAR path from all available sources
#'
#' @description
#' Resolves the JAR path in priority order:
#' \enumerate{
#'   \item User-defined path via \code{options(rTASSEL.java.path = ...)}
#'   \item Maven cache (from \code{\link{setupTASSEL}})
#'   \item Bundled \code{inst/java/} (legacy fallback)
#' }
#'
#' @param pkgname Package name (used for bundled path lookup).
#' @param libname Library location (used for bundled path lookup).
#'
#' @return A list with elements \code{path} (character or \code{NULL}) and
#'   \code{source} (\code{"option"}, \code{"maven cache"},
#'   \code{"bundled"}, or \code{NULL}).
#'
#' @keywords internal
resolveJarPath <- function(pkgname = "rTASSEL", libname = NULL) {
    # 1. User-defined path via option
    jarPath <- getOption("rTASSEL.java.path", default = NULL)
    if (!is.null(jarPath)) {
        return(list(path = jarPath, source = "option"))
    }

    # 2. Maven cache (from setupTASSEL())
    jarPath <- getTASSELJarPath(getActiveTASSELVersion())
    if (!is.null(jarPath)) {
        return(list(path = jarPath, source = "maven cache"))
    }

    # 3. Bundled inst/java/ (legacy fallback)
    bundledPath <- system.file("java", package = pkgname, lib.loc = libname)
    if (dir.exists(bundledPath) && length(list.files(bundledPath, "\\.jar$")) > 0) {
        return(list(path = bundledPath, source = "bundled"))
    }

    list(path = NULL, source = NULL)
}


## ----
# Download a file with a styled CLI progress bar
#
# @description
# Downloads a file from a URL using chunked binary reads and displays
# a \code{cli} progress bar with spinner, visual bar, percentage,
# human-readable sizes, and elapsed time.
#
# @param url
# The URL to download from.
# @param destfile
# Local path to write the file to.
# @param estimatedSizeMb
# Fallback estimate (MB) used when the server does not report
# \code{Content-Length}.
#
# @return The \code{destfile} path (invisibly).
#
# @keywords internal
downloadWithProgress <- function(url, destfile, estimatedSizeMb = 70) {
    # -- resolve total size from HTTP headers (fall back to estimate) --
    totalBytes <- as.numeric(estimatedSizeMb) * 1024^2
    tryCatch({
        headers <- curlGetHeaders(url)
        clLines <- grep("^content-length:", headers,
                        ignore.case = TRUE, value = TRUE)
        if (length(clLines) > 0) {
            parsed <- as.numeric(
                trimws(sub("^content-length:\\s*", "",
                           clLines[length(clLines)],
                           ignore.case = TRUE))
            )
            if (!is.na(parsed) && parsed > 0) totalBytes <- parsed
        }
    }, error = function(e) NULL)

    # -- helper: human-readable byte sizes --
    fmtBytes <- function(b) {
        if (b >= 1024^3)      sprintf("%.1f GB", b / 1024^3)
        else if (b >= 1024^2) sprintf("%.1f MB", b / 1024^2)
        else if (b >= 1024)   sprintf("%.0f KB", b / 1024)
        else                  sprintf("%d B",    as.integer(b))
    }

    # -- open connections --
    conIn  <- url(url, "rb")
    on.exit(close(conIn), add = TRUE)
    conOut <- file(destfile, "wb")
    on.exit(close(conOut), add = TRUE)

    # -- chunked download with CLI progress bar --
    chunkSize  <- 64L * 1024L
    downloaded <- 0

    cli::cli_progress_bar(
        format = paste0(
            "{cli::pb_spin} Downloading ",
            "{cli::pb_bar} {cli::pb_percent} ",
            "| {fmtBytes(cli::pb_current)}/{fmtBytes(cli::pb_total)} ",
            "[{cli::pb_elapsed}]"
        ),
        format_done = paste0(
            "{cli::col_green(cli::symbol$tick)} Downloaded ",
            "{.strong {fmtBytes(cli::pb_current)}} ",
            "in {cli::pb_elapsed}"
        ),
        total = totalBytes,
        clear = FALSE
    )

    repeat {
        chunk <- readBin(conIn, "raw", n = chunkSize)
        if (length(chunk) == 0L) break
        writeBin(chunk, conOut)
        downloaded <- downloaded + length(chunk)
        cli::cli_progress_update(set = downloaded)
    }

    cli::cli_progress_done()
    invisible(destfile)
}


## ----
# Verify a downloaded TASSEL artifact against its published checksum
#
# @description
# The SHA-1 sidecar published alongside the artifact is authoritative, so
# any version can be verified rather than only the pinned default. For the
# pinned version the sidecar is additionally cross-checked against the
# checksum recorded in \code{TASSEL_MAVEN}, which would catch an artifact
# being replaced upstream.
#
# @param jarFile Path to the downloaded file.
# @param url URL the file was downloaded from.
# @param version TASSEL version string.
#
# @return \code{TRUE} if verification passed or was skipped.
#
# @keywords internal
verifyTASSELChecksum <- function(jarFile, url, version) {
    expected <- fetchArtifactSha1(url)
    isPinned <- identical(as.character(version), TASSEL_MAVEN$VERSION)

    if (isPinned) {
        if (is.null(expected)) {
            expected <- TASSEL_MAVEN$SHA1_CHECKSUM
        } else if (!identical(expected, TASSEL_MAVEN$SHA1_CHECKSUM)) {
            cli::cli_warn(c(
                "!" = "Published checksum for TASSEL {version} does not match the one recorded in rTASSEL",
                "i" = "Published: {expected}",
                "i" = "Recorded: {TASSEL_MAVEN$SHA1_CHECKSUM}"
            ))
        }
    }

    if (is.null(expected)) {
        cli::cli_alert_warning("No published checksum found; skipping verification")
        return(invisible(TRUE))
    }

    cli::cli_alert_info("Verifying SHA-1 checksum...")
    actualSha1 <- digest::digest(jarFile, algo = "sha1", file = TRUE)

    if (!identical(actualSha1, expected)) {
        file.remove(jarFile)
        cli::cli_abort(c(
            "SHA-1 checksum verification failed",
            "x" = "Expected: {expected}",
            "x" = "Got: {actualSha1}",
            "i" = "The download may be corrupted. Try again with {.code setupTASSEL(force = TRUE)}"
        ))
    }

    cli::cli_alert_success("Checksum verified")

    invisible(TRUE)
}


## ----
# Download a single TASSEL JAR from Maven Central
#
# @param version TASSEL version string.
# @param jarName Filename of the artifact to download.
# @param jarDir Destination directory.
#
# @return The downloaded file path, invisibly.
#
# @keywords internal
downloadTASSELJar <- function(version, jarName, jarDir) {
    jarFile <- file.path(jarDir, jarName)
    url     <- tasselArtifactUrl(version, jarName)

    cli::cli_alert_info("Downloading TASSEL {version} from Maven Central...")
    cli::cli_alert("URL: {.url {url}}")

    tryCatch({
        downloadWithProgress(url, jarFile)
    }, error = function(e) {
        # Clean up partial download
        if (file.exists(jarFile)) file.remove(jarFile)
        cli::cli_abort(c(
            "Failed to download TASSEL from Maven Central",
            "x" = conditionMessage(e),
            "i" = "Check your internet connection and try again"
        ))
    })

    verifyTASSELChecksum(jarFile, url, version)

    invisible(jarFile)
}


## ----
#' @title Download and configure TASSEL JAR files from Maven Central
#'
#' @description
#' Downloads TASSEL from Maven Central and caches it locally. This only
#' needs to be run once per TASSEL version; subsequent package loads use
#' the cached JARs automatically.
#'
#' @param version
#' TASSEL version to download. Defaults to the version pinned by this
#' release of rTASSEL.
#' @param force
#' If \code{TRUE}, discard any cached copy and re-download.
#' Defaults to \code{FALSE}.
#'
#' @details
#' TASSEL is published in one of two layouts, and the appropriate one is
#' detected automatically:
#'
#' \describe{
#'   \item{Fat JAR}{Releases up to 5.2.96 ship a single
#'     \code{jar-with-dependencies} archive containing every dependency.}
#'   \item{Thin JAR}{Later releases ship only TASSEL's own classes, so its
#'     dependencies are resolved from the POM and downloaded alongside it.
#'     This pulls artifacts from Maven Central, SciJava, and JBoss, since
#'     no single repository serves the whole graph.}
#' }
#'
#' Files are cached under the standard R user cache directory
#' (see \code{\link[tools]{R_user_dir}}) at
#' \code{~/.cache/R/rTASSEL/java/<version>/} (Linux),
#' \code{~/Library/Caches/org.R-project.R/R/rTASSEL/java/<version>/} (macOS),
#' or the equivalent on Windows.
#'
#' Every downloaded artifact is verified against the SHA-1 checksum
#' published next to it.
#'
#' On success the version is recorded as active, so the next
#' \code{library(rTASSEL)} loads it.
#'
#' @return The path to the JAR cache directory (invisibly).
#'
#' @examples
#' \dontrun{
#' ## Download the default TASSEL version
#' setupTASSEL()
#'
#' ## Install a specific version
#' setupTASSEL(version = "5.2.96")
#'
#' ## Force re-download
#' setupTASSEL(force = TRUE)
#' }
#'
#' @export
setupTASSEL <- function(version = TASSEL_MAVEN$VERSION, force = FALSE) {
    version <- as.character(version)
    jarDir  <- getTASSELCacheDir(version)

    if (!is.null(getTASSELJarPath(version)) && !force) {
        cli::cli_alert_info("TASSEL {version} JARs already cached at {.path {jarDir}}")
        cli::cli_alert_info("Use {.code setupTASSEL(force = TRUE)} to re-download")
        setActiveTASSELVersion(version)
        return(invisible(jarDir))
    }

    # A forced re-install must not leave dependency JARs from a previous
    # attempt behind, since they would stay on the classpath.
    if (force && dir.exists(jarDir)) {
        unlink(jarDir, recursive = TRUE)
    }

    dir.create(jarDir, recursive = TRUE, showWarnings = FALSE)

    layout <- probeArtifactLayout(version)

    if (identical(layout, "none")) {
        newest <- tryCatch(latestInstallableTASSEL(), error = function(e) NULL)

        cli::cli_abort(c(
            "TASSEL {version} is not available in a form rTASSEL can install",
            "x" = "Maven Central has neither a {TASSEL_MAVEN$CLASSIFIER} JAR nor a JAR and POM pair for this version",
            "i" = if (is.null(newest)) {
                "Check {.url https://central.sonatype.com/artifact/net.maizegenetics/tassel} for available versions"
            } else {
                "The newest installable version is {.val {newest$latest}}"
            }
        ))
    }

    if (identical(layout, "fat")) {
        downloadTASSELJar(version, getTASSELJarName(version), jarDir)
    } else {
        downloadTASSELJar(version, getTASSELThinJarName(version), jarDir)

        cli::cli_alert_info(
            "TASSEL {version} is published without bundled dependencies"
        )

        deps <- resolveDependencies(
            TASSEL_MAVEN$GROUP_ID,
            TASSEL_MAVEN$ARTIFACT_ID,
            version
        )
        downloadArtifacts(deps, jarDir)
    }

    setActiveTASSELVersion(version)

    cli::cli_alert_success("TASSEL {version} cached at {.path {jarDir}}")
    cli::cli_alert_info("Restart R or call {.code library(rTASSEL)} to load the new JARs")

    invisible(jarDir)
}
