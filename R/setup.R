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
#' @title Get the local cache directory for TASSEL JARs
#'
#' @param version TASSEL version string.
#'
#' @return A character string path to the versioned cache directory.
#'
#' @keywords internal
getTASSELCacheDir <- function(version = TASSEL_MAVEN$VERSION) {
    file.path(tools::R_user_dir("rTASSEL", "cache"), "java", version)
}


## ----
#' @title Get path to cached TASSEL JAR files
#'
#' @description
#' Returns the file path to the cached TASSEL JAR directory, or \code{NULL}
#' if no cached JARs are found. JARs are cached per-version in the
#' standard R user cache directory.
#'
#' @param version
#' TASSEL version string. Defaults to the version bundled with this
#' release of rTASSEL.
#'
#' @return A character string path to the JAR cache directory, or
#'   \code{NULL} if no cached JARs exist.
#'
#' @keywords internal
getTASSELJarPath <- function(version = TASSEL_MAVEN$VERSION) {
    jarDir  <- getTASSELCacheDir(version)
    jarFile <- file.path(jarDir, getTASSELJarName(version))

    if (file.exists(jarFile)) return(jarDir)

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
    jarPath <- getTASSELJarPath()
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
#' @title Download and configure TASSEL JAR files from Maven Central
#'
#' @description
#' Downloads the TASSEL fat JAR (with all dependencies) from Maven
#' Central and caches it locally. This only needs to be run once per
#' TASSEL version. Subsequent package loads will use the cached JAR
#' automatically.
#'
#' @param version
#' TASSEL version to download. Defaults to \code{"5.2.96"}.
#' @param force
#' If \code{TRUE}, re-download even if a cached version exists.
#' Defaults to \code{FALSE}.
#'
#' @details
#' The JAR is downloaded from Maven Central at:
#' \url{https://mvnrepository.com/artifact/net.maizegenetics/tassel}
#'
#' Files are cached under the standard R user cache directory
#' (see \code{\link[tools]{R_user_dir}}) at
#' \code{~/.cache/R/rTASSEL/java/<version>/} (Linux),
#' \code{~/Library/Caches/org.R-project.R/R/rTASSEL/java/<version>/} (macOS),
#' or the equivalent on Windows.
#'
#' A SHA-1 checksum is verified after download to ensure file integrity.
#'
#' @return The path to the JAR cache directory (invisibly).
#'
#' @examples
#' \dontrun{
#' ## Download default TASSEL version
#' setupTASSEL()
#'
#' ## Force re-download
#' setupTASSEL(force = TRUE)
#' }
#'
#' @export
setupTASSEL <- function(version = TASSEL_MAVEN$VERSION, force = FALSE) {
    jarDir  <- getTASSELCacheDir(version)
    jarName <- getTASSELJarName(version)
    jarFile <- file.path(jarDir, jarName)

    if (file.exists(jarFile) && !force) {
        cli::cli_alert_info("TASSEL {version} JARs already cached at {.path {jarDir}}")
        cli::cli_alert_info("Use {.code setupTASSEL(force = TRUE)} to re-download")
        return(invisible(jarDir))
    }

    dir.create(jarDir, recursive = TRUE, showWarnings = FALSE)

    url <- sprintf(
        "%s/%s/%s/%s/%s",
        TASSEL_MAVEN$BASE_URL,
        TASSEL_MAVEN$GROUP_PATH,
        TASSEL_MAVEN$ARTIFACT_ID,
        version,
        jarName
    )

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

    # Verify SHA-1 checksum
    if (version == TASSEL_MAVEN$VERSION) {
        cli::cli_alert_info("Verifying SHA-1 checksum...")
        actualSha1 <- digest::digest(jarFile, algo = "sha1", file = TRUE)

        if (actualSha1 != TASSEL_MAVEN$SHA1_CHECKSUM) {
            file.remove(jarFile)
            cli::cli_abort(c(
                "SHA-1 checksum verification failed",
                "x" = "Expected: {TASSEL_MAVEN$SHA1_CHECKSUM}",
                "x" = "Got: {actualSha1}",
                "i" = "The download may be corrupted. Try again with {.code setupTASSEL(force = TRUE)}"
            ))
        }
        cli::cli_alert_success("Checksum verified")
    }

    cli::cli_alert_success("TASSEL {version} cached at {.path {jarDir}}")
    cli::cli_alert_info("Restart R or call {.code library(rTASSEL)} to load the new JARs")

    invisible(jarDir)
}
