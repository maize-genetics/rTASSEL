## ----
# Get the HTTP status code for a URL
#
# @description
# Issues a header-only request and returns the numeric status code.
# Any transport-level failure (no network, DNS failure, timeout) yields
# \code{NA_integer_} rather than an error.
#
# @param url The URL to probe.
# @param timeout Request timeout in seconds.
#
# @return An integer status code, or \code{NA_integer_} on failure.
#
# @keywords internal
httpStatus <- function(url, timeout = TASSEL_UPDATE$TIMEOUT_SECS) {
    tryCatch({
        headers <- curlGetHeaders(url, timeout = as.integer(timeout))
        status  <- attr(headers, "status")
        if (is.null(status)) NA_integer_ else as.integer(status)
    }, error = function(e) NA_integer_)
}


## ----
# Test whether a URL resolves to an existing resource
#
# @param url The URL to probe.
# @param timeout Request timeout in seconds.
#
# @return \code{TRUE} if the server returned HTTP 200.
#
# @keywords internal
urlExists <- function(url, timeout = TASSEL_UPDATE$TIMEOUT_SECS) {
    isTRUE(httpStatus(url, timeout) == 200L)
}


## ----
# Read the contents of a URL as a single character string
#
# @param url The URL to read.
# @param timeout Read timeout in seconds.
# @param headers Named character vector of request headers, or \code{NULL}.
#
# @return A character string containing the full response body.
#
# @keywords internal
fetchUrlText <- function(
    url,
    timeout = TASSEL_UPDATE$TIMEOUT_SECS,
    headers = NULL
) {
    oldTimeout <- getOption("timeout")
    options(timeout = as.integer(timeout))
    on.exit(options(timeout = oldTimeout), add = TRUE)

    con <- base::url(url, open = "rb", headers = headers)
    on.exit(close(con), add = TRUE)

    paste(readLines(con, warn = FALSE), collapse = "\n")
}


## ----
# Extract the values of all occurrences of a simple XML tag
#
# @description
# A deliberately small XML reader that avoids taking on an \code{xml2}
# dependency. It handles the flat, attribute-free element structure used
# by Maven metadata and POM files.
#
# @param txt Character string of XML content.
# @param tag Tag name to extract, without angle brackets.
#
# @return A character vector of tag values, possibly empty.
#
# @keywords internal
extractXmlTags <- function(txt, tag) {
    pattern <- sprintf("<%s>\\s*(.*?)\\s*</%s>", tag, tag)
    hits    <- regmatches(txt, gregexpr(pattern, txt, perl = TRUE))[[1]]

    if (!length(hits)) return(character(0))

    sub(pattern, "\\1", hits, perl = TRUE)
}


## ----
# Sort version strings from newest to oldest
#
# @description
# Non-numeric versions (snapshots, release candidates) are dropped, since
# rTASSEL only ever installs released numeric versions.
#
# @param versions A character vector of version strings.
#
# @return A character vector sorted in descending version order.
#
# @keywords internal
sortVersionsDesc <- function(versions) {
    versions <- versions[grepl("^[0-9]+(\\.[0-9]+)*$", versions)]

    if (!length(versions)) return(character(0))

    versions[order(numeric_version(versions), decreasing = TRUE)]
}


## ----
# Fetch and parse the TASSEL maven-metadata.xml file
#
# @param timeout Request timeout in seconds.
#
# @return A list with elements \code{release} (character or \code{NULL})
#   and \code{versions} (character vector, newest first).
#
# @keywords internal
fetchMavenMetadata <- function(timeout = TASSEL_UPDATE$TIMEOUT_SECS) {
    url <- sprintf(
        "%s/%s/%s/maven-metadata.xml",
        TASSEL_MAVEN$BASE_URL,
        TASSEL_MAVEN$GROUP_PATH,
        TASSEL_MAVEN$ARTIFACT_ID
    )

    txt      <- fetchUrlText(url, timeout)
    release  <- extractXmlTags(txt, "release")
    versions <- sortVersionsDesc(extractXmlTags(txt, "version"))

    list(
        release  = if (length(release)) release[[1]] else NULL,
        versions = versions
    )
}


## ----
# Determine how a given TASSEL version is packaged on Maven Central
#
# @description
# TASSEL releases through 5.2.96 shipped a \code{jar-with-dependencies}
# fat JAR. Later releases switched to a Gradle publish that emits a thin
# JAR plus a POM listing runtime dependencies, which must be resolved
# separately. Some published versions have neither usable form.
#
# @param version TASSEL version string.
# @param timeout Request timeout in seconds.
#
# @return One of \code{"fat"}, \code{"thin"}, or \code{"none"}.
#
# @keywords internal
probeArtifactLayout <- function(version, timeout = TASSEL_UPDATE$TIMEOUT_SECS) {
    baseUrl <- sprintf(
        "%s/%s/%s/%s",
        TASSEL_MAVEN$BASE_URL,
        TASSEL_MAVEN$GROUP_PATH,
        TASSEL_MAVEN$ARTIFACT_ID,
        version
    )

    fatJar <- sprintf("%s/%s", baseUrl, getTASSELJarName(version))
    if (urlExists(fatJar, timeout)) return("fat")

    thinJar <- sprintf("%s/%s-%s.jar", baseUrl, TASSEL_MAVEN$ARTIFACT_ID, version)
    pomFile <- sprintf("%s/%s-%s.pom", baseUrl, TASSEL_MAVEN$ARTIFACT_ID, version)
    if (urlExists(thinJar, timeout) && urlExists(pomFile, timeout)) return("thin")

    "none"
}


## ----
# Find the newest TASSEL version that rTASSEL can actually install
#
# @description
# Walks published versions from newest to oldest and returns the first
# one with a usable artifact layout. Probing is bounded by
# \code{maxProbe} so that a run of broken publishes cannot stall package
# startup.
#
# @param timeout Request timeout in seconds.
# @param maxProbe Maximum number of versions to probe.
#
# @return A list with elements \code{latest} and \code{layout}, or
#   \code{NULL} if no installable version was found.
#
# @keywords internal
latestInstallableTASSEL <- function(
    timeout  = TASSEL_UPDATE$TIMEOUT_SECS,
    maxProbe = 5L
) {
    meta <- fetchMavenMetadata(timeout)

    if (!length(meta$versions)) return(NULL)

    for (version in utils::head(meta$versions, maxProbe)) {
        layout <- probeArtifactLayout(version, timeout)

        if (!identical(layout, "none")) {
            return(list(latest = version, layout = layout))
        }
    }

    NULL
}


## ----
# Compare two version strings
#
# @description
# Returns \code{FALSE} rather than erroring when either value cannot be
# parsed, so that a hand-set version never breaks package startup.
#
# @param candidate Version string to test.
# @param current Version string to compare against.
#
# @return \code{TRUE} if \code{candidate} is strictly newer.
#
# @keywords internal
isNewerVersion <- function(candidate, current) {
    tryCatch(
        isTRUE(numeric_version(candidate) > numeric_version(current)),
        error = function(e) FALSE
    )
}


## ----
# Path to the on-disk update check cache
#
# @param file Cache filename, which differs per release channel.
#
# @return A character string path to the cache file.
#
# @keywords internal
updateCheckCachePath <- function(file = TASSEL_UPDATE$CACHE_FILE) {
    file.path(tools::R_user_dir("rTASSEL", "cache"), file)
}


## ----
# Read the cached update check result
#
# @param file Cache filename, which differs per release channel.
#
# @return The cached list, or \code{NULL} if absent or unreadable.
#
# @keywords internal
readUpdateCheckCache <- function(file = TASSEL_UPDATE$CACHE_FILE) {
    cacheFile <- updateCheckCachePath(file)

    if (!file.exists(cacheFile)) return(NULL)

    tryCatch(readRDS(cacheFile), error = function(e) NULL)
}


## ----
# Write the update check result to the on-disk cache
#
# @param result A list to cache.
# @param file Cache filename, which differs per release channel.
#
# @return \code{result}, invisibly.
#
# @keywords internal
writeUpdateCheckCache <- function(result, file = TASSEL_UPDATE$CACHE_FILE) {
    cacheFile <- updateCheckCachePath(file)

    tryCatch({
        dir.create(dirname(cacheFile), recursive = TRUE, showWarnings = FALSE)
        saveRDS(result, cacheFile)
    }, error = function(e) NULL)

    invisible(result)
}


## ----
# Decide whether an automatic update check should run
#
# @description
# The check is suppressed in every non-interactive context so that
# \code{R CMD check}, continuous integration, and scripted runs never
# reach the network.
#
# @return \code{TRUE} if the automatic check may run.
#
# @keywords internal
updateCheckEnabled <- function() {
    if (!isTRUE(getOption("rTASSEL.check.updates", default = TRUE))) return(FALSE)
    if (nzchar(Sys.getenv("RTASSEL_NO_VERSION_CHECK")))                return(FALSE)
    if (nzchar(Sys.getenv("CI")))                                      return(FALSE)
    if (any(grepl("^_R_CHECK_", names(Sys.getenv()))))                 return(FALSE)
    if (!interactive())                                                return(FALSE)

    TRUE
}


## ----
#' @title Check for a newer version of TASSEL
#'
#' @description
#' Queries Maven Central for the newest release of TASSEL that rTASSEL is
#' able to install, and reports it alongside the version currently in use.
#'
#' With \code{channel = "nightly"}, the TASSEL GitHub releases page is
#' queried for the newest nightly build instead.
#'
#' This check also runs automatically when the package is attached, at
#' most once per day, against whichever channel the installed version came
#' from. The result is cached on disk so that repeated calls to
#' \code{library(rTASSEL)} do not repeatedly hit the network.
#'
#' @param force
#' If \code{TRUE}, ignore both the cached result and the settings that
#' normally suppress automatic checks, and query the release channel
#' directly. Defaults to \code{FALSE}.
#' @param channel
#' Which release channel to query: \code{"maven"} for Maven Central
#' releases (the default) or \code{"nightly"} for GitHub nightly builds.
#'
#' @details
#' The automatic check is skipped in non-interactive sessions, during
#' \code{R CMD check}, and on continuous integration. It can be disabled
#' outright with either of:
#'
#' \preformatted{
#' options(rTASSEL.check.updates = FALSE)
#' Sys.setenv(RTASSEL_NO_VERSION_CHECK = "true")
#' }
#'
#' Versions are probed newest-first and only reported when the artifacts
#' rTASSEL needs are actually present, since a version can appear in
#' Maven metadata without being installable.
#'
#' Network failures are never propagated: if the release channel cannot be
#' reached, the last cached result is returned, or \code{NULL} if there
#' is none.
#'
#' @return
#' A list with element \code{latest} (the newest installable version) and
#' \code{checkedAt} (the time of the check), or \code{NULL} if no result is
#' available. Maven results also carry \code{layout} (\code{"fat"} or
#' \code{"thin"}); nightly results carry \code{tag}, \code{publishedAt},
#' and \code{channel}.
#'
#' @examples
#' \dontrun{
#' ## Check for updates, bypassing the daily cache
#' checkForTASSELUpdate(force = TRUE)
#'
#' ## Check for a newer nightly build
#' checkForTASSELUpdate(channel = "nightly")
#' }
#'
#' @export
checkForTASSELUpdate <- function(force = FALSE, channel = c("maven", "nightly")) {
    channel <- match.arg(channel)

    if (!force && !updateCheckEnabled()) return(NULL)

    cacheFile <- if (identical(channel, "nightly")) {
        TASSEL_GITHUB$CACHE_FILE
    } else {
        TASSEL_UPDATE$CACHE_FILE
    }

    cached <- readUpdateCheckCache(cacheFile)

    if (!force && !is.null(cached)) {
        age <- as.numeric(
            difftime(Sys.time(), cached$checkedAt, units = "secs")
        )

        if (!is.na(age) && age >= 0 && age < TASSEL_UPDATE$MAX_AGE_SECS) {
            return(cached)
        }
    }

    result <- if (identical(channel, "nightly")) {
        latestNightlyTASSEL()
    } else {
        tryCatch(latestInstallableTASSEL(), error = function(e) NULL)
    }

    # Fall back to a stale result rather than nothing when offline
    if (is.null(result)) return(cached)

    result$checkedAt <- Sys.time()
    writeUpdateCheckCache(result, cacheFile)

    result
}
