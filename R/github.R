# Installing TASSEL from the standalone archives on GitHub.
#
# Nightly builds are published only as release assets on the TASSEL GitHub
# repository, never to Maven Central, so this is the only route to a
# bleeding-edge TASSEL. The same archives are published for tagged releases,
# which makes a single version spec cover both channels.
#
# Each archive unpacks to a flat tree holding 'sTASSEL.jar' next to a 'lib'
# directory of dependency JARs. Those are flattened into the per-version cache
# directory the Maven installer also writes to, so nothing about classpath
# construction has to know where an installation came from.


# Grammar of a nightly version: a release version with the build date of the
# 'develop' snapshot appended, as in "5.2.98-dev.20260801".
GITHUB_NIGHTLY_PATTERN <- "^([0-9]+(?:\\.[0-9]+)*)-dev\\.([0-9]{8})$"

# Specs that select the newest release of a channel rather than naming one.
GITHUB_NIGHTLY_ALIASES <- c("nightly", "dev")
GITHUB_RELEASE_ALIASES <- c("latest", "stable", "release")


## ----
# Build the headers for a GitHub API request
#
# @description
# A token is used whenever one is present in the environment. The anonymous
# rate limit is 60 requests per hour per address and is shared with every
# other tool on the same network, so it can already be exhausted on arrival.
#
# @return A named character vector of headers.
#
# @keywords internal
githubApiHeaders <- function() {
    headers <- c("Accept" = "application/vnd.github+json")

    token <- Sys.getenv(c("GITHUB_PAT", "GITHUB_TOKEN"))
    token <- token[nzchar(token)]

    if (length(token)) {
        headers <- c(headers, "Authorization" = paste("Bearer", token[[1]]))
    }

    headers
}


## ----
# Fetch and parse a GitHub API response
#
# @param path API path relative to the API root.
# @param timeout Request timeout in seconds.
#
# @return The parsed response, as nested lists.
#
# @keywords internal
fetchGitHubJson <- function(path, timeout = TASSEL_GITHUB$TIMEOUT_SECS) {
    url <- sprintf("%s/%s", TASSEL_GITHUB$API_BASE, sub("^/", "", path))
    txt <- fetchUrlText(url, timeout, headers = githubApiHeaders())

    jsonlite::fromJSON(txt, simplifyVector = FALSE)
}


## ----
# Recover the TASSEL version from a standalone archive filename
#
# @param assetName Asset filename.
#
# @return A character string version, or \code{NA_character_}.
#
# @keywords internal
standaloneAssetVersion <- function(assetName) {
    if (!startsWith(assetName, TASSEL_GITHUB$ASSET_PREFIX) ||
        !endsWith(assetName, TASSEL_GITHUB$ASSET_EXT)) {
        return(NA_character_)
    }

    version <- substr(
        assetName,
        nchar(TASSEL_GITHUB$ASSET_PREFIX) + 1L,
        nchar(assetName) - nchar(TASSEL_GITHUB$ASSET_EXT)
    )

    if (!nzchar(version)) NA_character_ else version
}


## ----
# Find the standalone archive among a release's assets
#
# @description
# Releases also carry platform installers and a zip of the same archive,
# neither of which rTASSEL can use, so assets are matched by name.
#
# @param release A parsed release entry.
#
# @return The matching asset entry, or \code{NULL}.
#
# @keywords internal
standaloneAsset <- function(release) {
    for (asset in release$assets) {
        name <- asset$name

        if (!is.character(name) || length(name) != 1L) next

        if (startsWith(name, TASSEL_GITHUB$ASSET_PREFIX) &&
            endsWith(name, TASSEL_GITHUB$ASSET_EXT)) {
            return(asset)
        }
    }

    NULL
}


## ----
# Reduce a release entry to the fields the installer needs
#
# @param release A parsed release entry.
#
# @return A list describing the release, or \code{NULL} if it carries no
#   installable archive.
#
# @keywords internal
asReleaseRecord <- function(release) {
    if (isTRUE(release$draft)) return(NULL)

    asset <- standaloneAsset(release)

    if (is.null(asset)) return(NULL)

    version <- standaloneAssetVersion(asset$name)

    if (is.na(version)) return(NULL)

    sha256 <- if (is.character(asset$digest) &&
                  grepl("^sha256:", asset$digest)) {
        sub("^sha256:", "", asset$digest)
    } else {
        NULL
    }

    list(
        tag         = as.character(release$tag_name),
        version     = version,
        publishedAt = as.character(release$published_at),
        prerelease  = isTRUE(release$prerelease),
        assetName   = as.character(asset$name),
        assetUrl    = as.character(asset$browser_download_url),
        sha256      = sha256,
        size        = if (is.null(asset$size)) NA_real_ else as.numeric(asset$size)
    )
}


## ----
# List the GitHub releases that ship an installable TASSEL archive
#
# @param timeout Request timeout in seconds.
#
# @return A list of release records, newest first.
#
# @keywords internal
githubReleaseCandidates <- function(timeout = TASSEL_GITHUB$TIMEOUT_SECS) {
    releases <- fetchGitHubJson(
        sprintf(
            "repos/%s/releases?per_page=%d",
            TASSEL_GITHUB$REPO,
            TASSEL_GITHUB$PER_PAGE
        ),
        timeout
    )

    records <- lapply(releases, asReleaseRecord)
    records <- records[!vapply(records, is.null, logical(1))]

    if (!length(records)) return(list())

    published <- vapply(records, function(x) x$publishedAt, character(1))

    records[order(published, decreasing = TRUE)]
}


## ----
# Resolve a version spec against the published GitHub releases
#
# @param spec
# One of \code{"nightly"} (newest prerelease), \code{"latest"} (newest
# tagged release), a release tag such as \code{"dev-20260801"}, or a
# version such as \code{"5.2.98-dev.20260801"}.
# @param timeout Request timeout in seconds.
#
# @return A release record.
#
# @keywords internal
resolveGitHubRelease <- function(spec = "nightly", timeout = TASSEL_GITHUB$TIMEOUT_SECS) {
    spec     <- as.character(spec)
    releases <- githubReleaseCandidates(timeout)

    if (!length(releases)) {
        cli::cli_abort(c(
            "No installable TASSEL archives found on GitHub",
            "i" = "Check {.url https://github.com/{TASSEL_GITHUB$REPO}/releases}"
        ))
    }

    isPrerelease <- vapply(releases, function(x) x$prerelease, logical(1))

    matched <- if (spec %in% GITHUB_NIGHTLY_ALIASES) {
        releases[isPrerelease]
    } else if (spec %in% GITHUB_RELEASE_ALIASES) {
        releases[!isPrerelease]
    } else {
        tags     <- vapply(releases, function(x) x$tag, character(1))
        versions <- vapply(releases, function(x) x$version, character(1))

        releases[tags == spec | versions == spec | paste0("v", versions) == spec]
    }

    if (!length(matched)) {
        available <- vapply(releases, function(x) x$version, character(1))

        cli::cli_abort(c(
            "No TASSEL archive on GitHub matches {.val {spec}}",
            "i" = "Available versions include {.val {utils::head(available, 5)}}",
            "i" = "Use {.val nightly} for the newest nightly build or {.val latest} for the newest tagged release"
        ))
    }

    matched[[1]]
}


## ----
# Find the newest nightly build published on GitHub
#
# @description
# Network and parsing failures yield \code{NULL} rather than an error, so a
# GitHub outage cannot break package startup.
#
# @param timeout Request timeout in seconds.
#
# @return A list with elements \code{latest}, \code{tag}, \code{publishedAt},
#   and \code{channel}, or \code{NULL}.
#
# @keywords internal
latestNightlyTASSEL <- function(timeout = TASSEL_GITHUB$TIMEOUT_SECS) {
    release <- tryCatch(
        resolveGitHubRelease("nightly", timeout),
        error = function(e) NULL
    )

    if (is.null(release)) return(NULL)

    list(
        latest      = release$version,
        tag         = release$tag,
        publishedAt = release$publishedAt,
        channel     = "nightly"
    )
}


## ----
# Test whether a version string names a nightly build
#
# @param version Version string.
#
# @return \code{TRUE} if the version follows the nightly grammar.
#
# @keywords internal
isNightlyVersion <- function(version) {
    if (!is.character(version) || length(version) != 1L || is.na(version)) {
        return(FALSE)
    }

    grepl(GITHUB_NIGHTLY_PATTERN, version, perl = TRUE)
}


## ----
# Split a nightly version into its release version and build date
#
# @param version Version string.
#
# @return A list with elements \code{base} and \code{date}, or \code{NULL}
#   if the version is not a nightly.
#
# @keywords internal
parseNightlyVersion <- function(version) {
    if (!isNightlyVersion(version)) return(NULL)

    list(
        base = sub(GITHUB_NIGHTLY_PATTERN, "\\1", version, perl = TRUE),
        date = sub(GITHUB_NIGHTLY_PATTERN, "\\2", version, perl = TRUE)
    )
}


## ----
# Compare two nightly versions
#
# @description
# \code{numeric_version} cannot parse the nightly suffix, so the release
# version and the build date are compared separately. Build dates are
# zero-padded and so compare correctly as strings.
#
# @param candidate Nightly version to test.
# @param current Nightly version to compare against.
#
# @return \code{TRUE} if \code{candidate} is a strictly newer nightly.
#
# @keywords internal
isNewerNightly <- function(candidate, current) {
    cand <- parseNightlyVersion(candidate)
    curr <- parseNightlyVersion(current)

    if (is.null(cand) || is.null(curr)) return(FALSE)

    if (!identical(cand$base, curr$base)) {
        return(isNewerVersion(cand$base, curr$base))
    }

    cand$date > curr$date
}


## ----
# Verify a downloaded archive against the checksum GitHub published for it
#
# @param archive Path to the downloaded archive.
# @param release The release record the archive was downloaded from.
#
# @return \code{TRUE} if verification passed or was skipped.
#
# @keywords internal
verifyGitHubChecksum <- function(archive, release) {
    if (is.null(release$sha256)) {
        cli::cli_alert_warning("No published checksum found; skipping verification")
        return(invisible(TRUE))
    }

    cli::cli_alert_info("Verifying SHA-256 checksum...")
    actual <- digest::digest(archive, algo = "sha256", file = TRUE)

    if (!identical(actual, release$sha256)) {
        file.remove(archive)
        cli::cli_abort(c(
            "SHA-256 checksum verification failed",
            "x" = "Expected: {release$sha256}",
            "x" = "Got: {actual}",
            "i" = "The download may be corrupted. Try again with {.code setupTASSEL(source = \"github\", force = TRUE)}"
        ))
    }

    cli::cli_alert_success("Checksum verified")

    invisible(TRUE)
}


## ----
# Locate the root of an extracted standalone archive
#
# @description
# Current archives unpack their contents directly into the extraction
# directory, but a single wrapping directory is also tolerated so that a
# change to how the archive is built does not break installation.
#
# @param exdir Directory the archive was extracted into.
#
# @return A character string path, or \code{NULL} if no root was found.
#
# @keywords internal
standaloneRoot <- function(exdir) {
    if (file.exists(file.path(exdir, TASSEL_GITHUB$MAIN_JAR))) return(exdir)

    for (candidate in list.dirs(exdir, recursive = FALSE)) {
        if (file.exists(file.path(candidate, TASSEL_GITHUB$MAIN_JAR))) {
            return(candidate)
        }
    }

    NULL
}


## ----
# Download and unpack a TASSEL standalone archive into the JAR cache
#
# @description
# The main JAR and every dependency JAR are copied into a single flat
# directory, matching the layout the Maven installer produces.
#
# @param release A release record from \code{\link{resolveGitHubRelease}}.
# @param jarDir Destination directory.
#
# @return A character vector of installed JAR paths, invisibly.
#
# @keywords internal
installGitHubStandalone <- function(release, jarDir) {
    # A directory of its own, so that neither the download nor the unpacked
    # tree can collide with anything else living in the session temp dir.
    workDir <- file.path(tempdir(), paste0("rtassel-standalone-", release$version))
    unlink(workDir, recursive = TRUE)
    dir.create(workDir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(workDir, recursive = TRUE), add = TRUE)

    archive <- file.path(workDir, release$assetName)
    exdir   <- file.path(workDir, "unpacked")

    cli::cli_alert_info("Downloading TASSEL {release$version} from GitHub...")
    cli::cli_alert("URL: {.url {release$assetUrl}}")

    sizeMb <- if (length(release$size) != 1L || is.na(release$size)) {
        70
    } else {
        release$size / 1024^2
    }

    tryCatch({
        downloadWithProgress(release$assetUrl, archive, estimatedSizeMb = sizeMb)
    }, error = function(e) {
        if (file.exists(archive)) file.remove(archive)
        cli::cli_abort(c(
            "Failed to download the TASSEL standalone archive",
            "x" = conditionMessage(e),
            "i" = "Check your internet connection and try again"
        ))
    })

    verifyGitHubChecksum(archive, release)

    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

    cli::cli_alert_info("Extracting standalone archive...")

    extracted <- tryCatch(
        utils::untar(archive, exdir = exdir),
        error = function(e) conditionMessage(e)
    )

    if (!isTRUE(extracted == 0L)) {
        cli::cli_abort(c(
            "Failed to extract the TASSEL standalone archive",
            "x" = "{.path {release$assetName}} could not be unpacked",
            "i" = "Try again with {.code setupTASSEL(source = \"github\", force = TRUE)}"
        ))
    }

    root <- standaloneRoot(exdir)

    if (is.null(root)) {
        cli::cli_abort(c(
            "The TASSEL standalone archive does not contain {.file {TASSEL_GITHUB$MAIN_JAR}}",
            "i" = "This release may not be usable by rTASSEL"
        ))
    }

    jars <- c(
        file.path(root, TASSEL_GITHUB$MAIN_JAR),
        list.files(
            file.path(root, TASSEL_GITHUB$LIB_DIR),
            pattern    = "\\.jar$",
            full.names = TRUE
        )
    )

    dir.create(jarDir, recursive = TRUE, showWarnings = FALSE)
    copied <- file.copy(jars, jarDir, overwrite = TRUE)

    if (!all(copied)) {
        cli::cli_abort(c(
            "Could not copy {sum(!copied)} JAR{?s} into the cache",
            "x" = "{basename(jars[!copied])}",
            "i" = "Check that {.path {jarDir}} is writable"
        ))
    }

    cli::cli_alert_success(
        "Installed {length(jars)} JAR{?s} ({.file {TASSEL_GITHUB$MAIN_JAR}} and {length(jars) - 1L} dependenc{?y/ies})"
    )

    invisible(file.path(jarDir, basename(jars)))
}
