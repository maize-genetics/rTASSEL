# A small Maven dependency resolver.
#
# TASSEL releases after 5.2.96 are published as a thin JAR plus a POM that
# lists runtime dependencies, rather than as a single fat JAR. Installing
# those releases means walking the dependency graph ourselves.
#
# This implements the subset of Maven's resolution semantics that TASSEL's
# graph actually exercises: parent POM inheritance, property interpolation,
# dependency management (including BOM imports), scope and optional
# filtering, exclusions, and nearest-wins version conflict resolution.
# Anything outside that subset degrades to a warning rather than an error.


# Seconds to allow for artifact downloads, which are far larger than the
# metadata requests governed by TASSEL_UPDATE$TIMEOUT_SECS.
MAVEN_FETCH_TIMEOUT <- 60L

# Scopes that never contribute to a runtime classpath.
MAVEN_SKIPPED_SCOPES <- c("test", "provided", "system")

# Artifacts that ship the same classes as a newer coordinate under a
# different name. Maven's nearest-wins rule cannot deduplicate these,
# because it keys on groupId:artifactId, so both would land on the
# classpath and whichever is read first would win. Gradle avoids this via
# capability conflict detection, which has no equivalent here.
#
# google-collections is the pre-1.0 incarnation of Guava and provides its
# own, older com.google.common.base.Preconditions. Leaving it alongside
# Guava causes NoSuchMethodError once TASSEL calls a newer overload.
MAVEN_SUPERSEDED <- c(
    "com.google.collections:google-collections" = "com.google.guava:guava"
)


## ----
# Build a "groupId:artifactId" conflict key
#
# @param groupId Maven group identifier.
# @param artifactId Maven artifact identifier.
#
# @return A character string key.
#
# @keywords internal
depKey <- function(groupId, artifactId) {
    paste0(groupId, ":", artifactId)
}


## ----
# Convert a dotted Maven group identifier to a repository path
#
# @param groupId Maven group identifier.
#
# @return A character string with dots replaced by slashes.
#
# @keywords internal
groupIdToPath <- function(groupId) {
    gsub(".", "/", groupId, fixed = TRUE)
}


## ----
# Build the repository-relative path for a Maven artifact
#
# @param groupId Maven group identifier.
# @param artifactId Maven artifact identifier.
# @param version Artifact version.
# @param ext File extension, without a leading dot.
# @param classifier Optional artifact classifier.
#
# @return A character string path.
#
# @keywords internal
mavenArtifactPath <- function(
    groupId,
    artifactId,
    version,
    ext        = "jar",
    classifier = NULL
) {
    suffix <- if (is.null(classifier) || !nzchar(classifier)) {
        ""
    } else {
        paste0("-", classifier)
    }

    sprintf(
        "%s/%s/%s/%s-%s%s.%s",
        groupIdToPath(groupId), artifactId, version,
        artifactId, version, suffix, ext
    )
}


## ----
# Locate an artifact across the configured Maven repositories
#
# @description
# Repositories are tried in the order given by
# \code{TASSEL_MAVEN$REPOS}; the first that serves the artifact wins.
#
# @param path Repository-relative artifact path.
# @param timeout Request timeout in seconds.
#
# @return The full URL of the artifact, or \code{NULL} if no repository
#   serves it.
#
# @keywords internal
findArtifactUrl <- function(path, timeout = TASSEL_UPDATE$TIMEOUT_SECS) {
    for (repo in TASSEL_MAVEN$REPOS) {
        url <- paste0(repo, "/", path)
        if (urlExists(url, timeout)) return(url)
    }

    NULL
}


## ----
# Extract the inner content of a non-nested XML element
#
# @param txt Character string of XML content.
# @param tag Tag name, without angle brackets.
#
# @return The element's inner content, or \code{""} when absent.
#
# @keywords internal
extractXmlBlock <- function(txt, tag) {
    pattern <- sprintf("(?s)<%s>(.*?)</%s>", tag, tag)
    match   <- regmatches(txt, regexpr(pattern, txt, perl = TRUE))

    if (!length(match)) return("")

    sub(pattern, "\\1", match, perl = TRUE)
}


## ----
# Extract every occurrence of a repeated XML element
#
# @param txt Character string of XML content.
# @param tag Tag name, without angle brackets.
#
# @return A character vector of inner contents, possibly empty.
#
# @keywords internal
extractXmlBlocks <- function(txt, tag) {
    pattern <- sprintf("(?s)<%s>(.*?)</%s>", tag, tag)
    hits    <- regmatches(txt, gregexpr(pattern, txt, perl = TRUE))[[1]]

    if (!length(hits)) return(character(0))

    sub(pattern, "\\1", hits, perl = TRUE)
}


## ----
# Read a single scalar XML value
#
# @param txt Character string of XML content.
# @param tag Tag name, without angle brackets.
# @param default Value to return when the tag is absent.
#
# @return The tag value, or \code{default}.
#
# @keywords internal
xmlValue <- function(txt, tag, default = NULL) {
    values <- extractXmlTags(txt, tag)

    if (!length(values)) return(default)

    values[[1]]
}


# Properties that describe the POM currently being parsed. They must never
# be inherited from a parent, or a coordinate such as ${project.version}
# would resolve to some distant ancestor's version.
MAVEN_BUILTIN_PROPS <- c(
    "project.version", "pom.version", "version",
    "project.groupId", "pom.groupId"
)


## ----
# Isolate the coordinate header of a POM
#
# @description
# Returns the portion of a POM that declares the project's own
# coordinates: the \code{<parent>} block is removed, and everything from
# the first container element onward is discarded. Without this, reading
# the project's own \code{<version>} would pick up the parent's version,
# or a version belonging to the first declared dependency.
#
# @param txt POM contents as a character string.
#
# @return The header portion of the POM.
#
# @keywords internal
pomHeader <- function(txt) {
    header <- sub("(?s)<parent>.*?</parent>", "", txt, perl = TRUE)

    containers <- paste0(
        "<(dependencyManagement|dependencies|properties|modules|profiles",
        "|build|reporting|distributionManagement)>"
    )
    cut <- regexpr(containers, header, perl = TRUE)

    if (cut > 0L) header <- substr(header, 1L, cut - 1L)

    header
}


## ----
# Substitute Maven property placeholders in a string
#
# @description
# Resolves \code{${...}} placeholders against the supplied property map.
# Properties may reference other properties, so substitution repeats
# until it reaches a fixed point or hits an iteration ceiling, which
# guards against circular definitions.
#
# @param x Character string possibly containing placeholders.
# @param properties Named character vector of property values.
#
# @return The interpolated string.
#
# @keywords internal
interpolateProperties <- function(x, properties) {
    if (is.null(x) || !length(x) || !nzchar(x)) return(x)

    for (i in seq_len(10L)) {
        if (!grepl("${", x, fixed = TRUE)) break

        placeholders <- regmatches(
            x,
            gregexpr("\\$\\{[^}]+\\}", x, perl = TRUE)
        )[[1]]

        if (!length(placeholders)) break

        before <- x
        for (placeholder in unique(placeholders)) {
            name <- sub("^\\$\\{(.*)\\}$", "\\1", placeholder)

            if (name %in% names(properties)) {
                x <- gsub(placeholder, properties[[name]], x, fixed = TRUE)
            }
        }

        # No placeholder could be resolved; further passes are pointless
        if (identical(before, x)) break
    }

    x
}


## ----
# Download and cache a POM file
#
# @description
# POMs are memoised for the lifetime of a resolution run, since a single
# dependency graph typically revisits the same parent and BOM POMs many
# times.
#
# @param groupId Maven group identifier.
# @param artifactId Maven artifact identifier.
# @param version Artifact version.
# @param cache An environment used to memoise results.
#
# @return The POM contents as a character string, or \code{NULL} if the
#   POM could not be retrieved from any repository.
#
# @keywords internal
fetchPom <- function(groupId, artifactId, version, cache = new.env(parent = emptyenv())) {
    cacheKey <- paste0(depKey(groupId, artifactId), ":", version)

    if (!is.null(cache[[cacheKey]])) {
        value <- cache[[cacheKey]]
        return(if (identical(value, "")) NULL else value)
    }

    path <- mavenArtifactPath(groupId, artifactId, version, ext = "pom")
    txt  <- NULL

    # Request each repository directly rather than probing first with a
    # header request: resolution walks hundreds of POMs, and halving the
    # round trips roughly halves the total resolution time.
    for (repo in TASSEL_MAVEN$REPOS) {
        # A miss in one repository is the normal path to the next, so the
        # connection warning it raises is not worth surfacing.
        txt <- suppressWarnings(tryCatch(
            fetchUrlText(paste0(repo, "/", path), MAVEN_FETCH_TIMEOUT),
            error = function(e) NULL
        ))
        if (!is.null(txt) && nzchar(txt)) break
        txt <- NULL
    }

    # Cache negative lookups as "" so a missing POM is only probed once
    cache[[cacheKey]] <- if (is.null(txt)) "" else txt

    txt
}


## ----
# Parse the exclusions declared on a dependency
#
# @param depTxt Inner XML of a \code{<dependency>} element.
#
# @return A character vector of "groupId:artifactId" exclusion keys.
#
# @keywords internal
parseExclusions <- function(depTxt) {
    exclusionsTxt <- extractXmlBlock(depTxt, "exclusions")

    if (!nzchar(exclusionsTxt)) return(character(0))

    blocks <- extractXmlBlocks(exclusionsTxt, "exclusion")

    vapply(
        blocks,
        function(block) {
            depKey(
                xmlValue(block, "groupId",    ""),
                xmlValue(block, "artifactId", "")
            )
        },
        character(1),
        USE.NAMES = FALSE
    )
}


## ----
# Parse the \code{<dependency>} elements out of a POM fragment
#
# @param txt Character string of XML containing dependency elements.
# @param properties Named character vector used for interpolation.
#
# @return A list of dependency records.
#
# @keywords internal
parseDependencyBlocks <- function(txt, properties) {
    blocks <- extractXmlBlocks(txt, "dependency")

    lapply(blocks, function(block) {
        # Exclusions are parsed separately; drop them so their nested
        # groupId and artifactId tags cannot be mistaken for the
        # dependency's own coordinates.
        exclusions <- parseExclusions(block)
        coords     <- gsub("(?s)<exclusions>.*?</exclusions>", "", block, perl = TRUE)

        list(
            groupId    = interpolateProperties(xmlValue(coords, "groupId", ""),    properties),
            artifactId = interpolateProperties(xmlValue(coords, "artifactId", ""), properties),
            version    = interpolateProperties(xmlValue(coords, "version", ""),    properties),
            scope      = xmlValue(coords, "scope", "compile"),
            type       = xmlValue(coords, "type", "jar"),
            classifier = xmlValue(coords, "classifier", ""),
            optional   = identical(tolower(xmlValue(coords, "optional", "false")), "true"),
            exclusions = exclusions
        )
    })
}


## ----
# Parse a POM into properties, managed versions, and dependencies
#
# @description
# Parent POMs are resolved recursively so that inherited properties and
# managed versions are available. Child declarations take precedence over
# inherited ones. BOM imports (\code{<scope>import</scope>}) contribute
# their own managed versions.
#
# @param txt POM contents as a character string.
# @param cache An environment used to memoise POM downloads.
# @param depth Current recursion depth, used to bound parent chains.
#
# @return A list with elements \code{properties}, \code{managed}, and
#   \code{dependencies}.
#
# @keywords internal
parsePom <- function(txt, cache = new.env(parent = emptyenv()), depth = 0L) {
    empty <- list(
        properties   = character(0),
        managed      = character(0),
        dependencies = list()
    )

    if (is.null(txt) || !nzchar(txt) || depth > 10L) return(empty)

    # Profiles and build plugins declare dependencies that are not part of
    # the runtime graph; remove them before anything else.
    txt <- gsub("(?s)<profiles>.*?</profiles>", "", txt, perl = TRUE)
    txt <- gsub("(?s)<build>.*?</build>",       "", txt, perl = TRUE)

    parentTxt <- extractXmlBlock(txt, "parent")
    parent    <- empty

    if (nzchar(parentTxt)) {
        parentPom <- fetchPom(
            xmlValue(parentTxt, "groupId",    ""),
            xmlValue(parentTxt, "artifactId", ""),
            xmlValue(parentTxt, "version",    ""),
            cache
        )
        parent <- parsePom(parentPom, cache, depth + 1L)
    }

    header     <- pomHeader(txt)
    ownVersion <- xmlValue(header, "version", xmlValue(parentTxt, "version", ""))
    ownGroup   <- xmlValue(header, "groupId", xmlValue(parentTxt, "groupId", ""))

    propsTxt <- extractXmlBlock(txt, "properties")
    ownProps <- character(0)

    if (nzchar(propsTxt)) {
        hits <- regmatches(
            propsTxt,
            gregexpr("<([A-Za-z0-9._-]+)>([^<]*)</\\1>", propsTxt, perl = TRUE)
        )[[1]]

        if (length(hits)) {
            names(hits) <- sub("^<([A-Za-z0-9._-]+)>.*$", "\\1", hits, perl = TRUE)
            ownProps    <- trimws(
                sub("^<[A-Za-z0-9._-]+>(.*)</[A-Za-z0-9._-]+>$", "\\1", hits, perl = TRUE)
            )
            names(ownProps) <- names(hits)
        }
    }

    inherited <- parent$properties[
        !names(parent$properties) %in% MAVEN_BUILTIN_PROPS
    ]

    builtIns <- c(
        "project.version" = ownVersion,
        "pom.version"     = ownVersion,
        "version"         = ownVersion,
        "project.groupId" = ownGroup,
        "pom.groupId"     = ownGroup
    )

    # Built-ins first so they take precedence, then own declarations, then
    # anything inherited from the parent chain.
    properties <- c(builtIns, ownProps, inherited)
    properties <- properties[!duplicated(names(properties))]

    # --- dependency management -------------------------------------
    mgmtTxt <- extractXmlBlock(txt, "dependencyManagement")
    managed <- parent$managed

    if (nzchar(mgmtTxt)) {
        mgmtDeps <- parseDependencyBlocks(mgmtTxt, properties)

        for (dep in mgmtDeps) {
            if (!nzchar(dep$groupId) || !nzchar(dep$artifactId)) next

            # A BOM import contributes the managed versions of another POM
            if (identical(dep$scope, "import") && identical(dep$type, "pom")) {
                bom <- parsePom(
                    fetchPom(dep$groupId, dep$artifactId, dep$version, cache),
                    cache,
                    depth + 1L
                )
                imported <- bom$managed[!names(bom$managed) %in% names(managed)]
                managed  <- c(managed, imported)
                next
            }

            if (nzchar(dep$version)) {
                managed[[depKey(dep$groupId, dep$artifactId)]] <- dep$version
            }
        }
    }

    # --- declared dependencies -------------------------------------
    withoutMgmt <- gsub(
        "(?s)<dependencyManagement>.*?</dependencyManagement>", "",
        txt, perl = TRUE
    )
    depsTxt <- extractXmlBlock(withoutMgmt, "dependencies")

    dependencies <- if (nzchar(depsTxt)) {
        parseDependencyBlocks(depsTxt, properties)
    } else {
        list()
    }

    # Inherit any version the POM did not state explicitly
    dependencies <- lapply(dependencies, function(dep) {
        if (!nzchar(dep$version)) {
            key <- depKey(dep$groupId, dep$artifactId)
            if (key %in% names(managed)) dep$version <- managed[[key]]
        }
        dep
    })

    list(
        properties   = properties,
        managed      = managed,
        dependencies = c(dependencies, parent$dependencies)
    )
}


## ----
# Decide whether a dependency belongs on the runtime classpath
#
# @param dep A dependency record.
# @param exclusions Character vector of inherited exclusion keys.
#
# @return \code{TRUE} if the dependency should be resolved.
#
# @keywords internal
isRuntimeDependency <- function(dep, exclusions = character(0)) {
    if (!nzchar(dep$groupId) || !nzchar(dep$artifactId)) return(FALSE)
    if (!nzchar(dep$version))                            return(FALSE)
    if (dep$optional)                                    return(FALSE)
    if (dep$scope %in% MAVEN_SKIPPED_SCOPES)             return(FALSE)
    if (!identical(dep$type, "jar"))                     return(FALSE)

    key <- depKey(dep$groupId, dep$artifactId)
    if (key %in% exclusions)                             return(FALSE)

    # Wildcard exclusions, as in <groupId>org.foo</groupId><artifactId>*</...>
    wildcards <- c(
        depKey(dep$groupId, "*"),
        depKey("*", dep$artifactId),
        "*:*"
    )
    if (any(wildcards %in% exclusions))                  return(FALSE)

    TRUE
}


## ----
# Drop artifacts that a newer coordinate already supersedes
#
# @param deps A data frame of resolved coordinates.
#
# @return \code{deps} without any superseded rows.
#
# @keywords internal
dropSupersededArtifacts <- function(deps) {
    if (!nrow(deps)) return(deps)

    keys        <- depKey(deps$groupId, deps$artifactId)
    replacement <- unname(MAVEN_SUPERSEDED[keys])
    obsolete    <- !is.na(replacement) & replacement %in% keys

    if (!any(obsolete)) return(deps)

    deps <- deps[!obsolete, , drop = FALSE]
    row.names(deps) <- NULL

    deps
}


## ----
#' @title Resolve the transitive runtime dependencies of a Maven artifact
#'
#' @description
#' Performs a breadth-first walk of the dependency graph rooted at the
#' given coordinates, applying Maven's nearest-wins conflict rule: the
#' version encountered at the shallowest depth is the one that is kept.
#'
#' @param groupId Maven group identifier.
#' @param artifactId Maven artifact identifier.
#' @param version Artifact version.
#' @param verbose Whether to report progress. Defaults to \code{TRUE}.
#'
#' @return
#' A data frame with columns \code{groupId}, \code{artifactId}, and
#' \code{version}, excluding the root artifact itself.
#'
#' @keywords internal
resolveDependencies <- function(groupId, artifactId, version, verbose = TRUE) {
    pomCache <- new.env(parent = emptyenv())
    rootPom  <- fetchPom(groupId, artifactId, version, pomCache)

    if (is.null(rootPom)) {
        cli::cli_abort(c(
            "Could not retrieve the POM for {groupId}:{artifactId}:{version}",
            "i" = "Checked: {.val {unname(TASSEL_MAVEN$REPOS)}}"
        ))
    }

    if (verbose) cli::cli_alert_info("Resolving dependencies...")

    resolved <- list()

    queue <- lapply(
        Filter(
            function(d) isRuntimeDependency(d),
            parsePom(rootPom, pomCache)$dependencies
        ),
        function(d) list(dep = d, depth = 1L, exclusions = d$exclusions)
    )

    while (length(queue)) {
        item  <- queue[[1]]
        queue <- queue[-1]

        dep <- item$dep
        key <- depKey(dep$groupId, dep$artifactId)

        # Nearest-wins: an entry recorded at an equal or shallower depth
        # already dictates the version for this coordinate.
        if (!is.null(resolved[[key]]) && resolved[[key]]$depth <= item$depth) next

        resolved[[key]] <- list(
            groupId    = dep$groupId,
            artifactId = dep$artifactId,
            version    = dep$version,
            depth      = item$depth
        )

        childPom <- fetchPom(dep$groupId, dep$artifactId, dep$version, pomCache)
        if (is.null(childPom)) next

        children   <- parsePom(childPom, pomCache)$dependencies
        exclusions <- unique(c(item$exclusions, dep$exclusions))

        for (child in children) {
            if (!isRuntimeDependency(child, exclusions)) next

            queue[[length(queue) + 1L]] <- list(
                dep        = child,
                depth      = item$depth + 1L,
                exclusions = unique(c(exclusions, child$exclusions))
            )
        }
    }

    if (!length(resolved)) {
        return(
            data.frame(
                groupId    = character(0),
                artifactId = character(0),
                version    = character(0),
                stringsAsFactors = FALSE
            )
        )
    }

    result <- data.frame(
        groupId    = vapply(resolved, `[[`, character(1), "groupId"),
        artifactId = vapply(resolved, `[[`, character(1), "artifactId"),
        version    = vapply(resolved, `[[`, character(1), "version"),
        stringsAsFactors = FALSE,
        row.names  = NULL
    )

    result <- result[order(result$groupId, result$artifactId), , drop = FALSE]
    row.names(result) <- NULL
    result <- dropSupersededArtifacts(result)

    if (verbose) {
        cli::cli_alert_success("Resolved {nrow(result)} dependenc{?y/ies}")
    }

    result
}


## ----
# Fetch the published SHA-1 checksum for an artifact
#
# @param url Full URL of the artifact.
#
# @return A 40-character SHA-1 string, or \code{NULL} if unavailable.
#
# @keywords internal
fetchArtifactSha1 <- function(url) {
    txt <- tryCatch(
        fetchUrlText(paste0(url, ".sha1"), MAVEN_FETCH_TIMEOUT),
        error = function(e) NULL
    )

    if (is.null(txt)) return(NULL)

    match <- regmatches(txt, regexpr("[0-9a-fA-F]{40}", txt))

    if (!length(match)) return(NULL)

    tolower(match)
}


## ----
# Download a URL to a local file without progress reporting
#
# @param url The URL to download from.
# @param destfile Local path to write to.
#
# @return The \code{destfile} path, invisibly.
#
# @keywords internal
downloadQuietly <- function(url, destfile) {
    oldTimeout <- getOption("timeout")
    options(timeout = MAVEN_FETCH_TIMEOUT)
    on.exit(options(timeout = oldTimeout), add = TRUE)

    utils::download.file(url, destfile, mode = "wb", quiet = TRUE)

    invisible(destfile)
}


## ----
# Download resolved dependency JARs into a directory
#
# @description
# Each artifact is verified against its published SHA-1 checksum. An
# artifact that cannot be found or fails verification produces a warning
# and is skipped, so that a single bad coordinate cannot abort an
# otherwise usable installation.
#
# @param deps A data frame as returned by \code{resolveDependencies}.
# @param destDir Directory to download into.
#
# @return A character vector of downloaded file paths, invisibly.
#
# @keywords internal
downloadArtifacts <- function(deps, destDir) {
    if (!nrow(deps)) return(invisible(character(0)))

    dir.create(destDir, recursive = TRUE, showWarnings = FALSE)

    downloaded <- character(0)
    failed     <- character(0)

    cli::cli_progress_bar(
        format = paste0(
            "{cli::pb_spin} Downloading dependencies ",
            "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} ",
            "[{cli::pb_elapsed}]"
        ),
        total  = nrow(deps),
        clear  = FALSE
    )

    for (i in seq_len(nrow(deps))) {
        dep      <- deps[i, ]
        coord    <- sprintf("%s:%s:%s", dep$groupId, dep$artifactId, dep$version)
        path     <- mavenArtifactPath(dep$groupId, dep$artifactId, dep$version)
        destFile <- file.path(destDir, basename(path))

        cli::cli_progress_update(set = i)

        if (file.exists(destFile)) {
            downloaded <- c(downloaded, destFile)
            next
        }

        url <- findArtifactUrl(path)

        if (is.null(url)) {
            failed <- c(failed, coord)
            next
        }

        ok <- tryCatch({
            downloadQuietly(url, destFile)

            expected <- fetchArtifactSha1(url)
            if (!is.null(expected)) {
                actual <- digest::digest(destFile, algo = "sha1", file = TRUE)
                if (!identical(actual, expected)) {
                    stop("checksum mismatch", call. = FALSE)
                }
            }

            TRUE
        }, error = function(e) FALSE)

        if (ok) {
            downloaded <- c(downloaded, destFile)
        } else {
            if (file.exists(destFile)) file.remove(destFile)
            failed <- c(failed, coord)
        }
    }

    cli::cli_progress_done()

    if (length(failed)) {
        cli::cli_warn(c(
            "!" = "Could not download {length(failed)} dependenc{?y/ies}",
            "x" = "{failed}",
            "i" = "TASSEL may fail if these are needed at runtime"
        ))
    }

    invisible(downloaded)
}
