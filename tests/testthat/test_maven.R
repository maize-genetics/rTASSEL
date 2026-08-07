# === Maven metadata, update checks, and dependency resolution =======

# --- HTTP helpers -------------------------------------------------

# A file:// URL for a local path, in the form base::url() expects on
# every platform
fileUrl <- function(path) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)

    paste0("file://", if (startsWith(path, "/")) "" else "/", path)
}

test_that("httpStatus reports the status carried by the response headers", {
    local_mocked_bindings(
        curlGetHeaders = function(url, ...) {
            structure("HTTP/1.1 200 OK", status = 200L)
        },
        .package = "base"
    )

    expect_equal(httpStatus("https://example.org/x.jar"), 200L)
})

test_that("httpStatus is NA when the response carries no status", {
    local_mocked_bindings(
        curlGetHeaders = function(...) "nothing useful",
        .package = "base"
    )

    expect_equal(httpStatus("https://example.org/x.jar"), NA_integer_)
})

test_that("httpStatus is NA rather than an error when the request fails", {
    local_mocked_bindings(
        curlGetHeaders = function(...) stop("could not resolve host"),
        .package = "base"
    )

    expect_equal(httpStatus("https://example.org/x.jar"), NA_integer_)
})

test_that("httpStatus passes the timeout on as an integer", {
    requested <- NULL

    local_mocked_bindings(
        curlGetHeaders = function(url, timeout, ...) {
            requested <<- timeout
            structure("", status = 404L)
        },
        .package = "base"
    )

    httpStatus("https://example.org/x.jar", timeout = 2.7)

    expect_identical(requested, 2L)
})

test_that("urlExists is TRUE only for HTTP 200", {
    local_mocked_bindings(httpStatus = function(...) 200L)

    expect_true(urlExists("https://example.org/x.jar"))
})

test_that("urlExists is FALSE for any other status", {
    local_mocked_bindings(httpStatus = function(...) 404L)

    expect_false(urlExists("https://example.org/x.jar"))
})

test_that("urlExists is FALSE when the status could not be determined", {
    local_mocked_bindings(httpStatus = function(...) NA_integer_)

    expect_false(urlExists("https://example.org/x.jar"))
})

test_that("fetchUrlText returns the whole response body as one string", {
    payload <- withr::local_tempfile(fileext = ".xml")
    writeLines(c("<metadata>", "  <release>5.2.96</release>", "</metadata>"), payload)

    expect_equal(
        fetchUrlText(fileUrl(payload)),
        "<metadata>\n  <release>5.2.96</release>\n</metadata>"
    )
})

test_that("fetchUrlText restores the caller's timeout option", {
    payload <- withr::local_tempfile(fileext = ".xml")
    writeLines("<metadata/>", payload)

    withr::local_options(timeout = 13L)

    fetchUrlText(fileUrl(payload), timeout = 1L)

    expect_equal(getOption("timeout"), 13L)
})


# --- XML helpers --------------------------------------------------

test_that("extractXmlTags pulls every occurrence of a tag", {
    txt <- "<versions><version>1.0</version><version>2.0</version></versions>"

    expect_equal(extractXmlTags(txt, "version"), c("1.0", "2.0"))
})

test_that("extractXmlTags returns empty for an absent tag", {
    expect_length(extractXmlTags("<a>1</a>", "version"), 0)
})

test_that("extractXmlTags does not confuse a tag with a longer one", {
    txt <- "<versioning><version>1.0</version></versioning>"

    expect_equal(extractXmlTags(txt, "version"), "1.0")
})

test_that("extractXmlBlock captures multi-line content", {
    txt <- "<parent>\n  <groupId>g</groupId>\n  <version>1</version>\n</parent>"

    expect_match(extractXmlBlock(txt, "parent"), "groupId")
    expect_equal(extractXmlBlock(txt, "missing"), "")
})

test_that("extractXmlBlocks captures repeated elements", {
    txt <- "<deps><dependency>\na\n</dependency><dependency>\nb\n</dependency></deps>"

    expect_length(extractXmlBlocks(txt, "dependency"), 2)
})

test_that("extractXmlBlocks returns empty for an absent element", {
    expect_length(extractXmlBlocks("<deps></deps>", "dependency"), 0)
})


# --- version comparison -------------------------------------------

test_that("sortVersionsDesc orders numerically, not lexically", {
    versions <- c("5.2.9", "5.2.100", "5.2.96")

    expect_equal(sortVersionsDesc(versions), c("5.2.100", "5.2.96", "5.2.9"))
})

test_that("sortVersionsDesc drops non-numeric versions", {
    versions <- c("5.2.96", "5.3-SNAPSHOT", "5.2.97-RC1")

    expect_equal(sortVersionsDesc(versions), "5.2.96")
})

test_that("sortVersionsDesc returns empty when nothing is a released version", {
    expect_length(sortVersionsDesc(c("5.3-SNAPSHOT", "5.2.97-RC1")), 0)
    expect_length(sortVersionsDesc(character(0)), 0)
})

test_that("isNewerVersion compares released versions", {
    expect_true(isNewerVersion("5.2.97", "5.2.96"))
    expect_false(isNewerVersion("5.2.96", "5.2.96"))
    expect_false(isNewerVersion("5.2.95", "5.2.96"))
})

test_that("isNewerVersion is FALSE rather than an error on junk input", {
    expect_false(isNewerVersion("not-a-version", "5.2.96"))
    expect_false(isNewerVersion("5.2.96", "not-a-version"))
})


# --- maven-metadata.xml -------------------------------------------

test_that("fetchMavenMetadata parses release and versions", {
    metadata <- paste(
        "<metadata><versioning>",
        "<latest>5.2.97</latest><release>5.2.97</release>",
        "<versions><version>5.2.95</version><version>5.2.96</version>",
        "<version>5.2.97</version></versions>",
        "</versioning></metadata>",
        sep = ""
    )

    local_mocked_bindings(fetchUrlText = function(...) metadata)

    result <- fetchMavenMetadata()

    expect_equal(result$release, "5.2.97")
    expect_equal(result$versions, c("5.2.97", "5.2.96", "5.2.95"))
})


# --- artifact layout probing --------------------------------------

test_that("probeArtifactLayout reports 'fat' when the fat JAR exists", {
    local_mocked_bindings(urlExists = function(url, ...) grepl("jar-with-dependencies", url))

    expect_equal(probeArtifactLayout("5.2.96"), "fat")
})

test_that("probeArtifactLayout reports 'thin' for a JAR and POM pair", {
    local_mocked_bindings(
        urlExists = function(url, ...) !grepl("jar-with-dependencies", url)
    )

    expect_equal(probeArtifactLayout("5.2.97"), "thin")
})

test_that("probeArtifactLayout reports 'none' when nothing usable exists", {
    local_mocked_bindings(urlExists = function(...) FALSE)

    expect_equal(probeArtifactLayout("9.9.9"), "none")
})

test_that("probeArtifactLayout requires a POM alongside a thin JAR", {
    # A thin JAR with no POM cannot have its dependencies resolved
    local_mocked_bindings(
        urlExists = function(url, ...) {
            grepl("\\.jar$", url) && !grepl("jar-with-dependencies", url)
        }
    )

    expect_equal(probeArtifactLayout("5.2.97"), "none")
})


# --- newest installable version -----------------------------------

test_that("latestInstallableTASSEL skips versions that cannot be installed", {
    local_mocked_bindings(
        fetchMavenMetadata  = function(...) list(
            release  = "5.2.98",
            versions = c("5.2.98", "5.2.97", "5.2.96")
        ),
        probeArtifactLayout = function(version, ...) {
            if (identical(version, "5.2.96")) "fat" else "none"
        }
    )

    result <- latestInstallableTASSEL()

    expect_equal(result$latest, "5.2.96")
    expect_equal(result$layout, "fat")
})

test_that("latestInstallableTASSEL returns NULL when metadata lists no versions", {
    probed <- FALSE

    local_mocked_bindings(
        fetchMavenMetadata  = function(...) list(
            release  = NULL,
            versions = character(0)
        ),
        probeArtifactLayout = function(...) {
            probed <<- TRUE
            "fat"
        }
    )

    expect_null(latestInstallableTASSEL())
    expect_false(probed)
})

test_that("latestInstallableTASSEL returns NULL when nothing is installable", {
    local_mocked_bindings(
        fetchMavenMetadata  = function(...) list(release = NULL, versions = "5.2.97"),
        probeArtifactLayout = function(...) "none"
    )

    expect_null(latestInstallableTASSEL())
})

test_that("latestInstallableTASSEL bounds how many versions it probes", {
    probed <- 0L

    local_mocked_bindings(
        fetchMavenMetadata  = function(...) list(
            release  = NULL,
            versions = sprintf("5.2.%d", 100:1)
        ),
        probeArtifactLayout = function(...) {
            probed <<- probed + 1L
            "none"
        }
    )

    latestInstallableTASSEL(maxProbe = 3L)

    expect_equal(probed, 3L)
})


# --- update check gating ------------------------------------------

test_that("updateCheckEnabled is FALSE in non-interactive sessions", {
    # testthat always runs non-interactively
    expect_false(updateCheckEnabled())
})

test_that("updateCheckEnabled honours the opt-out option", {
    rlang::local_interactive(TRUE)
    withr::local_options(rTASSEL.check.updates = FALSE)

    expect_false(updateCheckEnabled())
})

test_that("updateCheckEnabled honours the opt-out environment variable", {
    rlang::local_interactive(TRUE)
    withr::local_envvar(RTASSEL_NO_VERSION_CHECK = "true", CI = "")

    expect_false(updateCheckEnabled())
})

test_that("updateCheckEnabled is FALSE on continuous integration", {
    rlang::local_interactive(TRUE)
    withr::local_envvar(RTASSEL_NO_VERSION_CHECK = "", CI = "true")

    expect_false(updateCheckEnabled())
})

test_that("updateCheckEnabled is TRUE in an interactive session with no opt-outs", {
    rlang::local_interactive(TRUE)
    withr::local_options(rTASSEL.check.updates = TRUE)

    # R CMD check exports its own settings as environment variables, which
    # would otherwise suppress the check no matter what the test asks for
    envvars <- list(RTASSEL_NO_VERSION_CHECK = "", CI = "")
    for (name in grep("^_R_CHECK_", names(Sys.getenv()), value = TRUE)) {
        envvars[[name]] <- NA
    }
    withr::local_envvar(envvars)

    expect_true(updateCheckEnabled())
})


# --- update check caching -----------------------------------------

test_that("checkForTASSELUpdate returns NULL when checks are disabled", {
    local_mocked_bindings(updateCheckEnabled = function() FALSE)

    expect_null(checkForTASSELUpdate())
})

test_that("checkForTASSELUpdate uses a fresh cache without hitting the network", {
    called <- FALSE

    local_mocked_bindings(
        updateCheckEnabled     = function() TRUE,
        readUpdateCheckCache   = function(...) list(
            latest    = "5.2.96",
            layout    = "fat",
            checkedAt = Sys.time()
        ),
        latestInstallableTASSEL = function(...) {
            called <<- TRUE
            list(latest = "5.2.97", layout = "thin")
        }
    )

    result <- checkForTASSELUpdate()

    expect_false(called)
    expect_equal(result$latest, "5.2.96")
})

test_that("checkForTASSELUpdate refreshes a stale cache", {
    stale <- Sys.time() - (TASSEL_UPDATE$MAX_AGE_SECS + 60)
    written <- NULL

    local_mocked_bindings(
        updateCheckEnabled      = function() TRUE,
        readUpdateCheckCache    = function(...) list(
            latest    = "5.2.96",
            layout    = "fat",
            checkedAt = stale
        ),
        writeUpdateCheckCache   = function(result, ...) {
            written <<- result
            invisible(result)
        },
        latestInstallableTASSEL = function(...) list(latest = "5.2.97", layout = "thin")
    )

    result <- checkForTASSELUpdate()

    expect_equal(result$latest, "5.2.97")
    expect_equal(written$latest, "5.2.97")
    expect_s3_class(result$checkedAt, "POSIXct")
})

test_that("checkForTASSELUpdate falls back to a stale result when offline", {
    stale <- Sys.time() - (TASSEL_UPDATE$MAX_AGE_SECS + 60)

    local_mocked_bindings(
        updateCheckEnabled      = function() TRUE,
        readUpdateCheckCache    = function(...) list(
            latest    = "5.2.96",
            layout    = "fat",
            checkedAt = stale
        ),
        latestInstallableTASSEL = function(...) stop("no network")
    )

    result <- checkForTASSELUpdate()

    expect_equal(result$latest, "5.2.96")
})

test_that("checkForTASSELUpdate returns NULL when offline with no cache", {
    local_mocked_bindings(
        updateCheckEnabled      = function() TRUE,
        readUpdateCheckCache    = function(...) NULL,
        latestInstallableTASSEL = function(...) stop("no network")
    )

    expect_null(checkForTASSELUpdate())
})

test_that("updateCheckCachePath points into the user cache directory", {
    expect_equal(
        updateCheckCachePath(),
        file.path(
            tools::R_user_dir("rTASSEL", "cache"),
            TASSEL_UPDATE$CACHE_FILE
        )
    )
})

test_that("updateCheckCachePath names the file for the requested channel", {
    expect_equal(
        basename(updateCheckCachePath(TASSEL_GITHUB$CACHE_FILE)),
        TASSEL_GITHUB$CACHE_FILE
    )
})

test_that("readUpdateCheckCache returns NULL when no cache has been written", {
    cacheFile <- file.path(withr::local_tempdir(), "version-check.rds")

    local_mocked_bindings(updateCheckCachePath = function(...) cacheFile)

    expect_null(readUpdateCheckCache())
})

test_that("update check cache survives a write and read round trip", {
    cacheFile <- withr::local_tempfile(fileext = ".rds")

    local_mocked_bindings(updateCheckCachePath = function(...) cacheFile)

    payload <- list(latest = "5.2.96", layout = "fat", checkedAt = Sys.time())
    writeUpdateCheckCache(payload)

    expect_equal(readUpdateCheckCache()$latest, "5.2.96")
})

test_that("readUpdateCheckCache tolerates a corrupt cache file", {
    cacheFile <- withr::local_tempfile(fileext = ".rds")
    writeLines("not an rds file", cacheFile)

    local_mocked_bindings(updateCheckCachePath = function(...) cacheFile)

    expect_null(readUpdateCheckCache())
})


# --- update check channels -----------------------------------------

test_that("checkForTASSELUpdate queries GitHub for the nightly channel", {
    local_mocked_bindings(
        updateCheckEnabled      = function() TRUE,
        readUpdateCheckCache    = function(...) NULL,
        writeUpdateCheckCache   = function(result, ...) invisible(result),
        latestInstallableTASSEL = function(...) stop("maven must not be queried"),
        latestNightlyTASSEL     = function(...) list(
            latest      = "5.2.98-dev.20260801",
            tag         = "dev-20260801",
            publishedAt = "2026-08-01T07:45:34Z",
            channel     = "nightly"
        )
    )

    result <- checkForTASSELUpdate(channel = "nightly")

    expect_equal(result$latest, "5.2.98-dev.20260801")
    expect_equal(result$channel, "nightly")
    expect_s3_class(result$checkedAt, "POSIXct")
})

test_that("checkForTASSELUpdate caches each channel separately", {
    requested <- character(0)

    local_mocked_bindings(
        updateCheckEnabled      = function() TRUE,
        readUpdateCheckCache    = function(file = TASSEL_UPDATE$CACHE_FILE) {
            requested <<- c(requested, file)
            NULL
        },
        writeUpdateCheckCache   = function(result, ...) invisible(result),
        latestInstallableTASSEL = function(...) list(latest = "5.2.97", layout = "thin"),
        latestNightlyTASSEL     = function(...) list(latest = "5.2.98-dev.20260801")
    )

    checkForTASSELUpdate()
    checkForTASSELUpdate(channel = "nightly")

    expect_equal(requested, c(TASSEL_UPDATE$CACHE_FILE, TASSEL_GITHUB$CACHE_FILE))
})

test_that("checkForTASSELUpdate rejects an unknown channel", {
    expect_error(checkForTASSELUpdate(channel = "snapshot"))
})


# --- coordinate helpers -------------------------------------------

test_that("groupIdToPath converts dots to slashes", {
    expect_equal(groupIdToPath("net.maizegenetics"), "net/maizegenetics")
})

test_that("mavenArtifactPath builds a repository-relative path", {
    expect_equal(
        mavenArtifactPath("com.google.guava", "guava", "22.0"),
        "com/google/guava/guava/22.0/guava-22.0.jar"
    )
})

test_that("mavenArtifactPath honours extension and classifier", {
    expect_equal(
        mavenArtifactPath("net.maizegenetics", "tassel", "5.2.96", ext = "pom"),
        "net/maizegenetics/tassel/5.2.96/tassel-5.2.96.pom"
    )
    expect_equal(
        mavenArtifactPath(
            "net.maizegenetics", "tassel", "5.2.96",
            classifier = "jar-with-dependencies"
        ),
        paste0(
            "net/maizegenetics/tassel/5.2.96/",
            "tassel-5.2.96-jar-with-dependencies.jar"
        )
    )
})


# --- repository probing -------------------------------------------

test_that("findArtifactUrl returns the first repository that serves the artifact", {
    local_mocked_bindings(urlExists = function(url, ...) grepl("scijava", url))

    expect_equal(
        findArtifactUrl("cisd/jhdf5/19.04.0/jhdf5-19.04.0.jar"),
        paste0(
            TASSEL_MAVEN$REPOS[["scijava"]],
            "/cisd/jhdf5/19.04.0/jhdf5-19.04.0.jar"
        )
    )
})

test_that("findArtifactUrl prefers repositories in declaration order", {
    probed <- character(0)

    local_mocked_bindings(
        urlExists = function(url, ...) {
            probed <<- c(probed, url)
            TRUE
        }
    )

    result <- findArtifactUrl("x/y/1/y-1.jar")

    expect_length(probed, 1)
    expect_equal(result, paste0(TASSEL_MAVEN$REPOS[["central"]], "/x/y/1/y-1.jar"))
})

test_that("findArtifactUrl returns NULL when no repository serves the artifact", {
    local_mocked_bindings(urlExists = function(...) FALSE)

    expect_null(findArtifactUrl("x/y/1/y-1.jar"))
})


# --- POM fetching -------------------------------------------------

test_that("fetchPom requests the POM coordinate from Maven Central first", {
    requested <- character(0)

    local_mocked_bindings(
        fetchUrlText = function(url, ...) {
            requested <<- c(requested, url)
            "<project/>"
        }
    )

    expect_equal(fetchPom("com.google.guava", "guava", "22.0"), "<project/>")
    expect_equal(
        requested,
        paste0(
            TASSEL_MAVEN$REPOS[["central"]],
            "/com/google/guava/guava/22.0/guava-22.0.pom"
        )
    )
})

test_that("fetchPom falls through to the next repository on a miss", {
    local_mocked_bindings(
        fetchUrlText = function(url, ...) {
            if (!grepl("scijava", url)) stop("404")
            "<project>scijava</project>"
        }
    )

    expect_equal(fetchPom("cisd", "jhdf5", "19.04.0"), "<project>scijava</project>")
})

test_that("fetchPom treats an empty response as a miss", {
    local_mocked_bindings(fetchUrlText = function(...) "")

    expect_null(fetchPom("x", "y", "1"))
})

test_that("fetchPom returns NULL when no repository serves the POM", {
    local_mocked_bindings(fetchUrlText = function(...) stop("404"))

    expect_null(fetchPom("x", "y", "1"))
})

test_that("fetchPom memoises a POM it has already retrieved", {
    calls <- 0L
    cache <- new.env(parent = emptyenv())

    local_mocked_bindings(
        fetchUrlText = function(...) {
            calls <<- calls + 1L
            "<project/>"
        }
    )

    fetchPom("x", "y", "1", cache)
    fetchPom("x", "y", "1", cache)

    expect_equal(calls, 1L)
})

test_that("fetchPom memoises a missing POM so it is only probed once", {
    calls <- 0L
    cache <- new.env(parent = emptyenv())

    local_mocked_bindings(
        fetchUrlText = function(...) {
            calls <<- calls + 1L
            stop("404")
        }
    )

    expect_null(fetchPom("x", "y", "1", cache))
    expect_null(fetchPom("x", "y", "1", cache))

    expect_equal(calls, length(TASSEL_MAVEN$REPOS))
})

test_that("fetchPom keys its cache on the full coordinate", {
    cache <- new.env(parent = emptyenv())

    local_mocked_bindings(
        fetchUrlText = function(url, ...) sprintf("<project>%s</project>", basename(url))
    )

    expect_match(fetchPom("x", "y", "1", cache), "y-1.pom")
    expect_match(fetchPom("x", "y", "2", cache), "y-2.pom")
})


# --- property interpolation ---------------------------------------

test_that("interpolateProperties substitutes a placeholder", {
    expect_equal(
        interpolateProperties("${a}", c(a = "1.0")),
        "1.0"
    )
})

test_that("interpolateProperties resolves chained properties", {
    expect_equal(
        interpolateProperties("${a}", c(a = "${b}", b = "2.0")),
        "2.0"
    )
})

test_that("interpolateProperties leaves unknown placeholders alone", {
    expect_equal(interpolateProperties("${nope}", c(a = "1")), "${nope}")
})

test_that("interpolateProperties terminates on a circular reference", {
    expect_equal(
        interpolateProperties("${a}", c(a = "${b}", b = "${a}")),
        "${a}"
    )
})

test_that("interpolateProperties leaves a malformed placeholder untouched", {
    expect_equal(interpolateProperties("${}", c(a = "1")), "${}")
    expect_equal(interpolateProperties("1.0-${", c(a = "1")), "1.0-${")
})

test_that("interpolateProperties passes through values with nothing to do", {
    expect_null(interpolateProperties(NULL, c(a = "1")))
    expect_equal(interpolateProperties("", c(a = "1")), "")
    expect_equal(interpolateProperties("1.0", c(a = "1")), "1.0")
})


# --- POM header isolation -----------------------------------------

test_that("pomHeader ignores the parent block when reading own coordinates", {
    pom <- paste(
        "<project>",
        "<parent><groupId>g</groupId><artifactId>p</artifactId>",
        "<version>99</version></parent>",
        "<groupId>own</groupId><artifactId>a</artifactId><version>1.8</version>",
        "<dependencies><dependency><version>7.7</version></dependency></dependencies>",
        "</project>",
        sep = ""
    )

    header <- pomHeader(pom)

    expect_equal(xmlValue(header, "version"), "1.8")
    expect_equal(xmlValue(header, "groupId"), "own")
})

test_that("pomHeader excludes dependency versions", {
    pom <- paste(
        "<project><groupId>own</groupId><artifactId>a</artifactId>",
        "<dependencies><dependency><version>7.7</version></dependency></dependencies>",
        "</project>",
        sep = ""
    )

    expect_null(xmlValue(pomHeader(pom), "version"))
})


# --- POM parsing --------------------------------------------------

test_that("parsePom returns an empty result for a POM that could not be read", {
    for (txt in list(NULL, "")) {
        parsed <- parsePom(txt)

        expect_length(parsed$properties, 0)
        expect_length(parsed$managed, 0)
        expect_length(parsed$dependencies, 0)
    }
})

test_that("parsePom bounds how far it walks a parent chain", {
    # A POM naming itself as its own parent would otherwise recurse forever
    selfPom <- paste(
        "<project>",
        "<parent><groupId>g</groupId><artifactId>a</artifactId>",
        "<version>1</version></parent>",
        "<groupId>g</groupId><artifactId>a</artifactId>",
        "<dependencies><dependency><groupId>x</groupId>",
        "<artifactId>y</artifactId><version>1</version>",
        "</dependency></dependencies></project>",
        sep = ""
    )

    local_mocked_bindings(fetchPom = function(...) selfPom)

    # One dependency per visited POM: the root plus ten ancestors
    expect_length(parsePom(selfPom)$dependencies, 11)
})

test_that("parsePom reads declared dependencies", {
    pom <- paste(
        "<project><groupId>g</groupId><artifactId>a</artifactId>",
        "<version>1.0</version><dependencies>",
        "<dependency><groupId>x</groupId><artifactId>y</artifactId>",
        "<version>2.0</version></dependency>",
        "</dependencies></project>",
        sep = ""
    )

    parsed <- parsePom(pom)

    expect_length(parsed$dependencies, 1)
    expect_equal(parsed$dependencies[[1]]$artifactId, "y")
    expect_equal(parsed$dependencies[[1]]$version, "2.0")
})

test_that("parsePom interpolates dependency versions from properties", {
    pom <- paste(
        "<project><groupId>g</groupId><artifactId>a</artifactId>",
        "<version>1.0</version>",
        "<properties><yVersion>3.1</yVersion></properties>",
        "<dependencies><dependency><groupId>x</groupId>",
        "<artifactId>y</artifactId><version>${yVersion}</version>",
        "</dependency></dependencies></project>",
        sep = ""
    )

    expect_equal(parsePom(pom)$dependencies[[1]]$version, "3.1")
})

test_that("parsePom fills a missing version from dependencyManagement", {
    pom <- paste(
        "<project><groupId>g</groupId><artifactId>a</artifactId>",
        "<version>1.0</version>",
        "<dependencyManagement><dependencies><dependency>",
        "<groupId>x</groupId><artifactId>y</artifactId><version>4.2</version>",
        "</dependency></dependencies></dependencyManagement>",
        "<dependencies><dependency><groupId>x</groupId>",
        "<artifactId>y</artifactId></dependency></dependencies></project>",
        sep = ""
    )

    expect_equal(parsePom(pom)$dependencies[[1]]$version, "4.2")
})

test_that("parsePom ignores dependencies declared inside build and profiles", {
    pom <- paste(
        "<project><groupId>g</groupId><artifactId>a</artifactId>",
        "<version>1.0</version>",
        "<build><dependencies><dependency><groupId>bad</groupId>",
        "<artifactId>bad</artifactId><version>1</version>",
        "</dependency></dependencies></build>",
        "<profiles><dependencies><dependency><groupId>worse</groupId>",
        "<artifactId>worse</artifactId><version>1</version>",
        "</dependency></dependencies></profiles>",
        "</project>",
        sep = ""
    )

    expect_length(parsePom(pom)$dependencies, 0)
})

test_that("parsePom keeps ${project.version} local to the POM being parsed", {
    # Regression: an ancestor's version must not leak downward. The real
    # case was doxia:1.8 -> maven-parent:30 -> apache:18, where the
    # grandparent's version was resolving managed doxia coordinates to 18.
    grandPom <- paste(
        "<project><groupId>g</groupId><artifactId>grand</artifactId>",
        "<version>99</version></project>",
        sep = ""
    )

    childPom <- paste(
        "<project>",
        "<parent><groupId>g</groupId><artifactId>grand</artifactId>",
        "<version>99</version></parent>",
        "<groupId>p</groupId><artifactId>child</artifactId><version>1.8</version>",
        "<dependencyManagement><dependencies><dependency>",
        "<groupId>p</groupId><artifactId>sink</artifactId>",
        "<version>${project.version}</version>",
        "</dependency></dependencies></dependencyManagement>",
        "</project>",
        sep = ""
    )

    local_mocked_bindings(fetchPom = function(...) grandPom)

    expect_equal(parsePom(childPom)$managed[["p:sink"]], "1.8")
})

test_that("parsePom imports managed versions from a BOM", {
    bomPom <- paste(
        "<project><groupId>b</groupId><artifactId>bom</artifactId>",
        "<version>1.0</version>",
        "<dependencyManagement><dependencies><dependency>",
        "<groupId>x</groupId><artifactId>y</artifactId><version>9.9</version>",
        "</dependency></dependencies></dependencyManagement></project>",
        sep = ""
    )

    pom <- paste(
        "<project><groupId>g</groupId><artifactId>a</artifactId>",
        "<version>1.0</version>",
        "<dependencyManagement><dependencies><dependency>",
        "<groupId>b</groupId><artifactId>bom</artifactId><version>1.0</version>",
        "<type>pom</type><scope>import</scope>",
        "</dependency></dependencies></dependencyManagement></project>",
        sep = ""
    )

    local_mocked_bindings(fetchPom = function(...) bomPom)

    expect_equal(parsePom(pom)$managed[["x:y"]], "9.9")
})

test_that("parsePom skips managed entries with incomplete coordinates", {
    pom <- paste(
        "<project><groupId>g</groupId><artifactId>a</artifactId>",
        "<version>1.0</version>",
        "<dependencyManagement><dependencies>",
        "<dependency><artifactId>no-group</artifactId>",
        "<version>1.0</version></dependency>",
        "<dependency><groupId>x</groupId><artifactId>y</artifactId>",
        "<version>4.2</version></dependency>",
        "</dependencies></dependencyManagement></project>",
        sep = ""
    )

    expect_equal(names(parsePom(pom)$managed), "x:y")
})

test_that("parsePom leaves a managed entry without a version unmanaged", {
    pom <- paste(
        "<project><groupId>g</groupId><artifactId>a</artifactId>",
        "<version>1.0</version>",
        "<dependencyManagement><dependencies><dependency>",
        "<groupId>x</groupId><artifactId>y</artifactId>",
        "<scope>runtime</scope></dependency>",
        "</dependencies></dependencyManagement></project>",
        sep = ""
    )

    expect_length(parsePom(pom)$managed, 0)
})

test_that("parseExclusions reads exclusion coordinates", {
    dep <- paste(
        "<groupId>x</groupId><artifactId>y</artifactId>",
        "<exclusions><exclusion><groupId>bad</groupId>",
        "<artifactId>thing</artifactId></exclusion></exclusions>",
        sep = ""
    )

    expect_equal(parseExclusions(dep), "bad:thing")
})

test_that("parseDependencyBlocks does not mistake exclusions for coordinates", {
    txt <- paste(
        "<dependency><groupId>x</groupId><artifactId>y</artifactId>",
        "<version>1.0</version>",
        "<exclusions><exclusion><groupId>bad</groupId>",
        "<artifactId>thing</artifactId></exclusion></exclusions>",
        "</dependency>",
        sep = ""
    )

    dep <- parseDependencyBlocks(txt, character(0))[[1]]

    expect_equal(dep$groupId, "x")
    expect_equal(dep$artifactId, "y")
    expect_equal(dep$exclusions, "bad:thing")
})


# --- runtime dependency filtering ---------------------------------

makeDep <- function(...) {
    utils::modifyList(
        list(
            groupId    = "x",
            artifactId = "y",
            version    = "1.0",
            scope      = "compile",
            type       = "jar",
            classifier = "",
            optional   = FALSE,
            exclusions = character(0)
        ),
        list(...)
    )
}

test_that("isRuntimeDependency accepts a plain compile dependency", {
    expect_true(isRuntimeDependency(makeDep()))
})

test_that("isRuntimeDependency rejects non-runtime scopes", {
    for (scope in c("test", "provided", "system")) {
        expect_false(isRuntimeDependency(makeDep(scope = scope)))
    }
})

test_that("isRuntimeDependency rejects optional dependencies", {
    expect_false(isRuntimeDependency(makeDep(optional = TRUE)))
})

test_that("isRuntimeDependency rejects non-JAR artifacts", {
    expect_false(isRuntimeDependency(makeDep(type = "pom")))
})

test_that("isRuntimeDependency rejects dependencies without a version", {
    expect_false(isRuntimeDependency(makeDep(version = "")))
})

test_that("isRuntimeDependency rejects incomplete coordinates", {
    expect_false(isRuntimeDependency(makeDep(groupId = "")))
    expect_false(isRuntimeDependency(makeDep(artifactId = "")))
})

test_that("isRuntimeDependency honours exclusions", {
    expect_false(isRuntimeDependency(makeDep(), exclusions = "x:y"))
})

test_that("isRuntimeDependency honours wildcard exclusions", {
    expect_false(isRuntimeDependency(makeDep(), exclusions = "x:*"))
    expect_false(isRuntimeDependency(makeDep(), exclusions = "*:*"))
})


# --- transitive resolution ----------------------------------------

mockPomSet <- function(poms) {
    function(groupId, artifactId, version, ...) {
        poms[[paste(groupId, artifactId, version, sep = ":")]]
    }
}

simplePom <- function(groupId, artifactId, version, deps = "") {
    sprintf(
        "<project><groupId>%s</groupId><artifactId>%s</artifactId>%s%s</project>",
        groupId, artifactId,
        sprintf("<version>%s</version>", version),
        if (nzchar(deps)) sprintf("<dependencies>%s</dependencies>", deps) else ""
    )
}

depXml <- function(groupId, artifactId, version, extra = "") {
    sprintf(
        "<dependency><groupId>%s</groupId><artifactId>%s</artifactId><version>%s</version>%s</dependency>",
        groupId, artifactId, version, extra
    )
}

test_that("resolveDependencies walks transitive dependencies", {
    poms <- list(
        "r:root:1"  = simplePom("r", "root", "1", depXml("x", "a", "1")),
        "x:a:1"     = simplePom("x", "a", "1", depXml("x", "b", "1")),
        "x:b:1"     = simplePom("x", "b", "1")
    )

    local_mocked_bindings(fetchPom = mockPomSet(poms))

    result <- resolveDependencies("r", "root", "1", verbose = FALSE)

    expect_equal(result$artifactId, c("a", "b"))
})

test_that("resolveDependencies applies nearest-wins conflict resolution", {
    poms <- list(
        "r:root:1" = simplePom(
            "r", "root", "1",
            paste0(depXml("x", "a", "1"), depXml("x", "c", "1"))
        ),
        "x:a:1"    = simplePom("x", "a", "1"),
        "x:c:1"    = simplePom("x", "c", "1", depXml("x", "a", "2")),
        "x:a:2"    = simplePom("x", "a", "2")
    )

    local_mocked_bindings(fetchPom = mockPomSet(poms))

    result <- resolveDependencies("r", "root", "1", verbose = FALSE)

    expect_equal(result$version[result$artifactId == "a"], "1")
})

test_that("resolveDependencies excludes the root artifact", {
    poms <- list("r:root:1" = simplePom("r", "root", "1", depXml("x", "a", "1")),
                 "x:a:1"    = simplePom("x", "a", "1"))

    local_mocked_bindings(fetchPom = mockPomSet(poms))

    result <- resolveDependencies("r", "root", "1", verbose = FALSE)

    expect_false("root" %in% result$artifactId)
})

test_that("resolveDependencies propagates exclusions to transitives", {
    poms <- list(
        "r:root:1" = simplePom(
            "r", "root", "1",
            depXml(
                "x", "a", "1",
                "<exclusions><exclusion><groupId>x</groupId><artifactId>b</artifactId></exclusion></exclusions>"
            )
        ),
        "x:a:1"    = simplePom("x", "a", "1", depXml("x", "b", "1")),
        "x:b:1"    = simplePom("x", "b", "1")
    )

    local_mocked_bindings(fetchPom = mockPomSet(poms))

    result <- resolveDependencies("r", "root", "1", verbose = FALSE)

    expect_equal(result$artifactId, "a")
})

test_that("resolveDependencies survives a dependency with an unreachable POM", {
    poms <- list("r:root:1" = simplePom("r", "root", "1", depXml("x", "a", "1")))

    local_mocked_bindings(fetchPom = mockPomSet(poms))

    result <- resolveDependencies("r", "root", "1", verbose = FALSE)

    expect_equal(result$artifactId, "a")
})

test_that("resolveDependencies aborts when the root POM is unavailable", {
    local_mocked_bindings(fetchPom = function(...) NULL)

    expect_error(
        resolveDependencies("r", "root", "1", verbose = FALSE),
        "Could not retrieve the POM"
    )
})

test_that("resolveDependencies returns an empty frame for a leaf artifact", {
    poms <- list("r:root:1" = simplePom("r", "root", "1"))

    local_mocked_bindings(fetchPom = mockPomSet(poms))

    result <- resolveDependencies("r", "root", "1", verbose = FALSE)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 0)
})

test_that("resolveDependencies reports what it resolved when verbose", {
    poms <- list(
        "r:root:1" = simplePom("r", "root", "1", depXml("x", "a", "1")),
        "x:a:1"    = simplePom("x", "a", "1")
    )

    local_mocked_bindings(fetchPom = mockPomSet(poms))

    messages <- paste(
        capture_messages(resolveDependencies("r", "root", "1", verbose = TRUE)),
        collapse = ""
    )

    expect_match(messages, "Resolving dependencies")
    expect_match(messages, "Resolved 1 dependency")
})

test_that("resolveDependencies is silent when not verbose", {
    poms <- list(
        "r:root:1" = simplePom("r", "root", "1", depXml("x", "a", "1")),
        "x:a:1"    = simplePom("x", "a", "1")
    )

    local_mocked_bindings(fetchPom = mockPomSet(poms))

    expect_silent(resolveDependencies("r", "root", "1", verbose = FALSE))
})


# --- superseded artifacts -----------------------------------------

test_that("dropSupersededArtifacts removes google-collections beside guava", {
    deps <- data.frame(
        groupId    = c("com.google.collections", "com.google.guava"),
        artifactId = c("google-collections", "guava"),
        version    = c("1.0", "22.0"),
        stringsAsFactors = FALSE
    )

    result <- dropSupersededArtifacts(deps)

    expect_equal(result$artifactId, "guava")
})

test_that("dropSupersededArtifacts keeps google-collections when guava is absent", {
    deps <- data.frame(
        groupId    = "com.google.collections",
        artifactId = "google-collections",
        version    = "1.0",
        stringsAsFactors = FALSE
    )

    expect_equal(nrow(dropSupersededArtifacts(deps)), 1)
})

test_that("dropSupersededArtifacts accepts an empty dependency set", {
    deps <- data.frame(
        groupId    = character(0),
        artifactId = character(0),
        version    = character(0),
        stringsAsFactors = FALSE
    )

    expect_equal(nrow(dropSupersededArtifacts(deps)), 0)
})

test_that("dropSupersededArtifacts leaves unrelated artifacts alone", {
    deps <- data.frame(
        groupId    = c("org.x", "org.y"),
        artifactId = c("a", "b"),
        version    = c("1", "2"),
        stringsAsFactors = FALSE
    )

    expect_equal(nrow(dropSupersededArtifacts(deps)), 2)
})

test_that("resolveDependencies drops a superseded artifact from the graph", {
    poms <- list(
        "r:root:1" = simplePom(
            "r", "root", "1",
            paste0(
                depXml("com.google.guava", "guava", "22.0"),
                depXml("com.google.collections", "google-collections", "1.0")
            )
        ),
        "com.google.guava:guava:22.0"                     = simplePom("com.google.guava", "guava", "22.0"),
        "com.google.collections:google-collections:1.0"   = simplePom("com.google.collections", "google-collections", "1.0")
    )

    local_mocked_bindings(fetchPom = mockPomSet(poms))

    result <- resolveDependencies("r", "root", "1", verbose = FALSE)

    expect_equal(result$artifactId, "guava")
})


# --- artifact download --------------------------------------------

# A dependency set in the shape resolveDependencies returns
depFrame <- function(...) {
    coords <- list(...)

    data.frame(
        groupId    = vapply(coords, `[[`, character(1), 1L),
        artifactId = vapply(coords, `[[`, character(1), 2L),
        version    = vapply(coords, `[[`, character(1), 3L),
        stringsAsFactors = FALSE
    )
}

# Stand in for a repository that serves every artifact asked of it
localRepo <- function(payload) {
    list(
        findArtifactUrl = function(path, ...) paste0("https://example.org/", path),
        downloadQuietly = function(url, destfile) {
            file.copy(payload, destfile, overwrite = TRUE)
            invisible(destfile)
        }
    )
}

test_that("downloadArtifacts is a no-op for an empty dependency set", {
    empty <- data.frame(
        groupId    = character(0),
        artifactId = character(0),
        version    = character(0)
    )

    expect_length(downloadArtifacts(empty, tempdir()), 0)
})

test_that("downloadQuietly writes the URL to the requested destination", {
    requested <- NULL

    local_mocked_bindings(
        download.file = function(url, destfile, ...) {
            requested <<- list(url = url, destfile = destfile)
            0L
        },
        .package = "utils"
    )

    destFile <- withr::local_tempfile(fileext = ".jar")

    expect_equal(downloadQuietly("https://example.org/x.jar", destFile), destFile)
    expect_equal(requested$url, "https://example.org/x.jar")
    expect_equal(requested$destfile, destFile)
})

test_that("downloadQuietly widens the timeout and then restores it", {
    duringDownload <- NULL

    local_mocked_bindings(
        download.file = function(...) {
            duringDownload <<- getOption("timeout")
            0L
        },
        .package = "utils"
    )

    withr::local_options(timeout = 13L)

    downloadQuietly("https://example.org/x.jar", withr::local_tempfile())

    expect_equal(duringDownload, MAVEN_FETCH_TIMEOUT)
    expect_equal(getOption("timeout"), 13L)
})

test_that("downloadArtifacts downloads every resolved artifact", {
    payload <- withr::local_tempfile(fileext = ".jar")
    writeLines("jar", payload)

    # A nested directory confirms the destination is created as needed
    destDir <- file.path(withr::local_tempdir(), "5.2.97")
    repo    <- localRepo(payload)

    local_mocked_bindings(
        findArtifactUrl   = repo$findArtifactUrl,
        downloadQuietly   = repo$downloadQuietly,
        fetchArtifactSha1 = function(...) NULL
    )

    result <- downloadArtifacts(depFrame(c("x", "a", "1"), c("x", "b", "2")), destDir)

    expect_equal(basename(result), c("a-1.jar", "b-2.jar"))
    expect_true(all(file.exists(result)))
})

test_that("downloadArtifacts keeps an artifact that is already cached", {
    destDir <- withr::local_tempdir()
    file.create(file.path(destDir, "a-1.jar"))

    local_mocked_bindings(
        findArtifactUrl = function(...) stop("no repository should be queried"),
        downloadQuietly = function(...) stop("nothing should be downloaded")
    )

    result <- downloadArtifacts(depFrame(c("x", "a", "1")), destDir)

    expect_equal(basename(result), "a-1.jar")
})

test_that("downloadArtifacts verifies an artifact against its published checksum", {
    payload <- withr::local_tempfile(fileext = ".jar")
    writeLines("jar", payload)
    sha1 <- digest::digest(payload, algo = "sha1", file = TRUE)

    destDir <- withr::local_tempdir()
    repo    <- localRepo(payload)

    local_mocked_bindings(
        findArtifactUrl   = repo$findArtifactUrl,
        downloadQuietly   = repo$downloadQuietly,
        fetchArtifactSha1 = function(...) sha1
    )

    result <- downloadArtifacts(depFrame(c("x", "a", "1")), destDir)

    expect_equal(basename(result), "a-1.jar")
})

test_that("downloadArtifacts discards an artifact that fails verification", {
    payload <- withr::local_tempfile(fileext = ".jar")
    writeLines("jar", payload)

    destDir <- withr::local_tempdir()
    repo    <- localRepo(payload)

    local_mocked_bindings(
        findArtifactUrl   = repo$findArtifactUrl,
        downloadQuietly   = repo$downloadQuietly,
        fetchArtifactSha1 = function(...) strrep("0", 40)
    )

    expect_warning(
        result <- downloadArtifacts(depFrame(c("x", "a", "1")), destDir),
        "Could not download 1 dependency"
    )

    expect_length(result, 0)
    expect_false(file.exists(file.path(destDir, "a-1.jar")))
})

test_that("downloadArtifacts warns about artifacts no repository serves", {
    local_mocked_bindings(findArtifactUrl = function(...) NULL)

    expect_warning(
        result <- downloadArtifacts(
            depFrame(c("x", "a", "1")),
            withr::local_tempdir()
        ),
        "Could not download 1 dependency"
    )

    expect_length(result, 0)
})

test_that("downloadArtifacts keeps going after a failed download", {
    payload <- withr::local_tempfile(fileext = ".jar")
    writeLines("jar", payload)

    destDir <- withr::local_tempdir()
    repo    <- localRepo(payload)

    local_mocked_bindings(
        findArtifactUrl   = repo$findArtifactUrl,
        downloadQuietly   = function(url, destfile) {
            if (grepl("a-1", url)) stop("connection reset")
            repo$downloadQuietly(url, destfile)
        },
        fetchArtifactSha1 = function(...) NULL
    )

    expect_warning(
        result <- downloadArtifacts(
            depFrame(c("x", "a", "1"), c("x", "b", "2")),
            destDir
        ),
        "x:a:1"
    )

    expect_equal(basename(result), "b-2.jar")
})

test_that("fetchArtifactSha1 extracts a checksum from sidecar text", {
    local_mocked_bindings(
        fetchUrlText = function(...) "9320966721A12741DA2A60F02FD3830639058D63  tassel.jar"
    )

    expect_equal(
        fetchArtifactSha1("https://example.org/x.jar"),
        "9320966721a12741da2a60f02fd3830639058d63"
    )
})

test_that("fetchArtifactSha1 returns NULL when no checksum is present", {
    local_mocked_bindings(fetchUrlText = function(...) "not a checksum")

    expect_null(fetchArtifactSha1("https://example.org/x.jar"))
})

test_that("fetchArtifactSha1 returns NULL when the sidecar is unreachable", {
    local_mocked_bindings(fetchUrlText = function(...) stop("404"))

    expect_null(fetchArtifactSha1("https://example.org/x.jar"))
})


# --- active version tracking --------------------------------------

test_that("getActiveTASSELVersion prefers the user option", {
    withr::local_options(rTASSEL.tassel.version = "1.2.3")

    expect_equal(getActiveTASSELVersion(), "1.2.3")
})

test_that("getActiveTASSELVersion falls back to the pinned version", {
    withr::local_options(rTASSEL.tassel.version = NULL)

    javaDir <- withr::local_tempdir()
    local_mocked_bindings(getTASSELJavaDir = function() javaDir)

    expect_equal(getActiveTASSELVersion(), TASSEL_MAVEN$VERSION)
})

test_that("getActiveTASSELVersion reads a recorded version", {
    withr::local_options(rTASSEL.tassel.version = NULL)

    javaDir <- withr::local_tempdir()
    dir.create(file.path(javaDir, "9.9.9"))

    local_mocked_bindings(getTASSELJavaDir = function() javaDir)

    setActiveTASSELVersion("9.9.9")

    expect_equal(getActiveTASSELVersion(), "9.9.9")
})

test_that("getActiveTASSELVersion ignores a recorded version with no cache dir", {
    withr::local_options(rTASSEL.tassel.version = NULL)

    javaDir <- withr::local_tempdir()
    local_mocked_bindings(getTASSELJavaDir = function() javaDir)

    setActiveTASSELVersion("9.9.9")

    expect_equal(getActiveTASSELVersion(), TASSEL_MAVEN$VERSION)
})


# --- JAR path resolution ------------------------------------------

test_that("getTASSELJarPath finds a fat JAR layout", {
    javaDir <- withr::local_tempdir()
    jarDir  <- file.path(javaDir, "5.2.96")
    dir.create(jarDir)
    file.create(file.path(jarDir, getTASSELJarName("5.2.96")))

    local_mocked_bindings(getTASSELJavaDir = function() javaDir)

    expect_equal(getTASSELJarPath("5.2.96"), jarDir)
})

test_that("getTASSELJarPath finds a thin JAR layout", {
    javaDir <- withr::local_tempdir()
    jarDir  <- file.path(javaDir, "5.2.97")
    dir.create(jarDir)
    file.create(file.path(jarDir, getTASSELThinJarName("5.2.97")))
    file.create(file.path(jarDir, "guava-22.0.jar"))

    local_mocked_bindings(getTASSELJavaDir = function() javaDir)

    expect_equal(getTASSELJarPath("5.2.97"), jarDir)
})

test_that("getTASSELJarPath returns NULL for a directory without TASSEL", {
    javaDir <- withr::local_tempdir()
    jarDir  <- file.path(javaDir, "5.2.97")
    dir.create(jarDir)
    file.create(file.path(jarDir, "guava-22.0.jar"))

    local_mocked_bindings(getTASSELJavaDir = function() javaDir)

    expect_null(getTASSELJarPath("5.2.97"))
})

test_that("getTASSELJarPath returns NULL when nothing is cached", {
    javaDir <- withr::local_tempdir()
    local_mocked_bindings(getTASSELJavaDir = function() javaDir)

    expect_null(getTASSELJarPath("5.2.96"))
})


# --- setupTASSEL guardrails ---------------------------------------

test_that("setupTASSEL aborts for a version with no usable artifacts", {
    javaDir <- withr::local_tempdir()

    local_mocked_bindings(
        getTASSELJavaDir        = function() javaDir,
        probeArtifactLayout     = function(...) "none",
        latestInstallableTASSEL = function(...) list(latest = "5.2.96", layout = "fat")
    )

    expect_error(
        setupTASSEL(version = "9.9.9"),
        "not available in a form"
    )
})

test_that("setupTASSEL reports an already cached version without downloading", {
    javaDir <- withr::local_tempdir()
    jarDir  <- file.path(javaDir, "5.2.96")
    dir.create(jarDir)
    file.create(file.path(jarDir, getTASSELJarName("5.2.96")))

    downloaded <- FALSE

    local_mocked_bindings(
        getTASSELJavaDir  = function() javaDir,
        downloadTASSELJar = function(...) {
            downloaded <<- TRUE
            invisible(NULL)
        }
    )

    expect_message(setupTASSEL(version = "5.2.96"), "already cached")
    expect_false(downloaded)
})
