# === Tests for GitHub release installation =========================


# --- Fixtures ------------------------------------------------------

# A trimmed copy of the releases API response, covering a nightly
# prerelease, an older nightly, a tagged release, and the installer-only
# release that carries no standalone archive.
mockReleases <- function() {
    list(
        list(
            tag_name     = "dev-20260801",
            prerelease   = TRUE,
            draft        = FALSE,
            published_at = "2026-08-01T07:45:34Z",
            assets = list(
                list(
                    name                 = "tassel-5-standalone-v5.2.98-dev.20260801.tar.gz",
                    size                 = 68859968,
                    digest               = "sha256:aaaa",
                    browser_download_url = "https://example.invalid/dev-20260801/tassel-5-standalone-v5.2.98-dev.20260801.tar.gz"
                ),
                list(
                    name                 = "tassel-5-standalone-v5.2.98-dev.20260801.zip",
                    size                 = 68868796,
                    digest               = "sha256:bbbb",
                    browser_download_url = "https://example.invalid/dev-20260801/tassel-5-standalone-v5.2.98-dev.20260801.zip"
                )
            )
        ),
        list(
            tag_name     = "dev-20260730",
            prerelease   = TRUE,
            draft        = FALSE,
            published_at = "2026-07-30T07:50:25Z",
            assets = list(
                list(
                    name                 = "tassel-5-standalone-v5.2.97-dev.20260730.tar.gz",
                    size                 = 66537136,
                    digest               = "sha256:cccc",
                    browser_download_url = "https://example.invalid/dev-20260730/tassel-5-standalone-v5.2.97-dev.20260730.tar.gz"
                )
            )
        ),
        list(
            tag_name     = "v5.2.97",
            prerelease   = FALSE,
            draft        = FALSE,
            published_at = "2026-07-19T17:55:22Z",
            assets = list(
                list(
                    name                 = "tassel-5-standalone-v5.2.97.tar.gz",
                    size                 = 67903842,
                    digest               = "sha256:dddd",
                    browser_download_url = "https://example.invalid/v5.2.97/tassel-5-standalone-v5.2.97.tar.gz"
                )
            )
        ),
        list(
            tag_name     = "main",
            prerelease   = TRUE,
            draft        = FALSE,
            published_at = "2026-07-19T16:52:46Z",
            assets = list(
                list(
                    name                 = "TASSEL.5.Installer-win-x64.exe",
                    size                 = 10506912,
                    browser_download_url = "https://example.invalid/main/TASSEL.5.Installer-win-x64.exe"
                )
            )
        )
    )
}

# Build a standalone archive with the layout the real ones use: the main
# JAR beside a 'lib' directory, both at the archive root.
makeStandaloneArchive <- function(
    version = "5.2.98-dev.20260801",
    libJars = c("guava-22.0.jar", "colt-1.2.0.jar"),
    mainJar = TASSEL_GITHUB$MAIN_JAR
) {
    stage   <- withr::local_tempdir(.local_envir = parent.frame())
    archive <- withr::local_tempfile(
        fileext      = ".tar.gz",
        .local_envir = parent.frame()
    )

    dir.create(file.path(stage, TASSEL_GITHUB$LIB_DIR), showWarnings = FALSE)

    if (!is.null(mainJar)) {
        writeLines("main jar", file.path(stage, mainJar))
    }

    for (jar in libJars) {
        writeLines(jar, file.path(stage, TASSEL_GITHUB$LIB_DIR, jar))
    }

    writeLines("#!/bin/sh", file.path(stage, "start_tassel.pl"))

    withr::with_dir(
        stage,
        utils::tar(archive, files = list.files("."), compression = "gzip")
    )

    list(
        path    = archive,
        version = version,
        release = list(
            tag         = "dev-20260801",
            version     = version,
            publishedAt = "2026-08-01T07:45:34Z",
            prerelease  = TRUE,
            assetName   = basename(archive),
            assetUrl    = paste0("file://", archive),
            sha256      = digest::digest(archive, algo = "sha256", file = TRUE),
            size        = file.size(archive)
        )
    )
}


# --- API request headers -------------------------------------------

test_that("githubApiHeaders always requests the versioned API media type", {
    withr::local_envvar(GITHUB_PAT = "", GITHUB_TOKEN = "")

    headers <- githubApiHeaders()

    expect_equal(headers[["Accept"]], "application/vnd.github+json")
    expect_false("Authorization" %in% names(headers))
})

test_that("githubApiHeaders uses a token from the environment", {
    withr::local_envvar(GITHUB_PAT = "secret", GITHUB_TOKEN = "")

    expect_equal(githubApiHeaders()[["Authorization"]], "Bearer secret")
})

test_that("githubApiHeaders falls back to GITHUB_TOKEN", {
    withr::local_envvar(GITHUB_PAT = "", GITHUB_TOKEN = "other")

    expect_equal(githubApiHeaders()[["Authorization"]], "Bearer other")
})


# --- asset parsing -------------------------------------------------

test_that("standaloneAssetVersion recovers the version from an asset name", {
    expect_equal(
        standaloneAssetVersion("tassel-5-standalone-v5.2.98-dev.20260801.tar.gz"),
        "5.2.98-dev.20260801"
    )
    expect_equal(
        standaloneAssetVersion("tassel-5-standalone-v5.2.97.tar.gz"),
        "5.2.97"
    )
})

test_that("standaloneAssetVersion rejects names it cannot parse", {
    expect_true(is.na(standaloneAssetVersion("tassel-5-standalone-v5.2.97.zip")))
    expect_true(is.na(standaloneAssetVersion("TASSEL.5.Installer-win-x64.exe")))
    expect_true(is.na(standaloneAssetVersion("tassel-5-standalone-v.tar.gz")))
})

test_that("standaloneAsset picks the tarball over the other assets", {
    asset <- standaloneAsset(mockReleases()[[1]])

    expect_equal(asset$name, "tassel-5-standalone-v5.2.98-dev.20260801.tar.gz")
})

test_that("standaloneAsset returns NULL when no archive is published", {
    expect_null(standaloneAsset(mockReleases()[[4]]))
})


# --- release records -----------------------------------------------

test_that("asReleaseRecord extracts the fields the installer needs", {
    record <- asReleaseRecord(mockReleases()[[1]])

    expect_equal(record$tag, "dev-20260801")
    expect_equal(record$version, "5.2.98-dev.20260801")
    expect_true(record$prerelease)
    expect_equal(record$sha256, "aaaa")
    expect_equal(record$size, 68859968)
})

test_that("asReleaseRecord skips releases without an installable archive", {
    expect_null(asReleaseRecord(mockReleases()[[4]]))
})

test_that("asReleaseRecord skips drafts", {
    release <- mockReleases()[[1]]
    release$draft <- TRUE

    expect_null(asReleaseRecord(release))
})

test_that("asReleaseRecord reports no checksum when none is published", {
    release <- mockReleases()[[1]]
    release$assets[[1]]$digest <- NULL

    expect_null(asReleaseRecord(release)$sha256)
})

test_that("asReleaseRecord ignores a digest in an unexpected algorithm", {
    release <- mockReleases()[[1]]
    release$assets[[1]]$digest <- "md5:aaaa"

    expect_null(asReleaseRecord(release)$sha256)
})


# --- release listing -----------------------------------------------

test_that("githubReleaseCandidates keeps only installable releases, newest first", {
    local_mocked_bindings(fetchGitHubJson = function(...) mockReleases())

    releases <- githubReleaseCandidates()

    expect_length(releases, 3)
    expect_equal(
        vapply(releases, function(x) x$version, character(1)),
        c("5.2.98-dev.20260801", "5.2.97-dev.20260730", "5.2.97")
    )
})

test_that("githubReleaseCandidates returns an empty list when nothing matches", {
    local_mocked_bindings(fetchGitHubJson = function(...) mockReleases()[4])

    expect_length(githubReleaseCandidates(), 0)
})

test_that("githubReleaseCandidates requests the TASSEL repository", {
    requested <- NULL

    local_mocked_bindings(fetchGitHubJson = function(path, ...) {
        requested <<- path
        list()
    })

    githubReleaseCandidates()

    expect_match(requested, TASSEL_GITHUB$REPO, fixed = TRUE)
    expect_match(requested, "releases", fixed = TRUE)
})


# --- spec resolution -----------------------------------------------

test_that("resolveGitHubRelease finds the newest nightly", {
    local_mocked_bindings(fetchGitHubJson = function(...) mockReleases())

    expect_equal(resolveGitHubRelease("nightly")$version, "5.2.98-dev.20260801")
    expect_equal(resolveGitHubRelease("dev")$version, "5.2.98-dev.20260801")
})

test_that("resolveGitHubRelease finds the newest tagged release", {
    local_mocked_bindings(fetchGitHubJson = function(...) mockReleases())

    expect_equal(resolveGitHubRelease("latest")$version, "5.2.97")
})

test_that("resolveGitHubRelease matches an exact tag", {
    local_mocked_bindings(fetchGitHubJson = function(...) mockReleases())

    expect_equal(resolveGitHubRelease("dev-20260730")$version, "5.2.97-dev.20260730")
})

test_that("resolveGitHubRelease matches a version, with or without a leading v", {
    local_mocked_bindings(fetchGitHubJson = function(...) mockReleases())

    expect_equal(resolveGitHubRelease("5.2.98-dev.20260801")$tag, "dev-20260801")
    expect_equal(resolveGitHubRelease("v5.2.97")$tag, "v5.2.97")
})

test_that("resolveGitHubRelease reports what is available when nothing matches", {
    local_mocked_bindings(fetchGitHubJson = function(...) mockReleases())

    expect_error(resolveGitHubRelease("9.9.9"), "9.9.9")
    expect_error(resolveGitHubRelease("9.9.9"), "5.2.98-dev.20260801")
})

test_that("resolveGitHubRelease errors when no archives are published", {
    local_mocked_bindings(fetchGitHubJson = function(...) list())

    expect_error(resolveGitHubRelease("nightly"), "No installable TASSEL archives")
})


# --- nightly version grammar ---------------------------------------

test_that("isNightlyVersion recognises the nightly grammar", {
    expect_true(isNightlyVersion("5.2.98-dev.20260801"))
    expect_true(isNightlyVersion("5.2-dev.20260801"))

    expect_false(isNightlyVersion("5.2.97"))
    expect_false(isNightlyVersion("5.2.98-dev"))
    expect_false(isNightlyVersion("5.2.98-dev.2026080"))
    expect_false(isNightlyVersion("5.2.98-snapshot.20260801"))
})

test_that("isNightlyVersion tolerates values that are not versions", {
    expect_false(isNightlyVersion(NULL))
    expect_false(isNightlyVersion(NA_character_))
    expect_false(isNightlyVersion(c("5.2.98-dev.20260801", "5.2.97")))
    expect_false(isNightlyVersion(5.2))
})

test_that("parseNightlyVersion splits a nightly into version and build date", {
    parsed <- parseNightlyVersion("5.2.98-dev.20260801")

    expect_equal(parsed$base, "5.2.98")
    expect_equal(parsed$date, "20260801")
    expect_null(parseNightlyVersion("5.2.97"))
})

test_that("isNewerNightly compares build dates within a release version", {
    expect_true(isNewerNightly("5.2.98-dev.20260802", "5.2.98-dev.20260801"))
    expect_false(isNewerNightly("5.2.98-dev.20260801", "5.2.98-dev.20260801"))
    expect_false(isNewerNightly("5.2.98-dev.20260731", "5.2.98-dev.20260801"))
})

test_that("isNewerNightly prefers the newer release version over the newer date", {
    expect_true(isNewerNightly("5.2.98-dev.20260701", "5.2.97-dev.20260801"))
    expect_false(isNewerNightly("5.2.97-dev.20260801", "5.2.98-dev.20260701"))
})

test_that("isNewerNightly refuses to compare non-nightly versions", {
    expect_false(isNewerNightly("5.2.98", "5.2.97-dev.20260801"))
    expect_false(isNewerNightly("5.2.98-dev.20260801", "5.2.97"))
})


# --- latest nightly lookup -----------------------------------------

test_that("latestNightlyTASSEL summarises the newest nightly", {
    local_mocked_bindings(fetchGitHubJson = function(...) mockReleases())

    result <- latestNightlyTASSEL()

    expect_equal(result$latest, "5.2.98-dev.20260801")
    expect_equal(result$tag, "dev-20260801")
    expect_equal(result$channel, "nightly")
})

test_that("latestNightlyTASSEL returns NULL rather than failing when offline", {
    local_mocked_bindings(fetchGitHubJson = function(...) stop("no network"))

    expect_null(latestNightlyTASSEL())
})


# --- checksum verification -----------------------------------------

test_that("verifyGitHubChecksum accepts a matching archive", {
    file <- withr::local_tempfile()
    writeLines("payload", file)

    release <- list(sha256 = digest::digest(file, algo = "sha256", file = TRUE))

    expect_true(suppressMessages(verifyGitHubChecksum(file, release)))
    expect_true(file.exists(file))
})

test_that("verifyGitHubChecksum discards an archive that fails to verify", {
    file <- withr::local_tempfile()
    writeLines("payload", file)

    expect_error(
        suppressMessages(verifyGitHubChecksum(file, list(sha256 = "deadbeef"))),
        "checksum verification failed"
    )
    expect_false(file.exists(file))
})

test_that("verifyGitHubChecksum skips verification when no checksum is published", {
    file <- withr::local_tempfile()
    writeLines("payload", file)

    expect_message(
        expect_true(verifyGitHubChecksum(file, list(sha256 = NULL))),
        "skipping verification"
    )
})


# --- archive layout ------------------------------------------------

test_that("standaloneRoot finds a flat archive root", {
    dir <- withr::local_tempdir()
    writeLines("main jar", file.path(dir, TASSEL_GITHUB$MAIN_JAR))

    expect_equal(standaloneRoot(dir), dir)
})

test_that("standaloneRoot finds a root nested one level down", {
    dir <- withr::local_tempdir()
    dir.create(file.path(dir, "tassel-5-standalone"))
    writeLines(
        "main jar",
        file.path(dir, "tassel-5-standalone", TASSEL_GITHUB$MAIN_JAR)
    )

    expect_equal(
        standaloneRoot(dir),
        file.path(dir, "tassel-5-standalone")
    )
})

test_that("standaloneRoot returns NULL when the main JAR is absent", {
    dir <- withr::local_tempdir()
    writeLines("not a jar", file.path(dir, "start_tassel.pl"))

    expect_null(standaloneRoot(dir))
})


# --- archive installation ------------------------------------------

test_that("installGitHubStandalone flattens the archive into one directory", {
    archive <- makeStandaloneArchive()
    jarDir  <- withr::local_tempdir()

    installed <- suppressMessages(
        installGitHubStandalone(archive$release, jarDir)
    )

    expect_setequal(
        basename(installed),
        c(TASSEL_GITHUB$MAIN_JAR, "guava-22.0.jar", "colt-1.2.0.jar")
    )
    expect_true(all(file.exists(installed)))

    # Only JARs belong on the classpath, so the launcher scripts and the
    # lib directory itself must not come along.
    expect_setequal(
        list.files(jarDir),
        c(TASSEL_GITHUB$MAIN_JAR, "guava-22.0.jar", "colt-1.2.0.jar")
    )
})

test_that("installGitHubStandalone leaves a directory getTASSELJarPath accepts", {
    archive <- makeStandaloneArchive()
    javaDir <- withr::local_tempdir()

    local_mocked_bindings(getTASSELJavaDir = function() javaDir)

    jarDir <- getTASSELCacheDir(archive$version)
    dir.create(jarDir, recursive = TRUE)

    suppressMessages(installGitHubStandalone(archive$release, jarDir))

    expect_equal(getTASSELJarPath(archive$version), jarDir)
})

test_that("installGitHubStandalone rejects an archive without the main JAR", {
    archive <- makeStandaloneArchive(mainJar = NULL)
    jarDir  <- withr::local_tempdir()

    expect_error(
        suppressMessages(installGitHubStandalone(archive$release, jarDir)),
        TASSEL_GITHUB$MAIN_JAR,
        fixed = TRUE
    )
})

test_that("installGitHubStandalone aborts on a checksum mismatch", {
    archive <- makeStandaloneArchive()
    jarDir  <- withr::local_tempdir()

    release <- archive$release
    release$sha256 <- "deadbeef"

    expect_error(
        suppressMessages(installGitHubStandalone(release, jarDir)),
        "checksum verification failed"
    )
    expect_length(list.files(jarDir), 0)
})

test_that("installGitHubStandalone reports a failed download", {
    release <- list(
        version   = "5.2.98-dev.20260801",
        assetName = "missing.tar.gz",
        assetUrl  = paste0("file://", file.path(tempdir(), "does-not-exist.tar.gz")),
        sha256    = NULL,
        size      = NA_real_
    )

    expect_error(
        suppressMessages(suppressWarnings(
            installGitHubStandalone(release, withr::local_tempdir())
        )),
        "Failed to download"
    )
})


# --- setupTASSEL dispatch ------------------------------------------

test_that("setupTASSEL installs a nightly build from GitHub", {
    archive <- makeStandaloneArchive()
    javaDir <- withr::local_tempdir()

    local_mocked_bindings(
        getTASSELJavaDir     = function() javaDir,
        resolveGitHubRelease = function(spec, ...) {
            expect_equal(spec, "nightly")
            archive$release
        }
    )

    jarDir <- suppressMessages(setupTASSEL(source = "github"))

    expect_equal(jarDir, file.path(javaDir, archive$version))
    expect_true(file.exists(file.path(jarDir, TASSEL_GITHUB$MAIN_JAR)))
    expect_equal(getActiveTASSELVersion(), archive$version)
    expect_equal(readInstallSource(archive$version), "github")
})

test_that("setupTASSEL passes an explicit spec through to GitHub", {
    archive <- makeStandaloneArchive()
    javaDir <- withr::local_tempdir()
    seen    <- NULL

    local_mocked_bindings(
        getTASSELJavaDir     = function() javaDir,
        resolveGitHubRelease = function(spec, ...) {
            seen <<- spec
            archive$release
        }
    )

    suppressMessages(setupTASSEL(version = "dev-20260801", source = "github"))

    expect_equal(seen, "dev-20260801")
})

test_that("setupTASSEL does not reach Maven when installing from GitHub", {
    archive <- makeStandaloneArchive()
    javaDir <- withr::local_tempdir()

    local_mocked_bindings(
        getTASSELJavaDir     = function() javaDir,
        resolveGitHubRelease = function(...) archive$release,
        probeArtifactLayout  = function(...) stop("Maven must not be probed")
    )

    expect_no_error(suppressMessages(setupTASSEL(source = "github")))
})

test_that("setupTASSEL reports a cached nightly without downloading again", {
    archive <- makeStandaloneArchive()
    javaDir <- withr::local_tempdir()

    local_mocked_bindings(
        getTASSELJavaDir     = function() javaDir,
        resolveGitHubRelease = function(...) archive$release
    )

    suppressMessages(setupTASSEL(source = "github"))

    expect_message(
        setupTASSEL(source = "github"),
        "already cached"
    )
})

test_that("setupTASSEL rejects an unknown source", {
    expect_error(setupTASSEL(source = "bitbucket"))
})
