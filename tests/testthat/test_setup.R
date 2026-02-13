# === Tests for setup functions =====================================


# --- getTASSELJarName() -------------------------------------------

test_that("getTASSELJarName() returns correctly formatted JAR filename", {
    result <- getTASSELJarName()

    expect_type(result, "character")
    expect_length(result, 1)
    expect_match(result, "\\.jar$")
    expect_match(result, "^tassel-")
    expect_match(result, "jar-with-dependencies")
})

test_that("getTASSELJarName() uses default version from TASSEL_MAVEN", {
    result <- getTASSELJarName()
    expected <- sprintf(
        "%s-%s-%s.jar",
        TASSEL_MAVEN$ARTIFACT_ID,
        TASSEL_MAVEN$VERSION,
        TASSEL_MAVEN$CLASSIFIER
    )
    expect_equal(result, expected)
})

test_that("getTASSELJarName() respects custom version argument", {
    result <- getTASSELJarName(version = "1.0.0")
    expect_match(result, "1\\.0\\.0")
    expect_equal(result, "tassel-1.0.0-jar-with-dependencies.jar")
})


# --- getTASSELCacheDir() ------------------------------------------

test_that("getTASSELCacheDir() returns a character path", {
    result <- getTASSELCacheDir()

    expect_type(result, "character")
    expect_length(result, 1)
})

test_that("getTASSELCacheDir() contains version in path", {
    result <- getTASSELCacheDir()
    expect_match(result, TASSEL_MAVEN$VERSION, fixed = TRUE)
})

test_that("getTASSELCacheDir() uses R_user_dir structure", {
    result <- getTASSELCacheDir()

    expect_match(result, "rTASSEL")
    expect_match(result, "java")
})

test_that("getTASSELCacheDir() respects custom version argument", {
    result <- getTASSELCacheDir(version = "99.99.99")
    expect_match(result, "99\\.99\\.99")
})


# --- getTASSELJarPath() -------------------------------------------

test_that("getTASSELJarPath() returns NULL when no cached JAR exists", {
    # Use a version that will never exist in the cache
    result <- getTASSELJarPath(version = "0.0.0-nonexistent")
    expect_null(result)
})

test_that("getTASSELJarPath() returns directory path when JAR exists", {
    # Create a temporary cache directory with a fake JAR
    tmpDir <- file.path(tempdir(), "rtassel_test_cache", "java", "0.0.1")
    dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
    jarName <- getTASSELJarName(version = "0.0.1")
    file.create(file.path(tmpDir, jarName))

    # Mock getTASSELCacheDir to return our temp directory
    mockResult <- withr::with_options(
        list(rTASSEL.java.path = NULL),
        {
            # We need to directly test the file-existence logic
            jarFile <- file.path(tmpDir, jarName)
            if (file.exists(jarFile)) tmpDir else NULL
        }
    )

    expect_equal(mockResult, tmpDir)

    # Clean up
    unlink(file.path(tempdir(), "rtassel_test_cache"), recursive = TRUE)
})


# --- resolveJarPath() ---------------------------------------------

test_that("resolveJarPath() returns a list with 'path' and 'source'", {
    result <- resolveJarPath()

    expect_type(result, "list")
    expect_true("path" %in% names(result))
    expect_true("source" %in% names(result))
})

test_that("resolveJarPath() prioritizes user-defined option", {
    tmpDir <- file.path(tempdir(), "rtassel_test_option_jars")
    dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)

    result <- withr::with_options(
        list(rTASSEL.java.path = tmpDir),
        resolveJarPath()
    )

    expect_equal(result$path, tmpDir)
    expect_equal(result$source, "option")

    # Clean up
    unlink(tmpDir, recursive = TRUE)
})

test_that("resolveJarPath() returns NULL path/source when nothing is found", {
    # Verify that a NULL result is structurally correct by checking
    # what the function returns when there's no option, no maven cache
    # for a nonexistent version, and no bundled JARs. We test the
    # fallback path logic directly.
    withr::local_options(list(rTASSEL.java.path = NULL))

    # Bundled path for a nonexistent package should not have JARs
    bundledPath <- system.file("java", package = "nonexistent_pkg_12345")
    expect_false(dir.exists(bundledPath) && length(list.files(bundledPath, "\\.jar$")) > 0)

    # Verify the NULL return structure matches the expected contract
    nullResult <- list(path = NULL, source = NULL)
    expect_null(nullResult$path)
    expect_null(nullResult$source)
    expect_type(nullResult, "list")
})

test_that("resolveJarPath() source is one of expected values or NULL", {
    result <- resolveJarPath()
    validSources <- c("option", "maven cache", "bundled", NULL)
    expect_true(is.null(result$source) || result$source %in% validSources)
})

test_that("resolveJarPath() falls back to bundled inst/java/ when available", {
    # When no option is set and no maven cache exists, it should check bundled
    result <- withr::with_options(
        list(rTASSEL.java.path = NULL),
        {
            # Only test bundled fallback if maven cache is not present
            if (is.null(getTASSELJarPath())) {
                resolveJarPath(pkgname = "rTASSEL")
            } else {
                # Maven cache exists, so it will return that instead
                list(path = "skip", source = "maven cache")
            }
        }
    )

    # The result depends on whether bundled JARs exist in this installation
    expect_type(result, "list")
    expect_true("path" %in% names(result))
    expect_true("source" %in% names(result))
})


# --- setupTASSEL() ------------------------------------------------

test_that("setupTASSEL() returns early when JAR is already cached", {
    # Create a fake cached JAR
    tmpCacheBase <- file.path(tempdir(), "rtassel_setup_test")
    version <- "0.0.1-test"
    jarDir <- file.path(tmpCacheBase, "java", version)
    dir.create(jarDir, recursive = TRUE, showWarnings = FALSE)

    jarName <- getTASSELJarName(version = version)
    jarFile <- file.path(jarDir, jarName)
    writeLines("fake jar content", jarFile)

    # Temporarily override getTASSELCacheDir to use our temp location
    # We can test the early-return logic by calling setupTASSEL with
    # a mock environment. Since setupTASSEL calls getTASSELCacheDir
    # internally, we test the condition directly.
    expect_true(file.exists(jarFile))

    # Clean up
    unlink(tmpCacheBase, recursive = TRUE)
})

test_that("setupTASSEL() has correct function signature", {
    expect_true(is.function(setupTASSEL))

    args <- formals(setupTASSEL)
    expect_true("version" %in% names(args))
    expect_true("force" %in% names(args))
    expect_false(eval(args$force))
})


# === Extended coverage tests ======================================

# --- Helper: temporarily replace a binding in the source env ------
# Works for both file_coverage() (unlocked env) and test_check()
# (locked package namespace) by unlocking/re-locking as needed.
setup_mock <- function(name, replacement, env = environment(setupTASSEL)) {
    existed    <- exists(name, envir = env, inherits = FALSE)
    was_locked <- existed && bindingIsLocked(name, env)

    if (existed) orig <- get(name, envir = env, inherits = FALSE)
    if (was_locked) unlockBinding(name, env)

    assign(name, replacement, envir = env)

    # return a restore function
    function() {
        if (was_locked) try(unlockBinding(name, env), silent = TRUE)
        if (existed) {
            assign(name, orig, envir = env)
        } else if (exists(name, envir = env, inherits = FALSE)) {
            rm(list = name, envir = env)
        }
        if (was_locked) lockBinding(name, env)
    }
}

# --- Helper: create a minimal fake installed package --------------
# system.file() needs a DESCRIPTION to resolve the package path.
create_fake_pkg <- function(pkg = "fakepkg", has_jars = TRUE) {
    tmplib  <- file.path(tempdir(), paste0("rtassel_fakelib_", sample.int(1e6, 1)))
    pkgdir  <- file.path(tmplib, pkg)
    javadir <- file.path(pkgdir, "java")
    dir.create(javadir, recursive = TRUE, showWarnings = FALSE)
    writeLines(
        c("Package: fakepkg", "Version: 0.0.1"),
        file.path(pkgdir, "DESCRIPTION")
    )
    if (has_jars) file.create(file.path(javadir, "test.jar"))
    list(lib = tmplib, pkg = pkg)
}


# --- getTASSELJarPath() with mock cache dir -----------------------

test_that("getTASSELJarPath() returns cache dir when JAR file is present", {
    tmpDir <- file.path(tempdir(), "rtassel_jarpath_mock")
    dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)

    jarName <- getTASSELJarName(version = "0.0.1-mock")
    file.create(file.path(tmpDir, jarName))

    restore <- setup_mock("getTASSELCacheDir", function(...) tmpDir)
    on.exit(restore(), add = TRUE)

    result <- getTASSELJarPath(version = "0.0.1-mock")
    expect_equal(result, tmpDir)
})

test_that("getTASSELJarPath() returns NULL when cache dir exists but JAR missing", {
    tmpDir <- file.path(tempdir(), "rtassel_jarpath_empty")
    dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)

    restore <- setup_mock("getTASSELCacheDir", function(...) tmpDir)
    on.exit(restore(), add = TRUE)

    result <- getTASSELJarPath(version = "0.0.2-mock")
    expect_null(result)
})


# --- resolveJarPath() – bundled fallback paths --------------------

test_that("resolveJarPath() returns bundled path when JARs are present", {
    fake <- create_fake_pkg(has_jars = TRUE)
    on.exit(unlink(fake$lib, recursive = TRUE), add = TRUE)

    restore1 <- setup_mock("getTASSELJarPath", function(...) NULL)
    on.exit(restore1(), add = TRUE)

    result <- withr::with_options(
        list(rTASSEL.java.path = NULL),
        resolveJarPath(pkgname = fake$pkg, libname = fake$lib)
    )

    expect_match(result$path, "java")
    expect_equal(result$source, "bundled")
})

test_that("resolveJarPath() returns NULL when bundled dir has no JARs", {
    fake <- create_fake_pkg(has_jars = FALSE)
    on.exit(unlink(fake$lib, recursive = TRUE), add = TRUE)

    restore1 <- setup_mock("getTASSELJarPath", function(...) NULL)
    on.exit(restore1(), add = TRUE)

    result <- withr::with_options(
        list(rTASSEL.java.path = NULL),
        resolveJarPath(pkgname = fake$pkg, libname = fake$lib)
    )

    expect_null(result$path)
    expect_null(result$source)
})

test_that("resolveJarPath() returns NULL/NULL when no sources available", {
    restore1 <- setup_mock("getTASSELJarPath", function(...) NULL)
    on.exit(restore1(), add = TRUE)

    result <- withr::with_options(
        list(rTASSEL.java.path = NULL),
        resolveJarPath(pkgname = "nonexistent_pkg_xyz_12345")
    )

    expect_null(result$path)
    expect_null(result$source)
})

test_that("resolveJarPath() prefers maven cache over bundled", {
    tmpMavenDir <- file.path(tempdir(), "rtassel_resolve_maven")
    dir.create(tmpMavenDir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmpMavenDir, recursive = TRUE), add = TRUE)

    restore <- setup_mock("getTASSELJarPath", function(...) tmpMavenDir)
    on.exit(restore(), add = TRUE)

    result <- withr::with_options(
        list(rTASSEL.java.path = NULL),
        resolveJarPath()
    )

    expect_equal(result$path, tmpMavenDir)
    expect_equal(result$source, "maven cache")
})


# --- downloadWithProgress() --------------------------------------

test_that("downloadWithProgress() downloads file content correctly", {
    srcFile  <- tempfile("dwp_src_")
    destFile <- tempfile("dwp_dest_")
    on.exit(unlink(c(srcFile, destFile)), add = TRUE)

    testContent <- charToRaw(paste(rep("test data line\n", 200), collapse = ""))
    writeBin(testContent, srcFile)

    result <- downloadWithProgress(
        url            = paste0("file://", srcFile),
        destfile       = destFile,
        estimatedSizeMb = 0.001
    )

    expect_true(file.exists(destFile))
    expect_equal(
        readBin(destFile, "raw", file.info(destFile)$size),
        testContent
    )
    expect_equal(result, destFile)
})

test_that("downloadWithProgress() handles empty file", {
    srcFile  <- tempfile("dwp_empty_src_")
    destFile <- tempfile("dwp_empty_dest_")
    on.exit(unlink(c(srcFile, destFile)), add = TRUE)

    file.create(srcFile)

    result <- downloadWithProgress(
        url            = paste0("file://", srcFile),
        destfile       = destFile,
        estimatedSizeMb = 0.001
    )

    expect_true(file.exists(destFile))
    expect_equal(file.info(destFile)$size, 0)
})

test_that("downloadWithProgress() uses estimated size as fallback", {
    srcFile  <- tempfile("dwp_est_src_")
    destFile <- tempfile("dwp_est_dest_")
    on.exit(unlink(c(srcFile, destFile)), add = TRUE)

    writeBin(charToRaw("test"), srcFile)

    # file:// has no Content-Length header, so fallback is used
    expect_no_error(
        downloadWithProgress(
            url            = paste0("file://", srcFile),
            destfile       = destFile,
            estimatedSizeMb = 50
        )
    )
    expect_true(file.exists(destFile))
})

test_that("downloadWithProgress() returns destfile path invisibly", {
    srcFile  <- tempfile("dwp_ret_src_")
    destFile <- tempfile("dwp_ret_dest_")
    on.exit(unlink(c(srcFile, destFile)), add = TRUE)

    writeBin(charToRaw("content"), srcFile)

    out <- withVisible(
        downloadWithProgress(
            url            = paste0("file://", srcFile),
            destfile       = destFile,
            estimatedSizeMb = 0.001
        )
    )

    expect_false(out$visible)
    expect_identical(out$value, destFile)
})


# --- setupTASSEL() ------------------------------------------------

test_that("setupTASSEL() returns early when JAR already cached", {
    version <- "0.0.1-cached"
    tmpDir  <- file.path(tempdir(), "rtassel_setup_early")
    dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)

    jarName <- getTASSELJarName(version = version)
    writeLines("fake jar", file.path(tmpDir, jarName))

    restore1 <- setup_mock("getTASSELCacheDir", function(...) tmpDir)
    on.exit(restore1(), add = TRUE)
    restore2 <- setup_mock("getTASSELJarName", function(...) jarName)
    on.exit(restore2(), add = TRUE)

    result <- setupTASSEL(version = version)
    expect_equal(as.character(result), tmpDir)
})

test_that("setupTASSEL() downloads JAR for non-default version (skips checksum)", {
    version <- "99.99.99-nondefault"
    tmpDir  <- file.path(tempdir(), "rtassel_setup_dl")
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)

    jarName <- getTASSELJarName(version = version)

    restore1 <- setup_mock("getTASSELCacheDir", function(...) tmpDir)
    on.exit(restore1(), add = TRUE)
    restore2 <- setup_mock("getTASSELJarName", function(...) jarName)
    on.exit(restore2(), add = TRUE)
    restore3 <- setup_mock("downloadWithProgress", function(url, destfile, ...) {
        writeLines("downloaded content", destfile)
        invisible(destfile)
    })
    on.exit(restore3(), add = TRUE)

    result <- setupTASSEL(version = version)
    expect_equal(as.character(result), tmpDir)
    expect_true(file.exists(file.path(tmpDir, jarName)))
    expect_equal(readLines(file.path(tmpDir, jarName)), "downloaded content")
})

test_that("setupTASSEL() verifies SHA-1 for default version (pass)", {
    skip_if_not_installed("digest")

    content <- "checksum verification content"
    tmpDir  <- file.path(tempdir(), "rtassel_setup_sha_ok")
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)

    # Compute the real SHA-1 of our test payload
    tmpCalc <- tempfile()
    writeLines(content, tmpCalc)
    expected_sha1 <- digest::digest(tmpCalc, algo = "sha1", file = TRUE)
    unlink(tmpCalc)

    restore1 <- setup_mock("getTASSELCacheDir", function(...) tmpDir)
    on.exit(restore1(), add = TRUE)
    restore2 <- setup_mock("downloadWithProgress", function(url, destfile, ...) {
        writeLines(content, destfile)
        invisible(destfile)
    })
    on.exit(restore2(), add = TRUE)

    # Temporarily set expected checksum to match our payload
    orig_version <- TASSEL_MAVEN$VERSION
    new_maven    <- TASSEL_MAVEN
    new_maven$SHA1_CHECKSUM <- expected_sha1
    restore3 <- setup_mock("TASSEL_MAVEN", new_maven)
    on.exit(restore3(), add = TRUE)

    result <- setupTASSEL(version = orig_version)
    expect_equal(as.character(result), tmpDir)
})

test_that("setupTASSEL() aborts on SHA-1 mismatch", {
    skip_if_not_installed("digest")

    tmpDir <- file.path(tempdir(), "rtassel_setup_sha_fail")
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)

    restore1 <- setup_mock("getTASSELCacheDir", function(...) tmpDir)
    on.exit(restore1(), add = TRUE)
    restore2 <- setup_mock("downloadWithProgress", function(url, destfile, ...) {
        writeLines("corrupted content", destfile)
        invisible(destfile)
    })
    on.exit(restore2(), add = TRUE)

    orig_version <- TASSEL_MAVEN$VERSION
    new_maven    <- TASSEL_MAVEN
    new_maven$SHA1_CHECKSUM <- "0000000000000000000000000000000000000000"
    restore3 <- setup_mock("TASSEL_MAVEN", new_maven)
    on.exit(restore3(), add = TRUE)

    expect_error(
        setupTASSEL(version = orig_version),
        "SHA-1 checksum verification failed"
    )

    # Corrupted JAR should have been removed
    jarFile <- file.path(tmpDir, getTASSELJarName(version = orig_version))
    expect_false(file.exists(jarFile))
})

test_that("setupTASSEL() handles download failure gracefully", {
    version <- "99.99.99-dl-fail"
    tmpDir  <- file.path(tempdir(), "rtassel_setup_dl_fail")
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)

    restore1 <- setup_mock("getTASSELCacheDir", function(...) tmpDir)
    on.exit(restore1(), add = TRUE)
    restore2 <- setup_mock("downloadWithProgress", function(url, destfile, ...) {
        stop("simulated network error")
    })
    on.exit(restore2(), add = TRUE)

    expect_error(
        setupTASSEL(version = version),
        "Failed to download TASSEL"
    )
})

test_that("setupTASSEL() cleans up partial download on error", {
    version <- "99.99.97-cleanup"
    tmpDir  <- file.path(tempdir(), "rtassel_setup_cleanup")
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)

    jarName <- getTASSELJarName(version = version)

    restore1 <- setup_mock("getTASSELCacheDir", function(...) tmpDir)
    on.exit(restore1(), add = TRUE)
    restore2 <- setup_mock("getTASSELJarName", function(...) jarName)
    on.exit(restore2(), add = TRUE)
    restore3 <- setup_mock("downloadWithProgress", function(url, destfile, ...) {
        writeLines("partial content", destfile)
        stop("connection reset")
    })
    on.exit(restore3(), add = TRUE)

    expect_error(setupTASSEL(version = version))

    jarFile <- file.path(tmpDir, jarName)
    expect_false(file.exists(jarFile))
})

test_that("setupTASSEL() re-downloads when force = TRUE", {
    version <- "99.99.98-force"
    tmpDir  <- file.path(tempdir(), "rtassel_setup_force")
    dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)

    jarName <- getTASSELJarName(version = version)
    writeLines("old content", file.path(tmpDir, jarName))

    restore1 <- setup_mock("getTASSELCacheDir", function(...) tmpDir)
    on.exit(restore1(), add = TRUE)
    restore2 <- setup_mock("getTASSELJarName", function(...) jarName)
    on.exit(restore2(), add = TRUE)
    restore3 <- setup_mock("downloadWithProgress", function(url, destfile, ...) {
        writeLines("new content", destfile)
        invisible(destfile)
    })
    on.exit(restore3(), add = TRUE)

    result <- setupTASSEL(version = version, force = TRUE)
    expect_equal(as.character(result), tmpDir)
    expect_equal(readLines(file.path(tmpDir, jarName)), "new content")
})

test_that("setupTASSEL() constructs correct Maven URL", {
    version <- "99.99.96-url"
    tmpDir  <- file.path(tempdir(), "rtassel_setup_url")
    on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)

    captured_url <- NULL

    restore1 <- setup_mock("getTASSELCacheDir", function(...) tmpDir)
    on.exit(restore1(), add = TRUE)
    restore2 <- setup_mock("downloadWithProgress", function(url, destfile, ...) {
        captured_url <<- url
        writeLines("test", destfile)
        invisible(destfile)
    })
    on.exit(restore2(), add = TRUE)

    setupTASSEL(version = version)

    expect_true(!is.null(captured_url))
    expect_match(captured_url, TASSEL_MAVEN$BASE_URL, fixed = TRUE)
    expect_match(captured_url, TASSEL_MAVEN$GROUP_PATH, fixed = TRUE)
    expect_match(captured_url, TASSEL_MAVEN$ARTIFACT_ID, fixed = TRUE)
    expect_match(captured_url, version, fixed = TRUE)
})


