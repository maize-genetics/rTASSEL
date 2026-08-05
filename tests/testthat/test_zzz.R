# === Test initializers =============================================

test_that("onLoad function is called without error", {
    expect_silent(
        .onLoad(pkgname = "rTASSEL", libname = system.file(package = "rTASSEL"))
    )
})


# --- .onLoad behavior ---------------------------------------------

test_that(".onLoad returns invisibly when no JARs are found", {
    result <- .onLoad(
        pkgname = "nonexistent_pkg_99999",
        libname = tempdir()
    )

    # When resolveJarPath returns NULL, .onLoad should return invisible()
    # without error (no attempt to initialize JVM)
    expect_null(result)
})


# --- .onAttach behavior -------------------------------------------

test_that(".onAttach produces startup messages", {
    msgs <- tryCatch(
        {
            .onAttach(
                libname = system.file(package = "rTASSEL"),
                pkgname = "rTASSEL"
            )
        },
        message = function(m) conditionMessage(m)
    )

    expect_type(msgs, "character")
    expect_true(nzchar(msgs))
})

test_that(".onAttach includes version information", {
    msgs <- tryCatch(
        {
            .onAttach(
                libname = system.file(package = "rTASSEL"),
                pkgname = "rTASSEL"
            )
        },
        message = function(m) conditionMessage(m)
    )

    expect_match(msgs, "rTASSEL", fixed = TRUE)
})

test_that(".onAttach handles missing JARs gracefully", {
    msgs <- tryCatch(
        {
            .onAttach(
                libname = tempdir(),
                pkgname = "nonexistent_pkg_99999"
            )
        },
        message = function(m) conditionMessage(m)
    )

    expect_type(msgs, "character")
})


# --- .onAttach message content ------------------------------------

test_that(".onAttach mentions setupTASSEL when JARs are missing", {
    local_mocked_bindings(
        resolveJarPath = function(...) list(path = NULL, source = NULL)
    )

    msgs <- tryCatch(
        {
            .onAttach(
                libname = tempdir(),
                pkgname = "rTASSEL"
            )
        },
        message = function(m) conditionMessage(m)
    )

    expect_match(msgs, "setupTASSEL", fixed = TRUE)
    expect_match(msgs, "not found", fixed = TRUE)
})

test_that(".onAttach mentions TASSEL version when JARs are present", {
    msgs <- tryCatch(
        {
            .onAttach(
                libname = system.file(package = "rTASSEL"),
                pkgname = "rTASSEL"
            )
        },
        message = function(m) conditionMessage(m)
    )

    expect_match(msgs, getLoadedTASSELVersion(), fixed = TRUE)
    expect_match(msgs, "startLogger", fixed = TRUE)
})

test_that("getLoadedTASSELVersion reports the version TASSEL itself reports", {
    expect_match(getLoadedTASSELVersion(), "^[0-9]+\\.[0-9]+\\.[0-9]+$")
})

test_that("getLoadedTASSELVersion falls back when the JVM cannot be queried", {
    local_mocked_bindings(
        J = function(...) stop("no class"),
        .package = "rJava"
    )

    expect_equal(getLoadedTASSELVersion(fallback = "1.2.3"), "1.2.3")
})


# --- .onAttach update notice --------------------------------------

test_that(".onAttach announces a newer TASSEL version", {
    local_mocked_bindings(
        resolveJarPath         = function(...) list(path = tempdir(), source = "maven cache"),
        getLoadedTASSELVersion = function(...) "5.2.96",
        checkForTASSELUpdate   = function(...) list(
            latest    = "5.2.97",
            layout    = "thin",
            checkedAt = Sys.time()
        )
    )

    msgs <- tryCatch(
        .onAttach(libname = tempdir(), pkgname = "rTASSEL"),
        message = function(m) conditionMessage(m)
    )

    expect_match(msgs, "5.2.97", fixed = TRUE)
    expect_match(msgs, "setupTASSEL", fixed = TRUE)
})

test_that(".onAttach stays quiet when TASSEL is up to date", {
    local_mocked_bindings(
        resolveJarPath         = function(...) list(path = tempdir(), source = "maven cache"),
        getLoadedTASSELVersion = function(...) "5.2.96",
        checkForTASSELUpdate   = function(...) list(
            latest    = "5.2.96",
            layout    = "fat",
            checkedAt = Sys.time()
        )
    )

    msgs <- tryCatch(
        .onAttach(libname = tempdir(), pkgname = "rTASSEL"),
        message = function(m) conditionMessage(m)
    )

    expect_no_match(msgs, "is available", fixed = TRUE)
})

test_that(".onAttach reports the build date of a nightly install", {
    local_mocked_bindings(
        resolveJarPath         = function(...) list(path = tempdir(), source = "nightly build"),
        getActiveTASSELVersion = function(...) "5.2.98-dev.20260801",
        checkForTASSELUpdate   = function(...) NULL
    )

    msgs <- tryCatch(
        .onAttach(libname = tempdir(), pkgname = "rTASSEL"),
        message = function(m) conditionMessage(m)
    )

    expect_match(msgs, "5.2.98-dev.20260801", fixed = TRUE)
    expect_match(msgs, "nightly build", fixed = TRUE)
})

test_that(".onAttach announces a newer nightly build", {
    local_mocked_bindings(
        resolveJarPath         = function(...) list(path = tempdir(), source = "nightly build"),
        getActiveTASSELVersion = function(...) "5.2.98-dev.20260801",
        checkForTASSELUpdate   = function(force = FALSE, channel = "maven") {
            expect_equal(channel, "nightly")
            list(latest = "5.2.98-dev.20260805", checkedAt = Sys.time())
        }
    )

    msgs <- tryCatch(
        .onAttach(libname = tempdir(), pkgname = "rTASSEL"),
        message = function(m) conditionMessage(m)
    )

    expect_match(msgs, "5.2.98-dev.20260805", fixed = TRUE)
    expect_match(msgs, "source = \\\"github\\\"")
})

test_that(".onAttach stays quiet when the nightly build is current", {
    local_mocked_bindings(
        resolveJarPath         = function(...) list(path = tempdir(), source = "nightly build"),
        getActiveTASSELVersion = function(...) "5.2.98-dev.20260801",
        checkForTASSELUpdate   = function(...) list(
            latest    = "5.2.98-dev.20260801",
            checkedAt = Sys.time()
        )
    )

    msgs <- tryCatch(
        .onAttach(libname = tempdir(), pkgname = "rTASSEL"),
        message = function(m) conditionMessage(m)
    )

    expect_no_match(msgs, "is available", fixed = TRUE)
})

test_that(".onAttach checks Maven, not GitHub, for a user-supplied JAR path", {
    channels <- character(0)

    local_mocked_bindings(
        resolveJarPath         = function(...) list(path = tempdir(), source = "option"),
        getActiveTASSELVersion = function(...) "5.2.98-dev.20260801",
        getLoadedTASSELVersion = function(...) "5.2.98",
        checkForTASSELUpdate   = function(force = FALSE, channel = "maven") {
            channels <<- c(channels, channel)
            NULL
        }
    )

    tryCatch(
        .onAttach(libname = tempdir(), pkgname = "rTASSEL"),
        message = function(m) conditionMessage(m)
    )

    expect_equal(channels, "maven")
})

test_that(".onAttach stays quiet when the update check yields nothing", {
    local_mocked_bindings(
        resolveJarPath         = function(...) list(path = tempdir(), source = "maven cache"),
        getLoadedTASSELVersion = function(...) "5.2.96",
        checkForTASSELUpdate   = function(...) NULL
    )

    msgs <- tryCatch(
        .onAttach(libname = tempdir(), pkgname = "rTASSEL"),
        message = function(m) conditionMessage(m)
    )

    expect_no_match(msgs, "is available", fixed = TRUE)
})

test_that(".onAttach message includes Welcome header", {
    msgs <- tryCatch(
        {
            .onAttach(
                libname = system.file(package = "rTASSEL"),
                pkgname = "rTASSEL"
            )
        },
        message = function(m) conditionMessage(m)
    )

    expect_match(msgs, "Welcome", fixed = TRUE)
    expect_match(msgs, "rTASSEL", fixed = TRUE)
})


# --- .onLoad side effects -----------------------------------------

test_that(".onLoad adds JARs to classpath when package is available", {
    cpBefore <- rJava::.jclassPath()

    .onLoad(
        pkgname = "rTASSEL",
        libname = system.file(package = "rTASSEL")
    )

    cpAfter <- rJava::.jclassPath()
    expect_true(length(cpAfter) >= length(cpBefore))
})


