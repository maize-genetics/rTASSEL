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


