# === Test initializers =============================================

test_that("onLoad function is called without error", {
    expect_silent(
        .onLoad(pkgname = "rTASSEL", libname = system.file(package = "rTASSEL"))
    )
})


