context("Logging utilities")

test_that("startLogger()", {
    expect_error(startLogger("~/test_dir"))
})

