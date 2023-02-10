# === Tests for logging methods =====================================

test_that("startLogger()", {
    expect_error(startLogger("~/test_dir"))
})

