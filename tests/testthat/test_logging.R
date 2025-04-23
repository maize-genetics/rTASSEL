# === Tests for logging methods =====================================

test_that("startLogger()", {
    startLogger(verbose = FALSE)
    expect_true(file.exists("rtassel_log.txt"))
})

