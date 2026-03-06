# === Tests for LDRegion S4 class ==========================================


## Constructor defaults ----------------------------------------------------

test_that("LDRegion constructor creates valid object with defaults", {
    r <- LDRegion(start = 100, end = 500)
    expect_s4_class(r, "LDRegion")
    expect_equal(r@start, 100)
    expect_equal(r@end, 500)
    expect_true(is.na(r@label))
    expect_equal(r@color, "black")
    expect_true(is.na(r@linewidth))
    expect_true(r@showSpan)
})


## Constructor with all slots ----------------------------------------------

test_that("LDRegion constructor accepts all slots", {
    r <- LDRegion(
        start     = 100,
        end       = 500,
        label     = "Block A",
        color     = "blue",
        linewidth = 1.5,
        showSpan  = FALSE
    )
    expect_equal(r@start, 100)
    expect_equal(r@end, 500)
    expect_equal(r@label, "Block A")
    expect_equal(r@color, "blue")
    expect_equal(r@linewidth, 1.5)
    expect_false(r@showSpan)
})

test_that("LDRegion coerces integer inputs to numeric", {
    r <- LDRegion(start = 100L, end = 500L)
    expect_type(r@start, "double")
    expect_type(r@end, "double")
})

test_that("LDRegion accepts hex color codes", {
    r <- LDRegion(start = 100, end = 500, color = "#FF5733")
    expect_equal(r@color, "#FF5733")
})

test_that("LDRegion accepts rgb() color strings", {
    rgbCol <- grDevices::rgb(0.5, 0.2, 0.8)
    r <- LDRegion(start = 100, end = 500, color = rgbCol)
    expect_equal(r@color, rgbCol)
})


## Validity: start / end --------------------------------------------------

test_that("LDRegion rejects end < start", {
    expect_error(LDRegion(start = 500, end = 100), "greater than or equal")
})

test_that("LDRegion allows start == end (zero-width region)", {
    r <- LDRegion(start = 300, end = 300)
    expect_equal(r@start, r@end)
})

test_that("LDRegion rejects non-finite start/end", {
    expect_error(LDRegion(start = NA, end = 100), "finite numeric")
    expect_error(LDRegion(start = Inf, end = 100), "finite numeric")
    expect_error(LDRegion(start = -Inf, end = 100), "finite numeric")
    expect_error(LDRegion(start = NaN, end = 100), "finite numeric")
    expect_error(LDRegion(start = 100, end = NA), "finite numeric")
    expect_error(LDRegion(start = 100, end = Inf), "finite numeric")
})


## Validity: color ---------------------------------------------------------

test_that("LDRegion rejects invalid color string", {
    expect_error(LDRegion(start = 100, end = 500, color = "notacolor"), "not a valid R color")
})

test_that("LDRegion rejects NA color", {
    expect_error(LDRegion(start = 100, end = 500, color = NA_character_), "non-NA")
})


## Validity: linewidth -----------------------------------------------------

test_that("LDRegion rejects non-positive linewidth", {
    expect_error(LDRegion(start = 100, end = 500, linewidth = 0), "positive")
    expect_error(LDRegion(start = 100, end = 500, linewidth = -1), "positive")
})

test_that("LDRegion accepts NA linewidth (auto)", {
    r <- LDRegion(start = 100, end = 500)
    expect_true(is.na(r@linewidth))
})

test_that("LDRegion accepts positive linewidth", {
    r <- LDRegion(start = 100, end = 500, linewidth = 0.5)
    expect_equal(r@linewidth, 0.5)
})


## Validity: showSpan ------------------------------------------------------

test_that("LDRegion rejects NA showSpan", {
    expect_error(LDRegion(start = 100, end = 500, showSpan = NA), "TRUE or FALSE")
})

test_that("LDRegion accepts TRUE/FALSE showSpan", {
    r1 <- LDRegion(start = 100, end = 500, showSpan = TRUE)
    r2 <- LDRegion(start = 100, end = 500, showSpan = FALSE)
    expect_true(r1@showSpan)
    expect_false(r2@showSpan)
})


## Show method -------------------------------------------------------------

test_that("LDRegion show method prints header and all fields", {
    r <- LDRegion(start = 100000, end = 500000, label = "Test")
    out <- cli::ansi_strip(capture.output(show(r)))

    expect_true(any(grepl("LDRegion", out)))
    expect_true(any(grepl("Range", out)))
    expect_true(any(grepl("100,000.*500,000.*bp", out)))
    expect_true(any(grepl("kbp", out)))
    expect_true(any(grepl("---", out)))
    expect_true(any(grepl("Label.*Test", out)))
    expect_true(any(grepl("Color.*black", out)))
    expect_true(any(grepl("Linewidth", out)))
    expect_true(any(grepl("Show span", out)))
})

test_that("LDRegion show method formats span in kbp/Mbp", {
    r_kbp <- LDRegion(start = 100, end = 500000)
    out_kbp <- cli::ansi_strip(capture.output(show(r_kbp)))
    expect_true(any(grepl("kbp", out_kbp)))

    r_mbp <- LDRegion(start = 100, end = 2000000)
    out_mbp <- cli::ansi_strip(capture.output(show(r_mbp)))
    expect_true(any(grepl("Mbp", out_mbp)))
})

test_that("LDRegion show method displays NA label as grey NA", {
    r <- LDRegion(start = 100, end = 500)
    out <- cli::ansi_strip(capture.output(show(r)))
    expect_true(any(grepl("Label.*NA", out)))
})

test_that("LDRegion show method displays linewidth auto vs explicit", {
    r_auto <- LDRegion(start = 100, end = 500)
    out_auto <- cli::ansi_strip(capture.output(show(r_auto)))
    expect_true(any(grepl("Linewidth.*auto", out_auto)))

    r_explicit <- LDRegion(start = 100, end = 500, linewidth = 2.0)
    out_explicit <- cli::ansi_strip(capture.output(show(r_explicit)))
    expect_true(any(grepl("Linewidth.*2.*mm", out_explicit)))
})

test_that("LDRegion show method displays showSpan value", {
    r_true <- LDRegion(start = 100, end = 500, showSpan = TRUE)
    out_t <- cli::ansi_strip(capture.output(show(r_true)))
    expect_true(any(grepl("Show span.*TRUE", out_t)))

    r_false <- LDRegion(start = 100, end = 500, showSpan = FALSE)
    out_f <- cli::ansi_strip(capture.output(show(r_false)))
    expect_true(any(grepl("Show span.*FALSE", out_f)))
})

test_that("LDRegion show method displays comma-formatted positions on one line", {
    r <- LDRegion(start = 1234567, end = 9876543)
    out <- cli::ansi_strip(capture.output(show(r)))
    expect_true(any(grepl("1,234,567 - 9,876,543 bp", out)))
})
