# === Tests for LDResults S4 class ==========================================


## Minimal valid data frame for constructing LDResults ---------------------

make_ld_df <- function(n = 3) {
    data.frame(
        Locus1          = rep("1", n),
        Position1       = as.character(seq(100, by = 100, length.out = n)),
        Site1           = seq_len(n),
        NumberOfStates1 = rep(2L, n),
        States1         = rep("A/C", n),
        Frequency1      = rep("0.5", n),
        Locus2          = rep("1", n),
        Position2       = as.character(seq(200, by = 100, length.out = n)),
        Site2           = seq_len(n) + 1L,
        NumberOfStates2 = rep(2L, n),
        States2         = rep("A/C", n),
        Frequency2      = rep("0.5", n),
        Dist_bp         = rep("100", n),
        `R^2`           = runif(n),
        DPrime          = runif(n),
        pDiseq          = runif(n),
        N               = rep(100L, n),
        stringsAsFactors = FALSE,
        check.names      = FALSE
    )
}


## Constructor / basic slots -----------------------------------------------

test_that("LDResults can be constructed with valid inputs", {
    df  <- make_ld_df()
    obj <- methods::new(
        "LDResults",
        results    = df,
        ldType     = "All",
        windowSize = -1,
        hetCalls   = "missing"
    )
    expect_s4_class(obj, "LDResults")
    expect_equal(obj@ldType, "All")
    expect_equal(obj@windowSize, -1)
    expect_equal(obj@hetCalls, "missing")
    expect_equal(nrow(obj@results), 3)
})

test_that("LDResults accepts SlidingWindow ldType", {
    df  <- make_ld_df()
    obj <- methods::new(
        "LDResults",
        results    = df,
        ldType     = "SlidingWindow",
        windowSize = 10,
        hetCalls   = "ignore"
    )
    expect_equal(obj@ldType, "SlidingWindow")
    expect_equal(obj@windowSize, 10)
    expect_equal(obj@hetCalls, "ignore")
})

test_that("LDResults accepts NA windowSize", {
    df  <- make_ld_df()
    obj <- methods::new(
        "LDResults",
        results    = df,
        ldType     = "All",
        windowSize = NA_real_,
        hetCalls   = "third"
    )
    expect_true(is.na(obj@windowSize))
})


## Validity: required columns ---------------------------------------------

test_that("LDResults rejects data.frame missing required columns", {
    df <- data.frame(Locus1 = "1", Position1 = "100")
    expect_error(
        methods::new(
            "LDResults",
            results    = df,
            ldType     = "All",
            windowSize = -1,
            hetCalls   = "missing"
        ),
        "missing required columns"
    )
})


## Validity: ldType --------------------------------------------------------

test_that("LDResults rejects invalid ldType", {
    df <- make_ld_df()
    expect_error(
        methods::new(
            "LDResults",
            results    = df,
            ldType     = "Everything",
            windowSize = -1,
            hetCalls   = "missing"
        ),
        "ldType"
    )
})


## Validity: hetCalls ------------------------------------------------------

test_that("LDResults rejects invalid hetCalls", {
    df <- make_ld_df()
    expect_error(
        methods::new(
            "LDResults",
            results    = df,
            ldType     = "All",
            windowSize = -1,
            hetCalls   = "bogus"
        ),
        "hetCalls"
    )
})


## Show method -------------------------------------------------------------

test_that("LDResults show method prints key fields", {
    df  <- make_ld_df(5)
    obj <- methods::new(
        "LDResults",
        results    = df,
        ldType     = "SlidingWindow",
        windowSize = 10,
        hetCalls   = "missing"
    )
    out <- cli::ansi_strip(capture.output(show(obj)))

    expect_true(any(grepl("LDResults", out)))
    expect_true(any(grepl("5 pairs", out)))
    expect_true(any(grepl("SlidingWindow", out)))
    expect_true(any(grepl("10", out)))
    expect_true(any(grepl("missing", out)))
})


## tableReport method ------------------------------------------------------

test_that("tableReport returns a tibble for LDResults", {
    df  <- make_ld_df()
    obj <- methods::new(
        "LDResults",
        results    = df,
        ldType     = "All",
        windowSize = -1,
        hetCalls   = "missing"
    )
    tbl <- tableReport(obj)
    expect_s3_class(tbl, "tbl_df")
    expect_equal(nrow(tbl), 3)
    expect_true("R^2" %in% colnames(tbl))
})
