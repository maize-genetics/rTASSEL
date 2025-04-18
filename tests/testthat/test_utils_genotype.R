
# subEllipsis tests

test_that("subEllipsis returns ellipsis without spaces when ind is NULL", {
    res <- rTASSEL:::subEllipsis(NULL)
    # Remove ANSI codes for comparison
    expect_true(grepl(cli::symbol$ellipsis, res))
})

test_that("subEllipsis returns ellipsis with correct spacing when ind is provided", {
    ind <- 3
    res <- subEllipsis(ind)
    expect_true(grepl(cli::symbol$ellipsis, res))
})

# truncateGtName tests

test_that("truncateGtName leaves short ids unchanged", {
    id <- "abc"
    expect_equal(truncateGtName(id, maxLength = 5), id)
})

test_that("truncateGtName truncates long ids and appends ellipsis", {
    id <- "abcdef"
    res <- truncateGtName(id, maxLength = 5)
    expect_equal(res, paste0(substr(id, 1, 4), cli::symbol$ellipsis))
})

# formatRefProb tests

test_that("formatRefProb errors on non-numeric input", {
    expect_error(formatRefProb("a"), "Input must be numeric")
})

test_that("formatRefProb errors on values outside range 0–1", {
    expect_error(formatRefProb(c(-0.1, 0.5)))
    expect_error(formatRefProb(c(0.5, 1.1)))
})


test_that("formatRefProb formats NA and numeric values correctly", {
    vals <- c(NA, 0, 0.21, 0.41, 0.61, 0.81)
    res <- formatRefProb(vals)

    validAnsiCols <- c(231, 195, 159, 87, 45)
    ansiTemplate <- "\033[30;48;5;%sm %s \033[0m"
    # Test for background colors
    expect_equal(res[1], " NA ")
    expect_equal(res[2], sprintf(ansiTemplate, validAnsiCols[1], "0.000"))
    expect_equal(res[3], sprintf(ansiTemplate, validAnsiCols[2], "0.210"))
    expect_equal(res[4], sprintf(ansiTemplate, validAnsiCols[3], "0.410"))
    expect_equal(res[5], sprintf(ansiTemplate, validAnsiCols[4], "0.610"))
    expect_equal(res[6], sprintf(ansiTemplate, validAnsiCols[5], "0.810"))
})

# ANSI style functions tests

test_that("boldStyle wraps allele with ANSI bold codes", {
    allele <- "X"
    expect_equal(boldStyle(allele), sprintf("\033[1m %s \033[22m", allele))
})

test_that("bgGreenBold wraps allele correctly", {
    allele <- "X"
    expect_equal(bgGreenBold(allele),
                 sprintf("\033[42m\033[37m\033[1m %s \033[22m\033[39m\033[49m", allele))
})

# styleCache consistency

test_that("styleCache entries match style functions", {
    expect_equal(styleCache[["N"]], boldStyle("N"))
    expect_equal(styleCache[["AMaj"]], bgYellowBold("A"))
    expect_equal(styleCache[["CMin"]], bgBlueWhiteBold("C"))
})

# formatAllele tests

test_that("formatAllele returns minor style when allele equals minAllele", {
    expect_equal(formatAllele("A", "A"), styleCache[["AMin"]])
})

test_that("formatAllele returns major style when allele differs from minAllele", {
    expect_equal(formatAllele("A", "C"), styleCache[["AMaj"]])
})

test_that("formatAllele returns subEllipsis for ellipsis input", {
    el <- cli::symbol$ellipsis
    expect_equal(formatAllele(el, ""), subEllipsis(1))
})

test_that("formatAllele defaults to bold for unrecognized allele", {
    res <- formatAllele("Z", "A")
    expect_equal(res, "\033[1m Z \033[22m")
})

# genMinorAlleles tests

gt_small <- list(
    numberOfSites = function() 3,
    minorAlleleAsString = function(i) letters[i + 1]
)
test_that("genMinorAlleles returns correct vector for small gt", {
    res <- genMinorAlleles(gt_small, nSites = 10)
    expect_equal(res, c("a", "b", "c"))
})

gt_equal <- list(
    numberOfSites = function() 11,
    minorAlleleAsString = function(i) as.character(i)
)
test_that("genMinorAlleles handles maxSite == nSites + 1", {
    res <- genMinorAlleles(gt_equal, nSites = 10)
    expect_equal(res, as.character(0:10))
})

gt_large <- list(
    numberOfSites = function() 13,
    minorAlleleAsString = function(i) as.character(i)
)
test_that("genMinorAlleles handles maxSite > nSites + 2", {
    res <- genMinorAlleles(gt_large, nSites = 10)
    expect_equal(length(res), 12)
    expect_equal(res[11], cli::symbol$ellipsis)
    expect_equal(res[12], "12")
})

# addGtRowColIds tests

gt_small2 <- list(numberOfSites = function() 3)
fgs <- c("foo", "bar")
test_that("addGtRowColIds small numeric = FALSE", {
    res <- addGtRowColIds(gt_small2, fgs, nSites = 10, numeric = FALSE)
    expect_true(grepl(" 0  1  2 ", res[1]))
    expect_equal(res[2], "foo")
    expect_equal(res[3], "bar")
})


