# === Tests for SNP density plot functions ===========================

## Input validation ----
test_that("plotSnpDensity rejects invalid inputs", {
    expect_error(
        plotSnpDensity(mtcars),
        "tasObj must be of class"
    )
    expect_error(
        plotSnpDensity(rtObjs$gt_hmp, windowSize = -100),
        "windowSize"
    )
    expect_error(
        plotSnpDensity(rtObjs$gt_hmp, windowSize = "abc"),
        "windowSize"
    )
    expect_error(
        plotSnpDensity(rtObjs$gt_hmp, colorOption = "neon_pink")
    )
})


## TasselGenotypePhenotype input ----
test_that("plotSnpDensity works with TasselGenotypePhenotype", {
    p <- plotSnpDensity(rtObjs$gt_hmp)
    expect_s3_class(p, "gg")

    pBld <- ggplot2::ggplot_build(p)
    expect_s3_class(p$layers[[1]]$geom, "GeomTile")
    expect_equal(p$labels$y, "Chromosome")
    expect_true(grepl("Position", p$labels$x))
})


## TasselGenotype input ----
test_that("plotSnpDensity works with TasselGenotype", {
    gtObj <- readGenotype(rtFiles$gt_hmp_path)
    expect_s4_class(gtObj, "TasselGenotype")

    p <- plotSnpDensity(gtObj)
    expect_s3_class(p, "gg")
    expect_s3_class(p$layers[[1]]$geom, "GeomTile")
})


## Color options ----
test_that("plotSnpDensity respects colorOption parameter", {
    pDefault <- plotSnpDensity(rtObjs$gt_hmp, colorOption = "viridis")
    pMagma   <- plotSnpDensity(rtObjs$gt_hmp, colorOption = "magma")

    expect_s3_class(pDefault, "gg")
    expect_s3_class(pMagma, "gg")

    bldDefault <- ggplot2::ggplot_build(pDefault)
    bldMagma   <- ggplot2::ggplot_build(pMagma)

    fillsDefault <- unique(bldDefault$data[[1]]$fill)
    fillsMagma   <- unique(bldMagma$data[[1]]$fill)
    expect_false(identical(fillsDefault, fillsMagma))
})


## Interactive output ----
test_that("plotSnpDensity returns plotly when interactive = TRUE", {
    p <- plotSnpDensity(rtObjs$gt_hmp, interactive = TRUE)
    expect_s3_class(p, "plotly")
    expect_s3_class(p, "htmlwidget")
})


## primeSnpDensityData structure ----
test_that("primeSnpDensityData returns correct data frame structure", {
    params <- list(
        "tasObj"     = rtObjs$gt_hmp,
        "windowSize" = 1e6
    )
    df <- primeSnpDensityData(params)

    expect_s3_class(df, "data.frame")
    expect_true(all(c("snpCount", "windowStart", "windowMid", "Chr") %in% colnames(df)))
    expect_s3_class(df$Chr, "factor")
    expect_true(all(df$snpCount >= 0))
    expect_true(all(df$windowStart >= 0))
})


## Window bins start at 0 per chromosome ----
test_that("primeSnpDensityData window bins start at 0 for each chromosome", {
    params <- list(
        "tasObj"     = rtObjs$gt_hmp,
        "windowSize" = 1e6
    )
    df <- primeSnpDensityData(params)

    minStarts <- tapply(df$windowStart, df$Chr, min)
    expect_true(all(minStarts == 0))
})


## Zero-count windows are present ----
test_that("primeSnpDensityData includes zero-count windows within range", {
    params <- list(
        "tasObj"     = rtObjs$gt_hmp,
        "windowSize" = 1e5
    )
    df <- primeSnpDensityData(params)

    expect_true(any(df$snpCount == 0))
})


## Window size affects binning ----
test_that("primeSnpDensityData produces different bins for different window sizes", {
    dfLarge <- primeSnpDensityData(list("tasObj" = rtObjs$gt_hmp, "windowSize" = 1e7))
    dfSmall <- primeSnpDensityData(list("tasObj" = rtObjs$gt_hmp, "windowSize" = 1e5))

    expect_true(nrow(dfSmall) > nrow(dfLarge))
})


## X axis scale adapts to genomic range ----
test_that("plotSnpDensityCore adapts axis label to genomic range", {
    pLarge <- plotSnpDensity(rtObjs$gt_hmp, windowSize = 1e6)
    expect_true(grepl("Mbp", pLarge$labels$x))
})


## Zero counts render as NA (grey) ----
test_that("plotSnpDensityCore converts zero counts to NA for grey fill", {
    p   <- plotSnpDensity(rtObjs$gt_hmp, windowSize = 1e5)
    bld <- ggplot2::ggplot_build(p)

    fills <- bld$data[[1]]$fill
    expect_true(any(fills == "grey95"))
})


## logNorm parameter ----
test_that("plotSnpDensity rejects invalid logNorm values", {
    expect_error(
        plotSnpDensity(rtObjs$gt_hmp, logNorm = "yes"),
        "logNorm"
    )
    expect_error(
        plotSnpDensity(rtObjs$gt_hmp, logNorm = c(TRUE, FALSE)),
        "logNorm"
    )
})

test_that("plotSnpDensity with logNorm = TRUE produces valid plot", {
    p <- plotSnpDensity(rtObjs$gt_hmp, logNorm = TRUE)
    expect_s3_class(p, "gg")
    expect_s3_class(p$layers[[1]]$geom, "GeomTile")
})

test_that("logNorm applies log10 transform to fill values", {
    pRaw <- plotSnpDensity(rtObjs$gt_hmp, windowSize = 1e5)
    pLog <- plotSnpDensity(rtObjs$gt_hmp, windowSize = 1e5, logNorm = TRUE)

    bldRaw <- ggplot2::ggplot_build(pRaw)
    bldLog <- ggplot2::ggplot_build(pLog)

    rawFills <- bldRaw$data[[1]]$fill
    logFills <- bldLog$data[[1]]$fill

    expect_false(identical(rawFills, logFills))
})

test_that("logNorm updates the legend label", {
    pRaw <- plotSnpDensity(rtObjs$gt_hmp)
    pLog <- plotSnpDensity(rtObjs$gt_hmp, logNorm = TRUE)

    rawScales <- pRaw$scales$scales
    logScales <- pLog$scales$scales

    rawFillScale <- Filter(function(s) "fill" %in% s$aesthetics, rawScales)[[1]]
    logFillScale <- Filter(function(s) "fill" %in% s$aesthetics, logScales)[[1]]

    expect_equal(rawFillScale$name, "SNP Count")
    expect_true(is.expression(logFillScale$name))
    expect_equal(as.character(logFillScale$name), "log[10](SNP ~ Count)")
})

test_that("logNorm works with interactive mode", {
    p <- plotSnpDensity(rtObjs$gt_hmp, logNorm = TRUE, interactive = TRUE)
    expect_s3_class(p, "plotly")
})


