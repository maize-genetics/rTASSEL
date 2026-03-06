# === Tests for plotLD helpers and plotLD ===============================

## Helper function tests (no Java required) ----------------------------

test_that("computeSizingParams returns correct structure", {
    params <- computeSizingParams(15)
    expect_type(params, "list")
    expected_names <- c(
        "scaleFactor", "labelGap", "borderLw",
        "idxTextSize", "snpTextSize", "blockTextSize", "snpIdHeight"
    )
    expect_named(params, expected_names)
    expect_equal(params$scaleFactor, 1)
})

test_that("computeSizingParams scales with nSites", {
    small <- computeSizingParams(5)
    large <- computeSizingParams(100)

    expect_gt(large$scaleFactor, small$scaleFactor)
    expect_lt(large$borderLw, small$borderLw)
    expect_lt(large$idxTextSize, small$idxTextSize)
    expect_lt(large$snpTextSize, small$snpTextSize)
})

test_that("computePhysicalPositions handles single chromosome", {
    sites <- data.frame(
        coord = c("1_100", "1_200", "1_300"),
        locus = c("1", "1", "1"),
        pos   = c(100, 200, 300),
        stringsAsFactors = FALSE
    )
    result <- computePhysicalPositions(sites, verbose = FALSE)
    expect_equal(result, c(100, 200, 300))
})

test_that("computePhysicalPositions handles multiple chromosomes cumulatively", {
    sites <- data.frame(
        coord = c("1_100", "1_200", "2_50", "2_150"),
        locus = c("1", "1", "2", "2"),
        pos   = c(100, 200, 50, 150),
        stringsAsFactors = FALSE
    )
    result <- suppressMessages(
        computePhysicalPositions(sites, verbose = TRUE)
    )
    expect_length(result, 4)
    expect_true(all(diff(result) >= 0))
    expect_equal(result[1], 0)
})


## Integration tests (require rTASSEL / Java) --------------------------

startLogger()

genoPathHMP <- system.file(
    "extdata",
    "mdp_genotype.hmp.txt",
    package = "rTASSEL"
)

tasGeno <- readGenotypeTableFromPath(path = genoPathHMP)

tasGenoSmall <- filterGenotypeTableSites(
    tasObj              = tasGeno,
    siteRangeFilterType = "sites",
    startSite           = 0,
    endSite             = 14
)

ldResSmall <- linkageDiseq(
    tasObj  = tasGenoSmall,
    ldType  = "All",
    verbose = FALSE
)


test_that("plotLD returns ggplot", {
    p <- plotLD(
        ldObj        = ldResSmall,
        genomicTrack = FALSE,
        verbose      = FALSE
    )
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("plotLD returns ggplot with genomic track", {
    p <- plotLD(
        ldObj        = ldResSmall,
        genomicTrack = TRUE,
        verbose      = FALSE
    )
    expect_s3_class(p, "gg")
    expect_s3_class(p, "ggplot")
})

test_that("plotLD supports different color schemes", {
    p_viridis <- plotLD(
        ldObj       = ldResSmall,
        colorScheme = "viridis",
        verbose     = FALSE
    )
    expect_s3_class(p_viridis, "ggplot")

    p_haplo <- plotLD(
        ldObj       = ldResSmall,
        colorScheme = "haploview",
        verbose     = FALSE
    )
    expect_s3_class(p_haplo, "ggplot")
})

test_that("plotLD handles ldBlocks parameter", {
    ldDF   <- ldResSmall@results
    chr    <- as.character(ldDF$Locus1[1])
    allPos <- sort(unique(c(
        as.numeric(ldDF$Position1),
        as.numeric(ldDF$Position2)
    )))

    blocks <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges   = IRanges::IRanges(
            start = min(allPos),
            end   = max(allPos)
        ),
        label = "TestBlock"
    )

    p <- plotLD(
        ldObj        = ldResSmall,
        ldBlocks     = blocks,
        genomicTrack = FALSE,
        verbose      = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotLD with ldBlocks and genomicTrack", {
    ldDF   <- ldResSmall@results
    chr    <- as.character(ldDF$Locus1[1])
    allPos <- sort(unique(c(
        as.numeric(ldDF$Position1),
        as.numeric(ldDF$Position2)
    )))

    blocks <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges   = IRanges::IRanges(
            start = min(allPos),
            end   = max(allPos)
        ),
        label = "TestBlock"
    )

    p <- plotLD(
        ldObj        = ldResSmall,
        ldBlocks     = blocks,
        genomicTrack = TRUE,
        verbose      = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotLD errors on non-LDResults input", {
    expect_error(
        plotLD(
            ldObj   = mtcars,
            verbose = FALSE
        ),
        "LDResults"
    )
})

test_that("plotLD errors on invalid ldBlocks", {
    expect_error(
        plotLD(
            ldObj    = ldResSmall,
            ldBlocks = data.frame(x = 1),
            verbose  = FALSE
        ),
        "LDRegion"
    )
})


## resolveBlocks helper tests (no Java required) -----------------------

test_that("resolveBlocks handles single LDRegion", {
    sites <- data.frame(
        coord = c("1_100", "1_200", "1_300"),
        locus = c("1", "1", "1"),
        pos   = c(100, 200, 300),
        stringsAsFactors = FALSE
    )
    r <- LDRegion(start = 100, end = 300, label = "A", color = "red", linewidth = 2)
    df <- resolveBlocks(r, sites)
    expect_equal(nrow(df), 1)
    expect_equal(df$chr, "1")
    expect_equal(df$start, 100)
    expect_equal(df$end, 300)
    expect_equal(df$label, "A")
    expect_equal(df$color, "red")
    expect_equal(df$linewidth, 2)
    expect_true(df$showSpan)
})

test_that("resolveBlocks handles list of LDRegion", {
    sites <- data.frame(
        coord = c("1_100", "1_200", "1_300"),
        locus = c("1", "1", "1"),
        pos   = c(100, 200, 300),
        stringsAsFactors = FALSE
    )
    blocks <- list(
        LDRegion(start = 100, end = 200, label = "A"),
        LDRegion(start = 200, end = 300, color = "blue", showSpan = FALSE)
    )
    df <- resolveBlocks(blocks, sites)
    expect_equal(nrow(df), 2)
    expect_equal(df$color, c("black", "blue"))
    expect_true(is.na(df$label[2]))
    expect_equal(df$showSpan, c(TRUE, FALSE))
})

test_that("resolveBlocks errors on multi-chromosome with LDRegion", {
    sites <- data.frame(
        coord = c("1_100", "2_200"),
        locus = c("1", "2"),
        pos   = c(100, 200),
        stringsAsFactors = FALSE
    )
    r <- LDRegion(start = 100, end = 200)
    expect_error(resolveBlocks(r, sites), "single-chromosome")
})


## Integration tests: plotLD with LDRegion (require rTASSEL / Java) ----

test_that("plotLD handles single LDRegion", {
    ldDF   <- ldResSmall@results
    allPos <- sort(unique(c(
        as.numeric(ldDF$Position1),
        as.numeric(ldDF$Position2)
    )))

    block <- LDRegion(
        start = min(allPos),
        end   = max(allPos),
        label = "TestBlock"
    )

    p <- plotLD(
        ldObj    = ldResSmall,
        ldBlocks = block,
        verbose  = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotLD handles list of LDRegion with custom colors", {
    ldDF   <- ldResSmall@results
    allPos <- sort(unique(c(
        as.numeric(ldDF$Position1),
        as.numeric(ldDF$Position2)
    )))
    mid <- median(allPos)

    blocks <- list(
        LDRegion(start = min(allPos), end = mid, label = "Left", color = "blue"),
        LDRegion(start = mid, end = max(allPos), label = "Right", color = "red")
    )

    p <- plotLD(
        ldObj    = ldResSmall,
        ldBlocks = blocks,
        verbose  = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotLD with LDRegion and genomicTrack", {
    ldDF   <- ldResSmall@results
    allPos <- sort(unique(c(
        as.numeric(ldDF$Position1),
        as.numeric(ldDF$Position2)
    )))

    block <- LDRegion(
        start = min(allPos),
        end   = max(allPos),
        label = "TestBlock",
        color = "#FF5733"
    )

    p <- plotLD(
        ldObj        = ldResSmall,
        ldBlocks     = block,
        genomicTrack = TRUE,
        verbose      = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotLD with LDRegion custom linewidth", {
    ldDF   <- ldResSmall@results
    allPos <- sort(unique(c(
        as.numeric(ldDF$Position1),
        as.numeric(ldDF$Position2)
    )))

    block <- LDRegion(
        start     = min(allPos),
        end       = max(allPos),
        label     = "Thick",
        linewidth = 2.5
    )

    p <- plotLD(
        ldObj    = ldResSmall,
        ldBlocks = block,
        verbose  = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotLD with LDRegion showSpan = FALSE suppresses span", {
    ldDF   <- ldResSmall@results
    allPos <- sort(unique(c(
        as.numeric(ldDF$Position1),
        as.numeric(ldDF$Position2)
    )))

    block <- LDRegion(
        start    = min(allPos),
        end      = max(allPos),
        label    = "NoSpan",
        showSpan = FALSE
    )

    p <- plotLD(
        ldObj    = ldResSmall,
        ldBlocks = block,
        verbose  = FALSE
    )
    expect_s3_class(p, "ggplot")
})

test_that("plotLD with LDRegion showSpan = TRUE and no label shows span only", {
    ldDF   <- ldResSmall@results
    allPos <- sort(unique(c(
        as.numeric(ldDF$Position1),
        as.numeric(ldDF$Position2)
    )))

    block <- LDRegion(
        start    = min(allPos),
        end      = max(allPos),
        showSpan = TRUE
    )

    p <- plotLD(
        ldObj    = ldResSmall,
        ldBlocks = block,
        verbose  = FALSE
    )
    expect_s3_class(p, "ggplot")
})
