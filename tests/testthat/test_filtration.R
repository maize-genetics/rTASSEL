# === Tests for filter functions ====================================

## Error tests ----
test_that("filterGenotypeTableSites returns error when parameters not specified", {

    # Load hapmap data
    startLogger()
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(path = genoPathHMP)

    expect_error(
        filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "sites",
            startSite = NULL,
            endSite = 300
        )
    )

    expect_error(
        filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "sites",
            startSite = 100,
            endSite = NULL
        )
    )

    expect_error(
        filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "sites",
            startSite = NULL,
            endSite = NULL
        )
    )

    expect_error(
        filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "positions",
            startChr = NULL,
            endChr = 10
        )
    )

    expect_error(
        filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "positions",
            startChr = 1,
            endChr = NULL
        )
    )

    expect_error(
        filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "positions",
            startChr = NULL,
            endChr = NULL
        )
    )

    expect_error(
        filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "positron",
            startChr = NULL,
            endChr = NULL
        )
    )
})
