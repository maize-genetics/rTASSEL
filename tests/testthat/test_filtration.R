# === Tests for filter functions ====================================

## Error tests ----
test_that("filterGenotypeTableSites returns error when parameters not specified", {

    # Load hapmap data
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
            siteRangeFilterType = "position",
            startChr = NULL,
            endChr = 10
        )
    )

    expect_error(
        filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "position",
            startChr = 1,
            endChr = NULL
        )
    )

    expect_error(
        filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "position",
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


## Equality tests ----
test_that("filterGenotypeTableSites returns correct positions.", {

    # Expected (observed) data
    obs_taxa <- 281
    obs_site <- 1486

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(path = genoPathHMP)

    # Filter
    tasFiltTests <- filterGenotypeTableSites(
        tasObj = tasGenoHMP,
        siteRangeFilterType = "position",
        startChr = 5,
        endChr = 10
    )

    # Test data
    test_taxa <- tasFiltTests@jTaxaList$numberOfTaxa()
    test_site <- tasFiltTests@jPositionList$numberOfSites()

    # Equality test
    expect_equal(
        object = c(test_taxa, test_site),
        expected = c(obs_taxa, obs_site)
    )

})

test_that("filterGenotypeTableSites returns correct positions.", {

    # Expected (observed) data
    obs_taxa <- 281
    obs_site <- 2896

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(path = genoPathHMP)

    # Filter
    tasFiltTests <- filterGenotypeTableSites(
        tasObj = tasGenoHMP,
        siteRangeFilterType = "position",
        startChr = 1,
        endChr = 10,
        startPos = 400000,
        endPos = 5000000
    )

    # Test data
    test_taxa <- tasFiltTests@jTaxaList$numberOfTaxa()
    test_site <- tasFiltTests@jPositionList$numberOfSites()

    # Equality test
    expect_equal(
        object = c(test_taxa, test_site),
        expected = c(obs_taxa, obs_site)
    )

})

test_that("filterGenotypeTableSites returns correct sites.", {

    # Expected (observed) data
    obs_taxa <- 281
    obs_site <- 300

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(path = genoPathHMP)

    # Filter
    tasFiltTests <- filterGenotypeTableSites(
        tasObj = tasGenoHMP,
        siteRangeFilterType = "sites",
        startSite = 1,
        endSite = 300
    )

    # Test data
    test_taxa <- tasFiltTests@jTaxaList$numberOfTaxa()
    test_site <- tasFiltTests@jPositionList$numberOfSites()

    # Equality test
    expect_equal(
        object = c(test_taxa, test_site),
        expected = c(obs_taxa, obs_site)
    )

})

test_that("filterGenotypeTableSites returns correct sites with no range filtration.", {

    # Expected (observed) data
    obs_taxa <- 281
    obs_site <- 3093

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(path = genoPathHMP)

    # Filter
    tasFiltTests <- filterGenotypeTableSites(
        tasObj = tasGenoHMP,
        siteRangeFilterType = "none",
        startSite = 1,
        endSite = 300
    )

    # Test data
    test_taxa <- tasFiltTests@jTaxaList$numberOfTaxa()
    test_site <- tasFiltTests@jPositionList$numberOfSites()

    # Equality test
    expect_equal(
        object = c(test_taxa, test_site),
        expected = c(obs_taxa, obs_site)
    )

})


