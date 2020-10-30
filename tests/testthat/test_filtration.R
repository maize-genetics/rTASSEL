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

    expect_error(
        filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            startChr = 5,
            endChr = 50
        )
    )

    expect_error(
        object = filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteMinCount = 282
        ),
        regexp = "Minimum number of taxa exceeds total number of taxa in genotype table."
    )

    expect_error(
        object = filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "position",
            startChr = 1,
            endChr = 1,
            startPos = 150,
            endPos = 140
        ),
        regexp = "Filtration paramaters outside acceptable range."
    )

    expect_error(
        object = filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "sites",
            startSite = 0,
            endSite = 4000
        ),
        regexp = "End site parameter exceeds total number of sites in genotype table."
    )

    expect_error(
        object = filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "sites",
            startSite = 500,
            endSite = 300
        ),
        regexp = "Start site cannot be larger than end site."
    )
})


## Equality tests ----
test_that("filterGenotypeTableSites returns correct positions.", {

    # Expected (observed) data
    obs_taxa <- 281
    obs_site <- 139

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
        siteMinAlleleFreq = 0,
        siteMaxAlleleFreq = 0
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
    obs_site <- 2957

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
        maxHeterozygous = 0.05
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
    obs_site <- 136

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
        minHeterozygous = 0.05
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
    obs_site <- 2559

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
        siteMinAlleleFreq = 0.05
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
    obs_site <- 297

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
        siteMinCount = 280
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
    obs_site <- 534

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
        siteMaxAlleleFreq = 0.05
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

test_that("filterGenotypeTableSites returns NA if filtration is too strict.", {

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
        minHeterozygous = 1,
        siteRangeFilterType = "position",
        startChr = 1,
        startPos = 150,
        endChr = 1,
        endPos = 200000
    )

    # Equality test
    expect_equal(
        object = tasFiltTests,
        expected = NA
    )
})

test_that("filterGenotypeTableSites correct chromosome with strings as parameters.", {

    # Expected (observed) data
    obs_taxa <- 1026
    obs_site <- 186

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "wheat_gbs_test.vcf",
        package = "rTASSEL"
    )
    tasWheatGeno <- readGenotypeTableFromPath(path = genoPathHMP)

    # Filter
    tasFiltTests <- filterGenotypeTableSites(
        tasObj = tasWheatGeno,
        siteRangeFilterType = "position",
        startChr = "1A",
        endChr = "3A"
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

test_that("filterGenotypeTableSites correct chromosome with strings as parameters.", {

    # Expected (observed) data
    obs_taxa <- 1026
    obs_site <- 159

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "wheat_gbs_test.vcf",
        package = "rTASSEL"
    )
    tasWheatGeno <- readGenotypeTableFromPath(path = genoPathHMP)

    # Filter
    tasFiltTests <- filterGenotypeTableSites(
        tasObj = tasWheatGeno,
        siteRangeFilterType = "position",
        startChr = "1A",
        endChr = "2D"
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

test_that("filterGenotypeTableSites correct chromosome with strings as parameters.", {

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
        startChr = "1",
        endChr = "10",
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

test_that("filterGenotypeTableSites properly filters by a GRanges object.", {

    # Expected (observed) data
    obs_taxa <- 281
    obs_site <- 3

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(path = genoPathHMP)

    # GRanges object
    gr <- GenomicRanges::GRanges(
        seqnames = c("1"),
        ranges = IRanges::IRanges(start = c(5353318), end = c(5562503))
    )

    # Filter
    tasFiltTests <- filterGenotypeTableSites(
        tasObj = tasGenoHMP,
        gRangesObj = gr
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

test_that("filterGenotypeTableSites properly filters by a GRanges object with other parameters.", {

    # Expected (observed) data
    obs_taxa <- 281
    obs_site <- 1

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(path = genoPathHMP)

    # GRanges object
    gr <- GenomicRanges::GRanges(
        seqnames = c("1"),
        ranges = IRanges::IRanges(start = c(5353318), end = c(5562503))
    )

    # Filter
    tasFiltTests <- filterGenotypeTableSites(
        tasObj = tasGenoHMP,
        siteMaxAlleleFreq = 0.05,
        gRangesObj = gr
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
