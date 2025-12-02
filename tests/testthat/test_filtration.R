# === Tests for filter functions ====================================

## Preamble - load data ----

### Start logging info
startLogger()

### Load hapmap data
genoPathHMP <- system.file(
    "extdata",
    "mdp_genotype.hmp.txt",
    package = "rTASSEL"
)

### Read data - need only non missing data!
phenoPathFast <- system.file(
    "extdata",
    "mdp_traits_nomissing.txt",
    package = "rTASSEL"
)

### Create rTASSEL phenotype only object
tasPheno <- readPhenotypeFromPath(
    path = phenoPathFast
)

### Create rTASSEL genotype only object
tasGeno <- readGenotypeTableFromPath(
    path = genoPathHMP
)

### Create rTASSEL object - use prior TASSEL genotype object
tasGenoPhenoFast <- readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoPathFast
)


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

    expect_error(filterGenotypeTableSites(mtcars))
    expect_error(filterGenotypeTableSites(tasPheno))
    expect_error(filterGenotypeTableSites(tasGeno, siteMinAlleleFreq = 5))
    expect_error(filterGenotypeTableSites(tasGeno, siteMinAlleleFreq = -1))
    expect_error(filterGenotypeTableSites(tasGeno, minHeterozygous = 5))
    expect_error(filterGenotypeTableSites(tasGeno, minHeterozygous = -1))
    expect_error(filterGenotypeTableSites(tasGeno, siteRangeFilterType = "noon"))
    expect_error(filterGenotypeTableSites(tasGeno, gRangesObj = mtcars))
    expect_error(filterGenotypeTableSites(tasGeno, siteMaxAlleleFreq = 2))
    expect_error(filterGenotypeTableSites(tasGeno, maxHeterozygous = 2))
    expect_error(filterGenotypeTableSites(tasGeno, maxHeterozygous = 2))
})

test_that("filterGenotypeTableSites rejects negative site indices", {

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(path = genoPathHMP)

    # Negative startSite
    expect_error(
        object = filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "sites",
            startSite = -21,
            endSite = 300
        ),
        regexp = "startSite and endSite must be non-negative integers."
    )

    # Negative endSite
    expect_error(
        object = filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "sites",
            startSite = 0,
            endSite = -5
        ),
        regexp = "startSite and endSite must be non-negative integers."
    )

    # Both negative
    expect_error(
        object = filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "sites",
            startSite = -10,
            endSite = -5
        ),
        regexp = "startSite and endSite must be non-negative integers."
    )
})

test_that("filterGenotypeTableSites rejects negative position values", {

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(path = genoPathHMP)

    # Negative startPos
    expect_error(
        object = filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "position",
            startChr = 1,
            endChr = 5,
            startPos = -100,
            endPos = 500000
        ),
        regexp = "startPos must be a non-negative integer."
    )

    # Negative endPos
    expect_error(
        object = filterGenotypeTableSites(
            tasObj = tasGenoHMP,
            siteRangeFilterType = "position",
            startChr = 1,
            endChr = 5,
            startPos = 100,
            endPos = -500000
        ),
        regexp = "endPos must be a non-negative integer."
    )
})

test_that("filterGenotypeTableSites handles same chromosome with NULL positions", {

    # Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )
    tasGenoHMP <- readGenotypeTableFromPath(path = genoPathHMP)

    # Same startChr and endChr with NULL positions should NOT error
    # (This was previously causing "missing value where TRUE/FALSE needed")
    tasFiltTests <- filterGenotypeTableSites(
        tasObj = tasGenoHMP,
        siteRangeFilterType = "position",
        startChr = 5,
        endChr = 5,
        startPos = NULL,
        endPos = NULL
    )

    # Should return a valid TasselGenotypePhenotype object
    expect_s4_class(tasFiltTests, "TasselGenotypePhenotype")

    # Should contain sites from chromosome 5
    expect_true(tasFiltTests@jPositionList$numberOfSites() > 0)
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

    expect_equal(
        object = length(taxaList(filterGenotypeTableSites(tasGenoPhenoFast))),
        expected = length(taxaList(tasGenoPhenoFast))
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

test_that("filterGenotypeTableTaxa() returns correct exceptions", {
    expect_error(filterGenotypeTableTaxa(mtcars))
    expect_error(filterGenotypeTableTaxa(tasPheno))
    expect_error(filterGenotypeTableTaxa(tasGeno, minNotMissing = 45))
    expect_error(filterGenotypeTableTaxa(tasGeno, minNotMissing = -2))

    expect_error(filterGenotypeTableTaxa(tasGeno, minHeterozygous = 45))
    expect_error(filterGenotypeTableTaxa(tasGeno, minHeterozygous = -2))

    expect_error(filterGenotypeTableTaxa(tasGeno, maxHeterozygous = 45))
    expect_error(filterGenotypeTableTaxa(tasGeno, maxHeterozygous = -2))

    expect_error(filterGenotypeTableTaxa(tasGeno, taxa = mtcars))
    expect_error(filterGenotypeTableTaxa(tasGeno, taxa = 1:50))
})

test_that("filterGenotypeTableTaxa() returns correct data", {
    expect_equal(
        object = length(taxaList(filterGenotypeTableTaxa(tasGenoPhenoFast))),
        expected = length(taxaList(tasGenoPhenoFast))
    )
})

test_that("filterGenotypeTableBySiteName() tests", {

    # positive control - known sites
    testSitesPos <- head(siteSummary(tasGeno)$Site_Name)

    # Known and not valid sites
    testSitesExp <- c(
        head(siteSummary(tasGeno)$Site_Name),
        "not_valid_a"
    )

    # negative control - not valid sites
    testSitesNeg <- c("not_valid_a", "not_valid_b")

    expect_null(filterGenotypeTableBySiteName(tasGeno, testSitesNeg))
    expect_error(object = filterGenotypeTableBySiteName(mtcars, testSitesPos))
    expect_equal(
        object = siteSummary(
            filterGenotypeTableBySiteName(tasGeno, testSitesPos)
        )$Site_Name,
        expected = testSitesPos
    )
    expect_equal(
        object = siteSummary(
            filterGenotypeTableBySiteName(tasGeno, testSitesExp)
        )$Site_Name,
        expected = testSitesPos
    )
})




























