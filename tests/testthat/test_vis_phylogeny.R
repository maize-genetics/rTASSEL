# === Tests for phylogenetic visualizations =========================

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

### Filter object for further tests
filterGenoObj <- filterGenotypeTableSites(
    tasObj = tasGeno,
    siteRangeFilterType = "sites",
    startSite = 0,
    endSite = 10
)
filterGenoObj <- filterGenotypeTableTaxa(
    tasObj = filterGenoObj,
    taxa = taxaList(tasGeno)[grep("^[0-9]|^A", taxaList(tasGeno))]
)


test_that("treeJavaApp() returns the correct exceptions", {
    expect_error(treeJavaApp(mtcars))
    expect_error(treeJavaApp(tasPheno))
    expect_error(treeJavaApp(filterGenoObj, clustMethod = "bad option"))

    # TODO - fix Forester exception
    expect_error(treeJavaApp(filterGenoObj))
    expect_error(treeJavaApp(filterGenoObj,clustMethod = "Neighbor_Joining"))
})


