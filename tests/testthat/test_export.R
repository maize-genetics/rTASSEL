# === Tests for export functions ====================================
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
tasGenoPheno <- readGenotypePhenotype(
    genoPathOrObj = genoPathHMP,
    phenoPathDFOrObj = phenoPathFast
)


## Error tests ----
test_that("exportGenotypeTable() returns errors", {
    expect_error(
        object = exportGenotypeTable(
            tasObj = mtcars,
            file = "test"
        ),
        regexp = "`tasObj` must be of class `TasselGenotypePhenotype`"
    )

    expect_error(
        object = exportGenotypeTable(
            tasObj = tasPheno,
            file = "test"
        ),
        regexp = "TASSEL genotype object not found"
    )

    expect_error(
        object = exportGenotypeTable(
            tasObj = tasGeno,
            file = ""
        ),
        regexp = "File name not specified."
    )

    expect_error(
        object = exportGenotypeTable(
            tasObj = tasGeno,
            file = "my_gt",
            format = "csv"
        ),
        regexp = "'arg' should be one of “vcf”, “hapmap”, “plink”, “flapjack”, “hdf5”"
    )
})


## Equality tests ----
test_that("exportGenotypeTable() writes correct file type.", {

    exportGenotypeTable(
        tasObj = tasGeno,
        file = "my_gt"
    )

    fileID <- "my_gt.vcf"

    expect_true(file.exists(fileID))

    file.remove(fileID)
})

test_that("exportGenotypeTable() writes correct file type.", {

    exportGenotypeTable(
        tasObj = tasGeno,
        file = "my_gt",
        format = "vcf"
    )

    fileID <- "my_gt.vcf"

    expect_true(file.exists(fileID))

    file.remove(fileID)
})

test_that("exportGenotypeTable() writes correct file type.", {

    exportGenotypeTable(
        tasObj = tasGeno,
        file = "my_gt",
        format = "hapmap"
    )

    fileID <- "my_gt.hmp.txt"

    expect_true(file.exists(fileID))

    file.remove(fileID)
})

# TODO - write better HDF5 test...

test_that("exportGenotypeTable() writes correct file type.", {

    exportGenotypeTable(
        tasObj = tasGeno,
        file = "my_gt",
        format = "plink"
    )

    fileID1 <- "my_gt.plk.map"
    fileID2 <- "my_gt.plk.ped"

    expect_true(all(file.exists(fileID1), file.exists(fileID2)))

    file.remove(fileID1)
    file.remove(fileID2)
})

test_that("exportGenotypeTable() writes correct file type.", {

    exportGenotypeTable(
        tasObj = tasGeno,
        file = "my_gt",
        format = "flapjack"
    )

    fileID1 <- "my_gt.flpjk.geno"
    fileID2 <- "my_gt.flpjk.map"

    expect_true(all(file.exists(fileID1), file.exists(fileID2)))

    file.remove(fileID1)
    file.remove(fileID2)
})


