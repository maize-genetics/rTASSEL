# === Tests for export functions ====================================
## Preamble - load data ----

### Start logging info
startLogger()

### Shared fixtures (see helper_vars.R)
tasPheno   <- rtObjs$ph_nomiss
tasGeno    <- rtObjs$gt_hmp
tasDataset <- rtObjs$ds_hmp_ph_nomiss


## Error tests ----
test_that("exportGenotypeTable() returns errors", {
    expect_error(
        object = exportGenotypeTable(
            tasObj = mtcars,
            file = "test"
        ),
        regexp = "Unsupported input object"
    )

    expect_error(
        object = exportGenotypeTable(
            tasObj = tasPheno,
            file = "test"
        ),
        regexp = "needs genotype data"
    )

    expect_error(
        object = exportGenotypeTable(
            tasObj = tasGeno,
            file = ""
        ),
        regexp = "File name not specified."
    )

    # TODO - fix Unix/Windows quote bug (Brandon)
    # expect_error(
    #     object = exportGenotypeTable(
    #         tasObj = tasGeno,
    #         file = "my_gt",
    #         format = "csv"
    #     ),
    #     regexp = "'arg' should be one of “vcf”, “hapmap”, “plink”, “flapjack”, “hdf5”"
    # )
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

test_that("exportGenotypeTable() writes a genomic dataset's genotype data.", {
    exportGenotypeTable(
        tasObj = tasDataset,
        file   = "my_ds",
        format = "hapmap"
    )

    fileID <- "my_ds.hmp.txt"

    expect_true(file.exists(fileID))

    file.remove(fileID)
})


## Back-compatibility ----
test_that("exportGenotypeTable() accepts a deprecated TasselGenotypePhenotype", {
    exportGenotypeTable(
        tasObj = rtObjsLegacy$gt_hmp_ph_nomiss,
        file   = "my_legacy_gt",
        format = "hapmap"
    )

    fileID <- "my_legacy_gt.hmp.txt"

    expect_true(file.exists(fileID))

    file.remove(fileID)
})


