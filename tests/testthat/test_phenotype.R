# === Tests for phenotype functions =================================

## Preamble - load data ----

### Start logging info
startLogger()
library(rJava)

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
test_that("readPhenotypeFromPath() throws general exceptions.", {
    badPath <- "does/not/exist"
    expandedPath <- normalizePath(badPath, mustWork = FALSE)
    message <- paste0("Cannot open file ", expandedPath, ": No such file or directory")

    expect_error(
        object = readPhenotypeFromPath(badPath),
        regexp = message,
        fixed = TRUE
    )
})

test_that("readPhenotypeFromPath() expands ~ in path before passing to Java", {
    skip_on_os("windows")

    # Point HOME at the directory containing the bundled test file so that
    # "~/<basename>" resolves to a real file.
    withr::local_envvar(c(HOME = dirname(phenoPathFast)))
    tildePath <- file.path("~", basename(phenoPathFast))

    expect_match(tildePath, "^~/")

    expect_message(
        ph <- readPhenotypeFromPath(tildePath),
        "deprecated"
    )
    expect_s4_class(ph, "TasselGenotypePhenotype")
})

test_that("readPhenotypeFromPath() reports expanded path on missing ~ file", {
    skip_on_os("windows")

    withr::local_envvar(c(HOME = tempdir()))
    badPath <- "~/__rtassel_definitely_missing__.txt"
    expandedPath <- normalizePath(badPath, mustWork = FALSE)

    expect_error(
        readPhenotypeFromPath(badPath),
        regexp = paste0("Cannot open file ", expandedPath, ": No such file or directory"),
        fixed = TRUE
    )
})

test_that("readGenotypePhenotype() expands ~ in phenotype path", {
    skip_on_os("windows")

    # Both bundled test files live in the same extdata directory, so a single
    # HOME override lets us refer to either via "~/<basename>".
    extDir <- dirname(genoPathHMP)
    expect_equal(extDir, dirname(phenoPathFast))

    withr::local_envvar(c(HOME = extDir))
    phenoTilde <- file.path("~", basename(phenoPathFast))

    expect_match(phenoTilde, "^~/")

    # Use the absolute geno path; the phenotype path goes through the previously
    # buggy PhenotypeBuilder$fromFile() branch that did not normalize the path.
    obj <- readGenotypePhenotype(
        genoPathOrObj    = genoPathHMP,
        phenoPathDFOrObj = phenoTilde
    )
    expect_s4_class(obj, "TasselGenotypePhenotype")
})

test_that("readPhenotypeFromDataFrame() throws general exceptions.", {

    phenoDF <- getPhenotypeDF(tasGenoPheno)
    phenoDF <- as.data.frame(phenoDF)

    errorMsg <- paste0(
        "Parameter `attributeTypes` contains incorrect attributes.\n",
        "Please select from the following:\n",
        "  taxa\n",
        "  factor\n",
        "  data\n",
        "  covariate\n"
    )

    expect_error(
        object = readPhenotypeFromDataFrame(
            phenotypeDF    = phenoDF,
            taxaID         = "Taxa",
            attributeTypes = "wrong"
        ),
        regexp = errorMsg
    )
})

test_that("getPhenotypeDF() throws general exceptions.", {
    expect_error(
        object = getPhenotypeDF(mtcars),
        regexp = "`tasObj` must be of class `TasselGenotypePhenotype`"
    )
    expect_error(
        object = getPhenotypeDF(tasGeno),
        regexp = "TASSEL phenotype object not found"
    )
})

# test_that("makeAttributeData validates attribute types correctly", {
#     phenoDF <- getPhenotypeDF(tasGenoPheno)
#     attrDF <- data.frame(
#         col_id = c("Taxa", "invalid"),
#         tassel_attr = c("taxa", "unknown")
#     )
#     expect_error(
#         makeAttributeData(tasGenoPheno@jPhenotypeTable, phenoDF, attrDF),
#         "Incorrect TASSEL attributes found"
#     )
#
#     # Test with valid attributes
#     validAttrDF <- data.frame(
#         col_id = c("Taxa", "EarHT"),
#         tassel_attr = c("taxa", "data")
#     )
#     result <- makeAttributeData(tasGenoPheno@jPhenotypeTable, phenoDF, validAttrDF)
#     expect_equal(nrow(result), 2)
#     expect_equal(ncol(result), 2)
# })

# test_that("validateAttributeDataFrame catches invalid inputs", {
#     # Missing required columns
#     badDF <- data.frame(col_id = "Taxa")
#     expect_error(
#         validateAttributeDataFrame(badDF),
#         "Incorrect column IDs used"
#     )
#
#     # Invalid data type
#     expect_error(
#         validateAttributeDataFrame("not a dataframe"),
#         "'attrDf' parameter needs to be of type 'data.frame'"
#     )
#
#     # Valid input should not error
#     validDF <- data.frame(
#         col_id = "Taxa",
#         tassel_attr = "taxa"
#     )
#     expect_silent(validateAttributeDataFrame(validDF))
# })

# test_that("validateTasselAttributes handles edge cases", {
#     # Test with missing taxa attribute
#     attrDF <- data.frame(
#         col_id = c("val1", "val2"),
#         tassel_attr = c("data", "data")
#     )
#     expect_error(
#         validateTasselAttributes(attrDF),
#         "No 'taxa' attribute found"
#     )
#
#     # Test with invalid attribute type
#     attrDF <- data.frame(
#         col_id = c("Taxa", "val1"),
#         tassel_attr = c("taxa", "invalid")
#     )
#     expect_error(
#         validateTasselAttributes(attrDF),
#         "Incorrect TASSEL attributes found"
#     )
#
#     # Valid case should work
#     validDF <- data.frame(
#         col_id = c("Taxa", "val1", "val2"),
#         tassel_attr = c("taxa", "data", "covariate")
#     )
#     expect_silent(validateTasselAttributes(validDF))
# })



## Return tests ----
test_that("readPhenotypeFromDataFrame() returns correct data.", {

    phenoDF <- getPhenotypeDF(tasGenoPheno)
    phenoDF <- as.data.frame(phenoDF)

    phenoDFTas <- readPhenotypeFromDataFrame(
        phenotypeDF = phenoDF,
        taxaID = "Taxa"
    )

    expect_true(class(phenoDFTas) == "TasselGenotypePhenotype")
})

test_that("getPhenotypeDF () returns correct data.", {

    phenoDF <- getPhenotypeDF(tasGenoPheno)

    expect_true(class(phenoDF) == "data.frame")

    expect_equal(
        object   = names(phenoDF),
        expected = c("Taxa", "EarHT", "dpoll", "EarDia")
    )
})








