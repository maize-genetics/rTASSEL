# === Tests for distance matrix methods =============================

## Preamble - load data ----

### Start logging info
startLogger()

### Load kinship data
kinshipPath <- system.file(
    "extdata",
    "mdp_kinship.txt",
    package = "rTASSEL"
)

### Create R-based kinship object (pairwise)
m <- 10
kinshipR <- diag(1, m, m)
colnames(kinshipR) <- rownames(kinshipR) <- paste0("s_", seq_len(m))


## Error tests ----
test_that("asTasselDistanceMatrix() throws general exceptions", {
    expect_error(
        object = asTasselDistanceMatrix(mtcars),
        regexp = "'m' parameter must be a matrix object"
    )
    expect_error(
        object = asTasselDistanceMatrix(matrix(data = 1, nrow = 10, ncol = 9)),
        regexp = "Matrix object must have equal rows and columns"
    )
    expect_error(
        object = asTasselDistanceMatrix(
            m = matrix(
                data = 1,
                nrow = 10,
                ncol = 10,
                dimnames = list(
                    paste0("r_", 1:10),
                    paste0("c_", 1:10)
                )
            )
        ),
        regexp = "Matrix object must have the same row and column name structure"
    )
})

test_that("readTasselDistanceMatrix() throws general exceptions", {
    expect_error(
        object = readTasselDistanceMatrix("bad/file/path"),
        regexp = "File does not exist."
    )
})


## Equality tests ----
test_that("readTasselDistanceMatrix() method returns correct properties", {
    tasselKinship <- readTasselDistanceMatrix(kinshipPath)

    expect_equal(
        object = class(tasselKinship)[1],
        expected = "TasselDistanceMatrix"
    )
    expect_equal(
        object = head(colnames(tasselKinship)),
        expected = c("33-16", "38-11", "4226", "4722","A188", "A214N")
    )
    expect_equal(
        object = dim(tasselKinship),
        expected = c(277, 277)
    )
    expect_equal(
        object = nrow(tasselKinship),
        expected = 277
    )
    expect_equal(
        object = ncol(tasselKinship),
        expected = 277
    )
})


test_that("asTasselDistanceMatrix() method returns correct properties", {
    tasselKinship <- asTasselDistanceMatrix(kinshipR)
    expect_equal(
        object = class(tasselKinship)[1],
        expected = "TasselDistanceMatrix"
    )
    expect_equal(
        object = head(colnames(tasselKinship)),
        expected = c("s_1", "s_2", "s_3", "s_4", "s_5", "s_6")
    )
    expect_equal(
        object = dim(tasselKinship),
        expected = c(10, 10)
    )
    expect_equal(
        object = nrow(tasselKinship),
        expected = 10
    )
    expect_equal(
        object = ncol(tasselKinship),
        expected = 10
    )
})
















