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

test_that("asTasselDistanceMatrix handles special matrix cases", {
    # Test empty matrix
    empty_mat <- matrix(nrow=0, ncol=0)
    expect_error(asTasselDistanceMatrix(empty_mat), "Matrix object must have equal rows and columns")
    
    # Test matrix with NA values
    m <- 3
    mat_with_na <- matrix(c(1,NA,0.5,NA,1,0.3,0.5,0.3,1), nrow=m)
    rownames(mat_with_na) <- colnames(mat_with_na) <- letters[1:m]
    expect_warning(asTasselDistanceMatrix(mat_with_na))
    
    # Test non-symmetric matrix
    non_sym <- matrix(runif(9), nrow=3)
    rownames(non_sym) <- colnames(non_sym) <- letters[1:3]
    expect_warning(asTasselDistanceMatrix(non_sym))
})

test_that("TasselDistanceMatrix methods work correctly", {
    m <- 4
    test_mat <- matrix(runif(m*m), nrow=m)
    test_mat[lower.tri(test_mat)] <- t(test_mat)[lower.tri(test_mat)]
    diag(test_mat) <- 1
    rownames(test_mat) <- colnames(test_mat) <- paste0("sample_", 1:m)
    
    dist_obj <- asTasselDistanceMatrix(test_mat)
    
    # Test dimensionality methods
    expect_equal(nrow(dist_obj), m)
    expect_equal(ncol(dist_obj), m)
    expect_equal(dim(dist_obj), c(m,m))
    
    # Test row/colnames methods
    expect_equal(rownames(dist_obj), paste0("sample_", 1:m))
    expect_equal(colnames(dist_obj), paste0("sample_", 1:m))
    
    # Test as.matrix conversion
    mat_back <- as.matrix(dist_obj)
    expect_equal(mat_back, test_mat)
})

test_that("readTasselDistanceMatrix handles various file formats", {
    # Test with missing file extension
    expect_error(readTasselDistanceMatrix("nonexistent"), "File does not exist")
    
    # Should handle both square and triangular matrix formats
    expect_true(is(readTasselDistanceMatrix(kinshipPath), "TasselDistanceMatrix"))
})
















