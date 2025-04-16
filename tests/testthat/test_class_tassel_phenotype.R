library(testthat)
library(rtassel)

# Test data setup
test_that("TasselPhenotype creation works", {
  # Create test data frame
  test_df <- data.frame(
    taxa_id = c("line_a", "line_b"),
    plant_height = c(12.3, 22.8),
    PC1 = c(0.5, -1.5),
    yield = c(2, 3)
  )
  
  # Create attribute metadata
  attr_df <- data.frame(
    col_id = c("taxa_id", "plant_height", "PC1", "yield"),
    tassel_attr = c("taxa", "data", "covariate", "data"),
    stringsAsFactors = FALSE
  )
  
  # Test phenotype object creation
  expect_no_error(phenotype <- readPhenotype(test_df, attr = attr_df))
  
  # Test object structure
  expect_s4_class(phenotype, "TasselPhenotype")
  expect_true(all(c("attrData", "attrSummary", "dispData", "rData", "jRefObj", "jMemAddress", "jClass") %in% slotNames(phenotype)))
})

test_that("TasselPhenotype methods work correctly", {
  # Setup test data
  test_df <- data.frame(
    taxa_id = c("line_a", "line_b"),
    plant_height = c(12.3, 22.8),
    PC1 = c(0.5, -1.5),
    yield = c(2, 3)
  )
  attr_df <- data.frame(
    col_id = c("taxa_id", "plant_height", "PC1", "yield"),
    tassel_attr = c("taxa", "data", "covariate", "data"),
    stringsAsFactors = FALSE
  )
  phenotype <- readPhenotype(test_df, attr = attr_df)
  
  # Test attributeData method
  expect_s3_class(attributeData(phenotype), "data.frame")
  expect_equal(nrow(attributeData(phenotype)), 4)
  expect_equal(ncol(attributeData(phenotype)), 2)
  
  # Test javaRefObj method
  expect_true(inherits(javaRefObj(phenotype), "jobjRef"))
})

test_that("TasselPhenotype validation works", {
  # Test invalid attribute metadata
  expect_error(
    readPhenotype(data.frame(x = 1), attr = NULL),
    "Phenotype objects evaluated from 'data.frame' need attribute metadata"
  )
  
  # Test invalid input type
  expect_error(
    readPhenotype(list(x = 1)),
    "Unsupported input type for 'x'"
  )
})

test_that("TasselPhenotype display methods work", {
  # Setup minimal test data
  test_df <- data.frame(
    taxa_id = "line_a",
    trait = 1
  )
  attr_df <- data.frame(
    col_id = c("taxa_id", "trait"),
    tassel_attr = c("taxa", "data"),
    stringsAsFactors = FALSE
  )
  phenotype <- readPhenotype(test_df, attr = attr_df)
  
  # Test that show method works without error
  expect_no_error(show(phenotype))
  
  # Test that the display data has correct attributes
  expect_true(all(c("nTaxa", "nTraits", "nDfRow", "nCap", "jMem") %in% names(attributes(phenotype@dispData))))
})