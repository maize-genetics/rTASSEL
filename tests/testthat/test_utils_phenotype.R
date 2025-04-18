# Sample attribute and data frames
attrDf <- tibble::tribble(
    ~col_id,       ~tassel_attr,
    "taxa_id",     "taxa",
    "plant_height","data",
    "PC1",         "covariate",
    "yield",       "data",
    "fct",         "covariate",
    "cov",         "factor"
)
df <- tibble::tribble(
    ~taxa_id, ~plant_height, ~PC1,   ~yield, ~fct,  ~cov,
    "line_a",   12.3,           0.5,    2,      0.04, "0.05",
    "line_b",   22.8,          -1.5,    3,      0.04, "0.1"
)

# Example data
phFromDf <- readPhenotype(df, attrDf)

test_that("vctr classes", {
    x <- rTASSEL:::factVctr(c("A", "A", "B", "C"))
    expect_s3_class(x, "fact")
    expect_equal(vctrs::vec_data(x), c("A", "A", "B", "C"))

    x <- rTASSEL:::covVctr(c(1, 2, 3))
    expect_s3_class(x, "cov")
    expect_equal(vctrs::vec_data(x), c(1, 2, 3))

    x <- rTASSEL:::dataVctr(c(1, 2, 3))
    expect_s3_class(x, "data")
    expect_equal(vctrs::vec_data(x), c(1, 2, 3))

    x <- rTASSEL:::taxaVctr(c("A", "A", "B", "C"))
    expect_s3_class(x, "taxa")
    expect_equal(vctrs::vec_data(x), c("A", "A", "B", "C"))
})

test_that("formatPhenotypeDisplay assigns correct vector classes", {
    # Prepare attrDf for formatPhenotypeDisplay
    fmtAttr <- tibble::tibble(
        trait_id   = names(df),
        trait_type = c("taxa", "data", "covariate", "data", "covariate", "factor")
    )
    tbl <- formatPhenotypeDisplay(df, fmtAttr, nCap = 5, nTaxa = nrow(df), jMem = "mem")
    expect_s3_class(tbl, "java_pheno_tbl")
    # Check each column class
    expect_s3_class(tbl$taxa_id,     "taxa")
    expect_s3_class(tbl$plant_height,"data")
    expect_s3_class(tbl$PC1,         "cov")      # covariate
    expect_s3_class(tbl$yield,       "data")
    expect_s3_class(tbl$fct,         "cov")      # covariate
    expect_s3_class(tbl$cov,         "fact")     # factor
})

test_that("selectTraitsCommon warns on missing traits and retains Taxa", {
    expect_warning(
        selectTraitsCommon(
            attributeData(phFromDf),
            c("Taxa", "plant_height", "plant_width"),
            javaRefObj(phFromDf)
        )
    )
})

test_that("selectTraitsFromFormula calls attributeData, parseFormula, selectTraits", {
    res <- selectTraitsFromFormula(phFromDf, plant_height ~ PC1)
    expect_true(is(res, "TasselPhenotype"))
    expect_equal(attributeData(res)$trait_id, c("Taxa", "plant_height", "PC1"))
    expect_true(is(javaRefObj(res), "jobjRef"))
})

test_that("selectTraits delegates to selectTraitsCommon with correct args", {
    out <- selectTraits(phFromDf, c("plant_height", "PC1"))
    expect_equal(attributeData(out)$trait_id, c("Taxa", "plant_height", "PC1"))
    expect_true(is(javaRefObj(out), "jobjRef"))
})

test_that("validateAttrDf errors if not a data.frame", {
    expect_error(
        validateAttrDf(list()),
        "'attrDf' parameter needs to be of type 'data.frame'"
    )
})

test_that("validateAttrDf errors if required columns missing", {
    bad <- data.frame(a = 1)
    expect_error(
        validateAttrDf(bad),
        "Incorrect column IDs used - must be of type 'col_id' and 'tassel_attr'"
    )
})

test_that("validateTasselAttributes errors on illegal attrs", {
    df_bad <- tibble::tibble(tassel_attr = c("data", "fake"))
    expect_error(
        validateTasselAttributes(df_bad, tibble::tibble(tassel_attr = "taxa")),
        "Illegal TASSEL attributes detected: fake"
    )
})

test_that("validateTasselAttributes errors when taxa count !=1", {
    df_good <- tibble::tibble(tassel_attr = c("taxa", "data"))
    # attrDf with zero taxa
    ad_zero <- tibble::tibble(tassel_attr = c("data", "covariate"))
    expect_error(
        validateTasselAttributes(df_good, ad_zero),
        "Exactly one 'taxa' attribute must be present in 'attrDf'"
    )
    # attrDf with two taxa
    ad_two <- tibble::tibble(tassel_attr = c("taxa", "taxa"))
    expect_error(
        validateTasselAttributes(df_good, ad_two),
        "Exactly one 'taxa' attribute must be present in 'attrDf'"
    )
})

wrongColDf <- df
wrongAttr <- attrDf
wrongAttr$col_id[1] <- "missing_col"
test_that("validateColumns errors on missing column", {
    expect_error(
        validateColumns(df, wrongAttr[1, ]),
        "Column 'missing_col' listed in 'attrDf' not found in 'df'"
    )
})

# non-numeric for data type
wrongAttr2 <- attrDf
wrongAttr2$tassel_attr[wrongAttr2$col_id == "taxa_id"] <- "data"
wrongDf2 <- df
# taxa_id is character
test_that("validateColumns errors on non-numeric data column", {
    expect_error(
        validateColumns(wrongDf2, wrongAttr2[1, ]),
        "Column 'taxa_id' is marked as 'data' but is not numeric"
    )
})

test_that("readPhenotypeFromFile aborts for non-existent file", {
    tmp <- tempfile()
    # ensure file does not exist
    unlink(tmp)
    expect_error(
        readPhenotypeFromFile(tmp),
        "The input path is not a valid file"
    )
})


