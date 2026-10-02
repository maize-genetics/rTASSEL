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
    x <- factVctr(c("A", "A", "B", "C"))
    expect_s3_class(x, "fact")
    expect_equal(vctrs::vec_data(x), c("A", "A", "B", "C"))

    x <- covVctr(c(1, 2, 3))
    expect_s3_class(x, "cov")
    expect_equal(vctrs::vec_data(x), c(1, 2, 3))

    x <- dataVctr(c(1, 2, 3))
    expect_s3_class(x, "data")
    expect_equal(vctrs::vec_data(x), c(1, 2, 3))

    x <- taxaVctr(c("A", "A", "B", "C"))
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

test_that("selectTraitsFromFormula calls attributeData, parseFormula, selectTraitsCommon", {
    res <- selectTraitsFromFormula(phFromDf, plant_height ~ PC1)
    expect_true(is(res, "TasselPhenotype"))
    expect_equal(attributeData(res)$trait_id, c("Taxa", "plant_height", "PC1"))
    expect_true(is(javaRefObj(res), "jobjRef"))
})

test_that("selectTraitsCommon keeps the named traits alongside Taxa", {
    out <- selectTraitsCommon(
        attributeData(phFromDf),
        c("plant_height", "PC1"),
        javaRefObj(phFromDf)
    )
    expect_equal(attributeData(out)$trait_id, c("Taxa", "plant_height", "PC1"))
    expect_true(is(javaRefObj(out), "jobjRef"))
})

# === normalizeAttrTypes ===========================================

test_that("normalizeAttrTypes reads a named character vector", {
    spec <- normalizeAttrTypes(c(taxa_id = "taxa", yield = "data"))

    expect_equal(names(spec), c("col_id", "tassel_attr"))
    expect_equal(spec$col_id, c("taxa_id", "yield"))
    expect_equal(spec$tassel_attr, c("taxa", "data"))
})

test_that("normalizeAttrTypes reads both data frame spellings alike", {
    colIdSpec <- tibble::tibble(
        col_id      = c("taxa_id", "yield"),
        tassel_attr = c("taxa", "data")
    )
    traitIdSpec <- tibble::tibble(
        trait_id   = c("taxa_id", "yield"),
        trait_type = c("taxa", "data")
    )

    expect_equal(normalizeAttrTypes(colIdSpec), colIdSpec)
    expect_equal(normalizeAttrTypes(traitIdSpec), colIdSpec)
})

test_that("normalizeAttrTypes ignores the extra columns attributeData() carries", {
    spec <- normalizeAttrTypes(attributeData(phFromDf))

    expect_equal(names(spec), c("col_id", "tassel_attr"))
    expect_true("Taxa" %in% spec$col_id)
})

test_that("normalizeAttrTypes rejects a vector it would have to read by position", {
    expect_error(
        normalizeAttrTypes(c("taxa", "data")),
        "needs a name"
    )
    expect_error(
        normalizeAttrTypes(c(taxa_id = "taxa", "data")),
        "needs a name"
    )
})

test_that("normalizeAttrTypes rejects other forms by name", {
    expect_error(
        normalizeAttrTypes(list(taxa_id = "taxa")),
        "must be a named character vector or a `data.frame`, not <list>"
    )
    expect_error(
        normalizeAttrTypes(data.frame(a = 1)),
        "does not name its columns and attribute types"
    )
})


# === validateAttrTypes ============================================

phSpec <- normalizeAttrTypes(attrDf)

test_that("validateAttrTypes passes a mapping that describes every column", {
    expect_silent(validateAttrTypes(df, phSpec))
})

test_that("validateAttrTypes errors on illegal attribute types", {
    bad <- phSpec
    bad$tassel_attr[2] <- "fake"

    expect_error(
        validateAttrTypes(df, bad),
        "Illegal TASSEL attributes detected: fake"
    )
})

test_that("validateAttrTypes errors when taxa count != 1", {
    zero <- phSpec
    zero$tassel_attr[1] <- "factor"
    expect_error(
        validateAttrTypes(df, zero),
        "Exactly one 'taxa' attribute must be present in `attrTypes`, not 0"
    )

    two <- phSpec
    two$tassel_attr[2] <- "taxa"
    expect_error(
        validateAttrTypes(df, two),
        "Exactly one 'taxa' attribute must be present in `attrTypes`, not 2"
    )
})

test_that("validateAttrTypes errors on a duplicated column ID", {
    dup <- rbind(phSpec, phSpec[2, ])

    expect_error(
        validateAttrTypes(df, dup),
        "describes a column more than once: plant_height"
    )
})

test_that("validateAttrTypes errors on a column it does not describe", {
    expect_error(
        validateAttrTypes(df, phSpec[-2, ]),
        "not described in `attrTypes`: plant_height"
    )
})

test_that("validateAttrTypes errors on an entry with no column", {
    extra <- phSpec
    extra$col_id[2] <- "missing_col"

    expect_error(
        validateAttrTypes(df, extra),
        "not a column of the data frame: missing_col"
    )
})

test_that("validateAttrTypes errors on a data frame with no observations", {
    expect_error(
        validateAttrTypes(df[0, ], phSpec),
        "has no observations"
    )
})


# === coercePhenotypeColumns =======================================

test_that("coercePhenotypeColumns hands Java only doubles and characters", {
    cols <- coercePhenotypeColumns(df, phSpec)

    expect_equal(names(cols), names(df))
    expect_type(cols$taxa_id, "character")
    expect_type(cols$plant_height, "double")
    expect_type(cols$PC1, "double")
    expect_type(cols$cov, "character")
})

test_that("coercePhenotypeColumns keeps the column order of the data frame", {
    reversed <- phSpec[rev(seq_len(nrow(phSpec))), ]

    expect_equal(names(coercePhenotypeColumns(df, reversed)), names(df))
})

test_that("coercePhenotypeColumns widens an integer numeric column", {
    intDf <- tibble::tibble(taxa_id = c("a", "b"), yield = c(1L, 2L))
    cols  <- coercePhenotypeColumns(
        intDf, normalizeAttrTypes(c(taxa_id = "taxa", yield = "data"))
    )

    expect_type(cols$yield, "double")
    expect_equal(cols$yield, c(1, 2))
})

test_that("coercePhenotypeColumns keeps the labels of an R factor", {
    fctDf <- tibble::tibble(taxa_id = c("a", "b"), loc = factor(c("hi", "lo")))
    cols  <- coercePhenotypeColumns(
        fctDf, normalizeAttrTypes(c(taxa_id = "taxa", loc = "factor"))
    )

    expect_equal(cols$loc, c("hi", "lo"))
})

test_that("coercePhenotypeColumns rejects a non-numeric numeric column", {
    charDf <- tibble::tibble(taxa_id = c("a", "b"), yield = c("1", "2"))
    expect_error(
        coercePhenotypeColumns(
            charDf, normalizeAttrTypes(c(taxa_id = "taxa", yield = "data"))
        ),
        "Column 'yield' is marked as 'data' but is <character>"
    )

    lglDf <- tibble::tibble(taxa_id = c("a", "b"), yield = c(TRUE, FALSE))
    expect_error(
        coercePhenotypeColumns(
            lglDf, normalizeAttrTypes(c(taxa_id = "taxa", yield = "covariate"))
        ),
        "Column 'yield' is marked as 'covariate' but is <logical>"
    )
})

test_that("coercePhenotypeColumns rejects a continuous column marked as a factor", {
    expect_error(
        coercePhenotypeColumns(
            df, normalizeAttrTypes(c(
                taxa_id = "taxa", plant_height = "factor", PC1 = "covariate",
                yield = "data", fct = "covariate", cov = "factor"
            ))
        ),
        "Column 'plant_height' is marked as 'factor' but is <double>"
    )
})

test_that("coercePhenotypeColumns rejects a non-character taxa column", {
    numTaxaDf <- tibble::tibble(taxa_id = c(1, 2), yield = c(1, 2))
    expect_error(
        coercePhenotypeColumns(
            numTaxaDf, normalizeAttrTypes(c(taxa_id = "taxa", yield = "data"))
        ),
        "Column 'taxa_id' is marked as 'taxa' but is <double>"
    )
})

test_that("coercePhenotypeColumns rejects missing taxa and factor values", {
    naSpec <- normalizeAttrTypes(c(taxa_id = "taxa", loc = "factor"))

    expect_error(
        coercePhenotypeColumns(
            tibble::tibble(taxa_id = c("a", NA), loc = c("hi", "lo")), naSpec
        ),
        "Column 'taxa_id' is marked as 'taxa' but has missing values"
    )
    expect_error(
        coercePhenotypeColumns(
            tibble::tibble(taxa_id = c("a", "  "), loc = c("hi", "lo")), naSpec
        ),
        "Column 'taxa_id' is marked as 'taxa' but has missing values"
    )
    expect_error(
        coercePhenotypeColumns(
            tibble::tibble(taxa_id = c("a", "b"), loc = c("hi", NA)), naSpec
        ),
        "Column 'loc' is marked as 'factor' but has missing values"
    )
})

test_that("coercePhenotypeColumns carries a missing numeric value through", {
    naDf <- tibble::tibble(taxa_id = c("a", "b"), yield = c(1, NA))
    cols <- coercePhenotypeColumns(
        naDf, normalizeAttrTypes(c(taxa_id = "taxa", yield = "data"))
    )

    expect_equal(cols$yield, c(1, NA_real_))
})

test_that("coercePhenotypeColumns rejects a column type TASSEL has no attribute for", {
    dateDf <- tibble::tibble(
        taxa_id = c("a", "b"),
        sown    = as.Date(c("2020-01-01", "2020-01-02"))
    )

    expect_error(
        coercePhenotypeColumns(
            dateDf, normalizeAttrTypes(c(taxa_id = "taxa", sown = "data"))
        ),
        "Column 'sown' is of an unsupported type <Date>"
    )
})


# === readPhenotypeFromDf ==========================================

test_that("readPhenotypeFromDf builds from a named character vector", {
    ph <- readPhenotype(
        tibble::tibble(taxa_id = c("a", "b"), yield = c(1, 2)),
        attrTypes = c(taxa_id = "taxa", yield = "data")
    )

    expect_s4_class(ph, "TasselPhenotype")
    expect_equal(traitNames(ph), "yield")
})

test_that("readPhenotypeFromDf reads a missing numeric value as TASSEL missing", {
    ph <- readPhenotype(
        tibble::tibble(taxa_id = c("a", "b"), yield = c(1, NA)),
        attrTypes = c(taxa_id = "taxa", yield = "data")
    )
    jPh <- javaRefObj(ph)

    expect_true(jPh$isMissing(1L, jPh$attributeIndexForName("yield")))
    expect_false(jPh$isMissing(0L, jPh$attributeIndexForName("yield")))
})

test_that("readPhenotypeFromDf keeps the trait order of the data frame", {
    ph <- readPhenotype(
        tibble::tibble(taxa_id = "a", zed = 1, alpha = 2),
        attrTypes = c(taxa_id = "taxa", zed = "data", alpha = "data")
    )

    expect_equal(traitNames(ph), c("zed", "alpha"))
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


