# === Tests for genomic dataset display utilities ====================

## Preamble ----
startLogger()

ds     <- rtObjs$ds_hmp_ph_nomiss
gt     <- rtObjs$gt_hmp
nSites <- 3093L

plainText <- function(x) cli::ansi_strip(paste(x, collapse = "\n"))


# /// genoVctr() /////////////////////////////////////////////////////

test_that("genoVctr() carries its display metadata", {
    gv <- genoVctr(
        x            = list(c("A", "C"), c("G", "N")),
        nSites       = 100L,
        nSitesShown  = 2L,
        minorAlleles = c("A", "G")
    )

    expect_s3_class(gv, "geno")
    expect_length(gv, 2)
    expect_equal(attr(gv, "nSites"), 100L)
    expect_equal(attr(gv, "nSitesShown"), 2L)
    expect_equal(attr(gv, "minorAlleles"), c("A", "G"))
    expect_false(attr(gv, "numeric"))
})

test_that("genoVctr() metadata survives subsetting", {
    gv <- genoVctr(
        x           = list(c("A", "C"), c("G", "N")),
        nSites      = 100L,
        nSitesShown = 2L
    )
    gvSub <- gv[1]

    expect_length(gvSub, 1)
    expect_equal(attr(gvSub, "nSites"), 100L)
})

test_that("genoVctr() abbreviates as <geno> in a tibble", {
    tbl <- tibble::tibble(
        gt = genoVctr(list("A"), nSites = 1L, nSitesShown = 1L)
    )

    expect_equal(vctrs::vec_ptype_abbr(tbl$gt), "geno")
})


# /// buildGenotypeColumn() //////////////////////////////////////////

test_that("buildGenotypeColumn() collects allele calls per observation", {
    gv <- buildGenotypeColumn(
        jGp   = ds@jRefObj,
        jGt   = genotype(ds)@jRefObj,
        nRows = 4
    )

    expect_s3_class(gv, "geno")
    expect_length(gv, 4)
    expect_equal(attr(gv, "nSites"), nSites)
    expect_equal(attr(gv, "nSitesShown"), 5L)
    expect_length(attr(gv, "minorAlleles"), 5)
    expect_equal(vctrs::vec_data(gv)[[1]], c("C", "C", "G", "T", "G"))
})

test_that("buildGenotypeColumn() honors the site limit", {
    gv <- buildGenotypeColumn(
        jGp    = ds@jRefObj,
        jGt    = genotype(ds)@jRefObj,
        nRows  = 2,
        nSites = 3
    )

    expect_equal(attr(gv, "nSitesShown"), 3L)
    expect_length(vctrs::vec_data(gv)[[1]], 3)
})

test_that("buildGenotypeColumn() never asks for more sites than exist", {
    dsSmall <- ds[, sites(1:2)]

    gv <- buildGenotypeColumn(
        jGp   = dsSmall@jRefObj,
        jGt   = genotype(dsSmall)@jRefObj,
        nRows = 2
    )

    expect_equal(attr(gv, "nSites"), 2L)
    expect_equal(attr(gv, "nSitesShown"), 2L)
})

test_that("buildGenotypeColumn() handles a single observation", {
    gv <- buildGenotypeColumn(
        jGp   = ds@jRefObj,
        jGt   = genotype(ds)@jRefObj,
        nRows = 1
    )

    expect_length(gv, 1)
    expect_length(vctrs::vec_data(gv)[[1]], 5)
})

test_that("buildGenotypeColumn() reports missing calls for union joins", {
    dsUnion <- readGenomicDataset(gt, rtObjs$ph_nomiss, join = "union")

    gv    <- buildGenotypeColumn(dsUnion@jRefObj, genotype(dsUnion)@jRefObj, 20)
    calls <- unlist(vctrs::vec_data(gv))

    # Taxa that only the phenotype knows about have no genotype calls
    expect_true("N" %in% calls)
})

test_that("buildGenotypeColumn() collects reference probabilities", {
    dsNum <- readGenomicDataset(
        rtMatrices$num_gt_sm,
        data.frame(Taxa = rownames(rtMatrices$num_gt_sm), yield = c(1, 2, 3)),
        attr = data.frame(
            col_id      = c("Taxa", "yield"),
            tassel_attr = c("taxa", "data")
        )
    )

    gv <- buildGenotypeColumn(dsNum@jRefObj, genotype(dsNum)@jRefObj, 3)

    expect_true(attr(gv, "numeric"))
    expect_equal(attr(gv, "minorAlleles"), character())
    expect_type(vctrs::vec_data(gv)[[1]], "double")
    expect_true(all(unlist(vctrs::vec_data(gv)) >= 0))
})


# /// pillar_shaft.geno() ////////////////////////////////////////////

test_that("pillar_shaft.geno() renders allele calls with an ellipsis", {
    gv <- genoVctr(
        x            = list(c("A", "C")),
        nSites       = 100L,
        nSitesShown  = 2L,
        minorAlleles = c("A", "C")
    )

    out <- plainText(format(pillar::pillar_shaft(gv), width = 60))

    expect_match(out, "A")
    expect_match(out, "C")
    expect_match(out, cli::symbol$ellipsis, fixed = TRUE)
})

test_that("pillar_shaft.geno() drops the ellipsis when nothing is hidden", {
    gv <- genoVctr(
        x            = list(c("A", "C")),
        nSites       = 2L,
        nSitesShown  = 2L,
        minorAlleles = c("A", "C")
    )

    out <- plainText(format(pillar::pillar_shaft(gv), width = 60))

    expect_false(grepl(cli::symbol$ellipsis, out, fixed = TRUE))
})

test_that("pillar_shaft.geno() renders reference probabilities", {
    gv <- genoVctr(
        x           = list(c(0.25, 0.75)),
        nSites      = 2L,
        nSitesShown = 2L,
        numeric     = TRUE
    )

    out <- plainText(format(pillar::pillar_shaft(gv), width = 60))

    expect_match(out, "0.250")
    expect_match(out, "0.750")
})


# /// formatColumnTypeSummary() //////////////////////////////////////

test_that("formatColumnTypeSummary() orders types canonically", {
    expect_equal(
        formatColumnTypeSummary(
            list(covariate = 3, data = 3, factor = 1, taxa = 1)
        ),
        "taxa: 1, factor: 1, data: 3, covariate: 3, genotype: 1"
    )
})

test_that("formatColumnTypeSummary() can omit the genotype column", {
    expect_equal(
        formatColumnTypeSummary(list(data = 2, taxa = 1), hasGenotype = FALSE),
        "taxa: 1, data: 2"
    )
    expect_equal(formatColumnTypeSummary(list(), hasGenotype = FALSE), "no columns")
})

test_that("formatColumnTypeSummary() keeps unrecognized types", {
    expect_match(
        formatColumnTypeSummary(list(taxa = 1, mystery = 2)),
        "mystery: 2$"
    )
})


# /// formatGenomicDatasetDisplay() //////////////////////////////////

test_that("formatGenomicDatasetDisplay() builds a java_geno_pheno_tbl", {
    dispData <- formatGenomicDatasetDisplay(ds)

    expect_s3_class(dispData, "java_geno_pheno_tbl")
    expect_s3_class(dispData, "tbl_df")
    expect_equal(
        colnames(dispData),
        c("Taxa", "EarHT", "dpoll", "EarDia", "Genotype")
    )
    expect_equal(nrow(dispData), 10L)
})

test_that("formatGenomicDatasetDisplay() sets the display attributes", {
    dispData <- formatGenomicDatasetDisplay(ds)

    expect_true(all(
        c("nTaxa", "nSites", "nSitesShown", "nCap", "nDfRow", "colTypes",
          "jMem", "pillar_focus") %in% names(attributes(dispData))
    ))
    expect_equal(attr(dispData, "nTaxa"), 278L)
    expect_equal(attr(dispData, "nSites"), nSites)
    expect_equal(attr(dispData, "pillar_focus"), "Genotype")
    expect_equal(attr(dispData, "jMem"), ds@jMemAddress)
})

test_that("formatGenomicDatasetDisplay() reuses the phenotype columns", {
    dispData <- formatGenomicDatasetDisplay(ds)

    expect_s3_class(dispData$Taxa, "taxa")
    expect_s3_class(dispData$EarHT, "data")
    expect_s3_class(dispData$Genotype, "geno")
})

test_that("formatGenomicDatasetDisplay() respects a site limit", {
    dispData <- formatGenomicDatasetDisplay(ds, nSites = 2)

    expect_equal(attr(dispData, "nSitesShown"), 2L)
})


# /// tbl_format_header() / tbl_format_footer() //////////////////////

test_that("tbl_format_header() reports the dataset dimensions", {
    hdr <- plainText(pillar::tbl_format_header(formatGenomicDatasetDisplay(ds)))

    expect_match(hdr, "TasselGenomicDataset", fixed = TRUE)
    expect_match(hdr, "278 taxa")
    expect_match(hdr, "3093 sites")
})

test_that("tbl_format_footer() reports rows, sites, types, and memory", {
    dispData <- formatGenomicDatasetDisplay(ds)
    setup    <- pillar::tbl_format_setup(dispData)
    ftr      <- plainText(pillar::tbl_format_footer(dispData, setup))

    expect_match(ftr, "showing the first 10 of 278 rows")
    expect_match(ftr, "Genotype: showing the first 5 of 3093 sites")
    expect_match(ftr, "Column types: taxa: 1, data: 3, genotype: 1")
    expect_match(ftr, ds@jMemAddress)
})
