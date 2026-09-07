# === Tests for the filter verbs on phenotype data ==================
#
# The taxa verbs address observations rather than sites here, and the
# trait verbs are the second axis, so most of what is asserted is that a
# verb agrees with the equivalent subset of the phenotype's own
# `data.frame`.

## Preamble ----
startLogger()

# `lifecycle::expect_deprecated()` matches on condition class, which
# needs the third edition of testthat
local_edition(3)

# 'mdp_phenotype.txt' holds 563 observations of 284 taxa and carries all
# three non-taxa attribute types, with missing values in `EarDia`
ph <- rtObjs$ph_full
ds <- rtObjs$ds_hmp_ph_full

phDf <- as.data.frame(ph)

nObs    <- function(x) nrow(as.data.frame(x))
nTaxa   <- function(x) getGenotypeTable(x)$numberOfTaxa()
nSites  <- function(x) getGenotypeTable(x)$numberOfSites()
obsTaxa <- function(x) unique(as.data.frame(x)$Taxa)


# /// filterTaxa() ///////////////////////////////////////////////////

test_that("filterTaxa() keeps the observations a data frame subset would", {
    verb <- ph |> filterTaxa(EarHT > 100)

    expect_s4_class(verb, "TasselPhenotype")
    expect_equal(
        as.data.frame(verb)$Taxa,
        phDf$Taxa[which(phDf$EarHT > 100)]
    )
})

test_that("filterTaxa() drops observations whose predicate is missing", {
    # `EarDia` is the one trait with missing values, so a predicate over
    # it is NA rather than FALSE for those observations
    expect_true(any(is.na(phDf$EarDia)))
    expect_equal(
        nObs(ph |> filterTaxa(EarDia > 0)),
        sum(phDf$EarDia > 0, na.rm = TRUE)
    )
})

test_that("filterTaxa() combines several predicates with &", {
    expect_equal(
        nObs(ph |> filterTaxa(location == "A", !is.na(EarDia))),
        sum(phDf$location == "A" & !is.na(phDf$EarDia))
    )
})

test_that("filterTaxa() exposes the taxa column as taxaId", {
    expect_equal(
        nObs(ph |> filterTaxa(startsWith(taxaId, "CML"))),
        sum(startsWith(phDf$Taxa, "CML"))
    )
})

test_that("filterTaxa() keeps every attribute and its TASSEL type", {
    verb <- ph |> filterTaxa(EarHT > 100)

    expect_equal(attributeData(verb)$trait_id, attributeData(ph)$trait_id)
    expect_equal(attributeData(verb)$trait_type, attributeData(ph)$trait_type)
})

test_that("filterTaxa() rejects a non-logical predicate", {
    expect_error(ph |> filterTaxa(EarHT), "logical vector")
})

test_that("filterTaxa() errors when nothing is left", {
    expect_error(ph |> filterTaxa(EarHT > 1e6), "No observations match")
})

test_that("filterTaxa() does not offer genotype metrics on a phenotype", {
    # `notMissing` and `het` describe genotype calls, so a phenotype on
    # its own has nothing to evaluate them against
    expect_error(ph |> filterTaxa(notMissing >= 0.9), "notMissing")
})


# /// selectTaxa() ///////////////////////////////////////////////////

test_that("selectTaxa() keeps every observation of a selected taxon", {
    verb <- ph |> selectTaxa("33-16", "38-11")

    expect_s4_class(verb, "TasselPhenotype")
    expect_setequal(obsTaxa(verb), c("33-16", "38-11"))
    expect_equal(
        nObs(verb),
        sum(phDf$Taxa %in% c("33-16", "38-11"))
    )
})

test_that("selectTaxa() accepts tidyselect helpers", {
    prefixed <- ph |> selectTaxa(starts_with("CML"))
    expect_true(all(startsWith(obsTaxa(prefixed), "CML")))
    expect_equal(
        sort(obsTaxa(prefixed)),
        sort(unique(phDf$Taxa[startsWith(phDf$Taxa, "CML")]))
    )

    expect_setequal(
        obsTaxa(ph |> selectTaxa(matches("^CML"))),
        obsTaxa(prefixed)
    )
})

test_that("selectTaxa() drops with a leading minus", {
    dropped <- c("33-16", "38-11")
    expect_setequal(
        obsTaxa(ph |> selectTaxa(-any_of(dropped))),
        setdiff(unique(phDf$Taxa), dropped)
    )
})

test_that("selectTaxa() errors when nothing is selected", {
    expect_error(ph |> selectTaxa(any_of("NOT_A_TAXON")), "No taxa match")
})


# /// sliceTaxa() ////////////////////////////////////////////////////

test_that("sliceTaxa() keeps taxa by position", {
    expect_setequal(
        obsTaxa(ph |> sliceTaxa(1:10)),
        taxaList(ph)[1:10]
    )
})

test_that("sliceTaxa() drops on negative positions", {
    expect_setequal(
        obsTaxa(ph |> sliceTaxa(-(1:10))),
        taxaList(ph)[-(1:10)]
    )
})

test_that("sliceTaxa() ignores positions past the end", {
    total <- length(taxaList(ph))
    expect_equal(nObs(ph |> sliceTaxa(-(total + 1))), nrow(phDf))
})


# /// selectTraits() /////////////////////////////////////////////////

test_that("selectTraits() keeps traits named literally", {
    verb <- ph |> selectTraits(EarHT, dpoll)

    expect_s4_class(verb, "TasselPhenotype")
    expect_equal(traitNames(verb), c("EarHT", "dpoll"))
})

test_that("selectTraits() always retains the taxa column", {
    verb <- ph |> selectTraits(EarHT)

    expect_equal(colnames(as.data.frame(verb)), c("Taxa", "EarHT"))
    expect_true("taxa" %in% attributeData(verb)$trait_type)
})

test_that("selectTraits() accepts tidyselect helpers", {
    traits <- c("EarHT", "dpoll")

    expect_equal(traitNames(ph |> selectTraits(all_of(traits))), traits)
    expect_equal(
        traitNames(ph |> selectTraits(any_of(c(traits, "NOT_A_TRAIT")))),
        traits
    )
    expect_equal(
        traitNames(ph |> selectTraits(starts_with("Q"))),
        c("Q1", "Q2", "Q3")
    )
})

test_that("selectTraits() selects on trait values as well as names", {
    # The traits are handed to tidyselect as columns, so `where()` sees
    # the values behind each name
    expect_equal(
        traitNames(ph |> selectTraits(where(is.numeric))),
        setdiff(traitNames(ph), "location")
    )
})

test_that("selectTraits() drops with a leading minus", {
    expect_equal(
        traitNames(ph |> selectTraits(-EarDia)),
        setdiff(traitNames(ph), "EarDia")
    )
})

test_that("selectTraits() errors on an absent all_of() trait", {
    expect_error(ph |> selectTraits(all_of("NOT_A_TRAIT")))
})

test_that("selectTraits() errors when nothing is selected", {
    expect_error(
        ph |> selectTraits(any_of("NOT_A_TRAIT")),
        "No traits match"
    )
})


# /// filterTraits() /////////////////////////////////////////////////

test_that("filterTraits() keeps traits by TASSEL attribute type", {
    verb <- ph |> filterTraits(traitType == "covariate")

    expect_s4_class(verb, "TasselPhenotype")
    expect_equal(traitNames(verb), c("Q1", "Q2", "Q3"))
    expect_equal(
        attributeData(verb)$trait_type,
        c("taxa", rep("covariate", 3))
    )
})

test_that("filterTraits() keeps traits by their missingness", {
    notMissing <- vapply(
        traitNames(ph),
        function(id) mean(!is.na(phDf[[id]])),
        numeric(1)
    )
    expect_true(any(notMissing < 0.95))

    expect_equal(
        traitNames(ph |> filterTraits(notMissing >= 0.95)),
        names(notMissing)[notMissing >= 0.95]
    )
})

test_that("filterTraits() combines several predicates with &", {
    # Of the three `data` traits, only `EarDia` misses more than 5% of
    # its observations
    expect_equal(
        traitNames(ph |> filterTraits(traitType == "data", notMissing >= 0.95)),
        c("EarHT", "dpoll")
    )
})

test_that("filterTraits() indexes traits from one, skipping taxa", {
    expect_equal(
        traitNames(ph |> filterTraits(traitIndex <= 2)),
        traitNames(ph)[1:2]
    )
})

test_that("filterTraits() rejects a non-logical predicate", {
    expect_error(ph |> filterTraits(traitId), "logical vector")
})

test_that("filterTraits() errors when nothing is left", {
    expect_error(
        ph |> filterTraits(traitType == "NOT_A_TYPE"),
        "No traits match"
    )
})


# /// sliceTraits() //////////////////////////////////////////////////

test_that("sliceTraits() keeps traits by position", {
    expect_equal(
        traitNames(ph |> sliceTraits(1:3)),
        traitNames(ph)[1:3]
    )
    expect_equal(
        traitNames(ph |> sliceTraits(c(1, 4))),
        traitNames(ph)[c(1, 4)]
    )
})

test_that("sliceTraits() drops on negative positions", {
    expect_equal(
        traitNames(ph |> sliceTraits(-1)),
        traitNames(ph)[-1]
    )
})

test_that("sliceTraits() ignores zeros and positions past the end", {
    total <- length(traitNames(ph))

    expect_equal(
        traitNames(ph |> sliceTraits(c(0, 1, total + 100))),
        traitNames(ph)[1]
    )
    expect_equal(
        traitNames(ph |> sliceTraits(-(total + 1))),
        traitNames(ph)
    )
})

test_that("sliceTraits() cannot mix signs", {
    expect_error(ph |> sliceTraits(c(1, -2)), "positive and negative")
})

test_that("sliceTraits() rejects non-numeric positions", {
    expect_error(ph |> sliceTraits("EarHT"), "must be numeric")
})


# /// No-ops /////////////////////////////////////////////////////////

test_that("phenotype verbs called with no arguments return the input", {
    expect_identical(filterTaxa(ph), ph)
    expect_identical(selectTaxa(ph), ph)
    expect_identical(sliceTaxa(ph), ph)
    expect_identical(filterTraits(ph), ph)
    expect_identical(selectTraits(ph), ph)
    expect_identical(sliceTraits(ph), ph)
})


# /// Genomic datasets ///////////////////////////////////////////////

test_that("trait verbs on a dataset leave the genotype table alone", {
    verb <- ds |> selectTraits(EarHT)

    expect_s4_class(verb, "TasselGenomicDataset")
    expect_equal(traitNames(verb), "EarHT")
    expect_equal(nTaxa(verb), nTaxa(ds))
    expect_equal(nSites(verb), nSites(ds))
    expect_equal(nObs(verb), nObs(ds))

    expect_equal(
        traitNames(ds |> filterTraits(traitType == "data")),
        traitNames(ds |> selectTraits(where(is.numeric), -starts_with("Q")))
    )
})

test_that("filterTaxa() on a dataset trims the genotype to surviving taxa", {
    verb <- ds |> filterTaxa(EarHT > 100)

    expect_s4_class(verb, "TasselGenomicDataset")
    expect_setequal(taxaList(verb), obsTaxa(verb))
    expect_true(nTaxa(verb) < nTaxa(ds))
    expect_equal(nSites(verb), nSites(ds))
})

test_that("filterTaxa() on a dataset keeps the genotype path for genotype metrics", {
    # A predicate naming no phenotype column must still reach the
    # genotype table, and so agree with the bracket form
    verb      <- ds |> filterTaxa(notMissing >= 0.9)
    bracketed <- ds[taxaWhere(notMissing >= 0.9), ]

    expect_equal(taxaList(verb), taxaList(bracketed))
})

test_that("filterTaxa() on a dataset mixes genotype and phenotype criteria", {
    verb <- ds |> filterTaxa(notMissing >= 0.9, location == "A")

    expect_s4_class(verb, "TasselGenomicDataset")
    expect_true(all(as.data.frame(verb)$location == "A"))
    expect_true(all(taxaList(verb) %in% taxaList(ds |> filterTaxa(notMissing >= 0.9))))
})

test_that("taxa verbs that select by ID or position stay on the genotype", {
    expect_equal(
        taxaList(ds |> selectTaxa(starts_with("CML"))),
        taxaList(ds[taxaWhere(startsWith(taxaId, "CML")), ])
    )
    expect_equal(taxaList(ds |> sliceTaxa(1:10)), taxaList(ds)[1:10])
})


# /// Input handling /////////////////////////////////////////////////

test_that("trait verbs reject input without phenotype data", {
    expect_error(
        filterTraits(rtObjs$gt_hmp, notMissing >= 0.5),
        "needs phenotype data"
    )
    expect_error(selectTraits(rtObjs$gt_hmp, EarHT), "needs phenotype data")
    expect_error(sliceTraits(rtObjs$gt_hmp, 1), "needs phenotype data")
    expect_error(filterTraits("not an rTASSEL object", 1), "Unsupported")
})

test_that("verbs warn but work for phenotype-only TasselGenotypePhenotype input", {
    legacy <- rtObjsLegacy$ph_nomiss

    lifecycle::expect_deprecated(
        rows <- filterTaxa(legacy, EarHT > 100),
        "TasselGenotypePhenotype"
    )
    lifecycle::expect_deprecated(
        traits <- selectTraits(legacy, EarHT),
        "TasselGenotypePhenotype"
    )

    expect_s4_class(rows, "TasselGenotypePhenotype")
    expect_s4_class(traits, "TasselGenotypePhenotype")

    # Read through the Java tables rather than `getPhenotypeDF()`, which
    # is itself deprecated and would warn inside this block
    jRows   <- getPhenotypeTable(rows)
    jTraits <- getPhenotypeTable(traits)

    expect_true(jRows$numberOfObservations() > 0)
    expect_equal(jTraits$getTableColumnNames(), c("Taxa", "EarHT"))
})
