# === Tests for bracket-based filtering of phenotype data ============
#
# The bracket grammar and the verbs are backed by the same machinery, so
# most of what is asserted is that `ph[i, j]` agrees either with the
# equivalent verb or with the equivalent subset of the phenotype's own
# `data.frame`.

## Preamble ----
startLogger()

# 'mdp_phenotype.txt' holds 563 observations of 284 taxa and carries all
# three non-taxa attribute types, with missing values in `EarDia`
ph <- rtObjs$ph_full
ds <- rtObjs$ds_hmp_ph_full

phDf <- as.data.frame(ph)

nObs   <- function(x) nrow(as.data.frame(x))
nTaxa  <- function(x) getGenotypeTable(x)$numberOfTaxa()
nSites <- function(x) getGenotypeTable(x)$numberOfSites()


# /// Selector construction tests ////////////////////////////////////

test_that("traits() creates a valid TraitSelector", {
    sel <- traits("EarHT", "dpoll")
    expect_s4_class(sel, "TraitSelector")
    expect_equal(sel@type, "names")
    expect_equal(sel@ids, c("EarHT", "dpoll"))
    expect_false(sel@negate)
})

test_that("traitsWhere() creates a valid predicate TraitSelector", {
    sel <- traitsWhere(traitType == "covariate")
    expect_s4_class(sel, "TraitSelector")
    expect_equal(sel@type, "predicate")
    expect_false(sel@negate)
})

test_that("! toggles negate on TraitSelector", {
    sel <- traits("EarHT")
    expect_true((!sel)@negate)
    expect_false((!!sel)@negate)
})

test_that("traits() rejects empty input", {
    expect_error(traits(), "At least one trait name")
    expect_error(traits(character(0)), "At least one trait name")
})


# /// Bracket filtering: observations /////////////////////////////////

test_that("[taxa()] keeps every observation of the named taxa", {
    sub <- ph[taxa("33-16", "38-11"), ]

    expect_s4_class(sub, "TasselPhenotype")
    expect_equal(
        as.data.frame(sub),
        as.data.frame(ph |> selectTaxa("33-16", "38-11"))
    )
    expect_equal(nObs(sub), sum(phDf$Taxa %in% c("33-16", "38-11")))
})

test_that("[character] selects taxa by plain character vector", {
    expect_equal(
        as.data.frame(ph[c("33-16", "38-11"), ]),
        as.data.frame(ph[taxa("33-16", "38-11"), ])
    )
})

test_that("[taxa()] skips IDs that name no taxon", {
    expect_equal(
        nObs(ph[taxa("33-16", "not-a-taxon"), ]),
        nObs(ph[taxa("33-16"), ])
    )
})

test_that("[!taxa()] keeps every taxon that was not named", {
    sub <- ph[!taxa("33-16"), ]

    expect_equal(nObs(sub), sum(phDf$Taxa != "33-16"))
    expect_false("33-16" %in% as.data.frame(sub)$Taxa)
})

test_that("[taxaWhere()] tests each observation on its own", {
    expect_equal(
        as.data.frame(ph[taxaWhere(EarHT > 100), ]),
        as.data.frame(ph |> filterTaxa(EarHT > 100))
    )
    expect_equal(
        as.data.frame(ph[taxaWhere(EarHT > 100), ])$Taxa,
        phDf$Taxa[which(phDf$EarHT > 100)]
    )
})

test_that("[taxaWhere()] drops observations whose predicate is missing", {
    expect_true(any(is.na(phDf$EarDia)))
    expect_equal(
        nObs(ph[taxaWhere(EarDia > 0), ]),
        sum(phDf$EarDia > 0, na.rm = TRUE)
    )
})

test_that("[taxaWhere()] exposes the taxa column as taxaId", {
    expect_equal(
        nObs(ph[taxaWhere(startsWith(taxaId, "CML")), ]),
        sum(startsWith(phDf$Taxa, "CML"))
    )
})

test_that("[taxaWhere()] combines criteria with & in one expression", {
    expect_equal(
        as.data.frame(ph[taxaWhere(location == "A" & !is.na(EarDia)), ]),
        as.data.frame(ph |> filterTaxa(location == "A", !is.na(EarDia)))
    )
})

test_that("[!taxaWhere()] keeps the observations the predicate did not", {
    expect_equal(
        as.data.frame(ph[!taxaWhere(location == "A"), ]),
        as.data.frame(ph |> filterTaxa(location != "A"))
    )
})

test_that("[taxaWhere()] does not offer genotype metrics on a phenotype", {
    expect_error(ph[taxaWhere(notMissing >= 0.9), ], "notMissing")
    expect_error(ph[taxaWhere(het <= 0.1), ], "het")
})

test_that("[i] keeps every attribute and its TASSEL type", {
    sub <- ph[taxaWhere(EarHT > 100), ]

    expect_equal(attributeData(sub)$trait_id, attributeData(ph)$trait_id)
    expect_equal(attributeData(sub)$trait_type, attributeData(ph)$trait_type)
})

test_that("[i] rejects a non-logical predicate", {
    expect_error(ph[taxaWhere(EarHT), ], "logical vector")
})

test_that("[i] errors when nothing is left", {
    expect_error(ph[taxaWhere(EarHT > 1e6), ], "No observations match")
    expect_error(ph[taxa("not-a-taxon"), ], "No taxa match")
    expect_error(ph[!taxaWhere(!is.na(Taxa)), ], "No observations match")
})


# /// Bracket filtering: traits //////////////////////////////////////

test_that("[traits()] keeps traits named literally", {
    sub <- ph[, traits("EarHT", "dpoll")]

    expect_s4_class(sub, "TasselPhenotype")
    expect_equal(traitNames(sub), c("EarHT", "dpoll"))
    expect_equal(
        as.data.frame(sub),
        as.data.frame(ph |> selectTraits(EarHT, dpoll))
    )
})

test_that("[character] selects traits by plain character vector", {
    expect_equal(
        as.data.frame(ph[, c("EarHT", "dpoll")]),
        as.data.frame(ph[, traits("EarHT", "dpoll")])
    )
})

test_that("[numeric] selects traits by 1-based position", {
    expect_equal(
        as.data.frame(ph[, 1:3]),
        as.data.frame(ph |> sliceTraits(1:3))
    )
    expect_equal(traitNames(ph[, 1:3]), traitNames(ph)[1:3])
})

test_that("[numeric] rejects positions outside the trait axis", {
    expect_error(ph[, 0], "1-based")
    expect_error(ph[, -1], "1-based")
    expect_error(ph[, length(traitNames(ph)) + 1L], "1-based")
})

test_that("[j] always retains the taxa column", {
    sub <- ph[, traits("EarHT")]

    expect_true("taxa" %in% attributeData(sub)$trait_type)
    expect_equal(colnames(as.data.frame(sub)), c("Taxa", "EarHT"))
})

test_that("[!traits()] keeps every trait that was not named", {
    expect_equal(
        as.data.frame(ph[, !traits("EarDia")]),
        as.data.frame(ph |> selectTraits(-EarDia))
    )
})

test_that("[traitsWhere()] keeps traits by TASSEL attribute type", {
    sub <- ph[, traitsWhere(traitType == "covariate")]

    expect_equal(
        as.data.frame(sub),
        as.data.frame(ph |> filterTraits(traitType == "covariate"))
    )
    expect_equal(
        unique(attributeData(sub)$trait_type),
        c("taxa", "covariate")
    )
})

test_that("[traitsWhere()] keeps traits by their missingness", {
    # `EarDia` is the one trait missing more than 5% of its observations
    keep <- colMeans(!is.na(phDf[traitNames(ph)])) >= 0.95

    expect_equal(
        traitNames(ph[, traitsWhere(notMissing >= 0.95)]),
        traitNames(ph)[keep]
    )
})

test_that("[traitsWhere()] indexes traits from one, skipping taxa", {
    expect_equal(
        traitNames(ph[, traitsWhere(traitIndex <= 3)]),
        traitNames(ph[, 1:3])
    )
})

test_that("[!traitsWhere()] keeps the traits the predicate did not", {
    expect_equal(
        traitNames(ph[, !traitsWhere(traitType == "covariate")]),
        traitNames(ph |> filterTraits(traitType != "covariate"))
    )
})

test_that("[j] rejects a non-logical predicate", {
    expect_error(ph[, traitsWhere(traitId)], "logical vector")
})

test_that("[j] errors when nothing is left", {
    expect_error(ph[, traits("not-a-trait")], "No traits match")
    expect_error(ph[, traitsWhere(notMissing > 1)], "No traits match")
    expect_error(ph[, !traits(traitNames(ph))], "No traits match")
})

test_that("[j] cannot name the taxa column", {
    expect_error(ph[, traits("Taxa")], "No traits match")
})


# /// Combined observations + traits /////////////////////////////////

test_that("bracket filters both axes simultaneously", {
    expect_equal(
        as.data.frame(ph[
            taxaWhere(location == "A" & !is.na(EarDia)),
            traits("EarHT", "dpoll", "EarDia")
        ]),
        as.data.frame(
            ph |>
                filterTaxa(location == "A", !is.na(EarDia)) |>
                selectTraits(EarHT, dpoll, EarDia)
        )
    )
})

test_that("a trait predicate sees the observations the row index left", {
    # Every observation missing `EarDia` is gone by the time the trait
    # index is evaluated, so the same threshold now keeps that trait
    expect_false(
        "EarDia" %in% traitNames(ph[, traitsWhere(notMissing >= 0.95)])
    )
    expect_true(
        "EarDia" %in% traitNames(
            ph[taxaWhere(!is.na(EarDia)), traitsWhere(notMissing >= 0.95)]
        )
    )
})

test_that("bracket indexes are equivalent to chained calls", {
    expect_equal(
        as.data.frame(ph[taxaWhere(EarHT > 100), traits("EarHT", "dpoll")]),
        as.data.frame(ph[taxaWhere(EarHT > 100), ][, traits("EarHT", "dpoll")])
    )
})

test_that("bracket with no indexes returns the input unchanged", {
    expect_equal(as.data.frame(ph[]), as.data.frame(ph))
    expect_equal(attributeData(ph[]), attributeData(ph))
})


# /// Bracket filtering on a genomic dataset /////////////////////////

test_that("bracket routes a phenotype predicate to the observations", {
    sub <- ds[taxaWhere(EarHT > 100), ]

    expect_s4_class(sub, "TasselGenomicDataset")
    expect_equal(
        as.data.frame(sub),
        as.data.frame(ds |> filterTaxa(EarHT > 100))
    )
    # The two halves are re-joined, so the genotype table loses the taxa
    # left without any observations
    expect_lt(nTaxa(sub), nTaxa(ds))
})

test_that("bracket keeps the genotype path for genotype metrics", {
    sub <- ds[taxaWhere(notMissing >= 0.9), ]

    expect_equal(nTaxa(sub), nTaxa(ds |> filterTaxa(notMissing >= 0.9)))
    expect_lt(nTaxa(sub), nTaxa(ds))
})

test_that("bracket mixes genotype and phenotype criteria in one call", {
    expect_equal(
        as.data.frame(ds[taxaWhere(notMissing >= 0.9 & location == "A"), ]),
        as.data.frame(ds |> filterTaxa(notMissing >= 0.9, location == "A"))
    )
})

test_that("a trait index on a dataset leaves the genotype table alone", {
    sub <- ds[, traits("EarHT")]

    expect_s4_class(sub, "TasselGenomicDataset")
    expect_equal(traitNames(sub), "EarHT")
    expect_equal(nSites(sub), nSites(ds))
    expect_equal(nTaxa(sub), nTaxa(ds))
    expect_equal(as.data.frame(sub), as.data.frame(ds |> selectTraits(EarHT)))
})

test_that("a trait predicate on a dataset reads the joined observations", {
    expect_equal(
        traitNames(ds[, traitsWhere(notMissing >= 0.95)]),
        traitNames(ds |> filterTraits(notMissing >= 0.95))
    )
})

test_that("a site index on a dataset still filters sites", {
    expect_equal(nSites(ds[, sites(1:10)]), 10)
    expect_equal(traitNames(ds[, sites(1:10)]), traitNames(ds))
})

test_that("sites and traits are subset by chaining two calls", {
    sub <- ds[, sites(1:10)][, traits("EarHT")]

    expect_equal(nSites(sub), 10)
    expect_equal(traitNames(sub), "EarHT")
})

test_that("both indexes can be given at once on a dataset", {
    sub <- ds[taxa("33-16"), traits("EarHT")]

    expect_equal(nTaxa(sub), 1)
    expect_equal(traitNames(sub), "EarHT")
})

test_that("bracket with no indexes leaves a dataset unchanged", {
    expect_equal(nSites(ds[]), nSites(ds))
    expect_equal(nTaxa(ds[]), nTaxa(ds))
    expect_equal(traitNames(ds[]), traitNames(ds))
})
