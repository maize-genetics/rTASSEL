## ----setup, include=FALSE-----------------------------------------------------
library(rTASSEL)

knitr::opts_chunk$set(
    fig.path='figure/graphics-',
    cache.path='cache/graphics-',
    fig.align='center',
    external=TRUE,
    echo=TRUE,
    warning=FALSE
    # fig.pos="H"
)




## ----eval=FALSE---------------------------------------------------------------
# myPheno[taxaWhere(EarHT > 100), ]
# myPheno |> filterTaxa(EarHT > 100)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
phenoPath <- system.file("extdata", "mdp_phenotype.txt", package = "rTASSEL")

myPheno <- readPhenotype(phenoPath)
myPheno


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
attributeData(myPheno)

traitNames(myPheno)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
length(taxaList(myPheno))

nrow(as.data.frame(myPheno))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[taxaWhere(EarHT > 100), ]

myPheno |> filterTaxa(EarHT > 100)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[, traits("EarHT")]

myPheno |> selectTraits(EarHT)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[taxaWhere(EarHT > 100), ]

myPheno |> filterTaxa(EarHT > 100)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[taxaWhere(EarHT >= 60 & EarHT <= 80), ]

myPheno |> filterTaxa(EarHT >= 60, EarHT <= 80)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTaxa(EarHT > 100 | dpoll > 80)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTaxa(EarHT > mean(EarHT, na.rm = TRUE))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTaxa(Q3 > 0.9)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[taxaWhere(location == "A"), ]

myPheno |> filterTaxa(location == "A")


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTaxa(EarDia > 35)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[taxaWhere(!is.na(EarDia)), ]

myPheno |> filterTaxa(!is.na(EarDia))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTaxa(!is.na(EarHT), !is.na(dpoll), !is.na(EarDia))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[taxa("33-16", "38-11"), ]

myPheno |> selectTaxa("33-16", "38-11")


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myTaxa <- c("33-16", "38-11", "4226")

myPheno |> selectTaxa(all_of(myTaxa))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> selectTaxa(starts_with("CML"))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[taxaWhere(taxaId %in% myTaxa), ]

myPheno |> filterTaxa(taxaId %in% myTaxa)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> sliceTaxa(1:3)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> sliceTaxa(-(1:3))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[taxaWhere(startsWith(taxaId, "CML")), ]

myPheno |> filterTaxa(startsWith(taxaId, "CML"))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTaxa(startsWith(taxaId, "CML"), location == "A", !is.na(EarDia))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
tallEars <- myPheno |> filterTaxa(EarHT > 100)

length(taxaList(tallEars))

nrow(as.data.frame(tallEars))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> selectTaxa(all_of(taxaList(tallEars)))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[, traits("EarHT", "dpoll")]

myPheno |> selectTraits(EarHT, dpoll)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myTraits <- c("EarHT", "dpoll")

myPheno |> selectTraits(all_of(myTraits))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> selectTraits(starts_with("Q"))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> selectTraits(EarHT:EarDia)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> selectTraits(EarHT, starts_with("Q"))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> selectTraits(where(is.numeric))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[, traitsWhere(traitType == "covariate")]

myPheno |> filterTraits(traitType == "covariate")


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTraits(traitType %in% c("data", "factor"))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTraits(traitType == "covariate") |> attributeData()


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTraits(traitId != "EarDia")


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
phenoDF <- as.data.frame(myPheno)

colMeans(!is.na(phenoDF[traitNames(myPheno)]))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[, traitsWhere(notMissing >= 0.95)]

myPheno |> filterTraits(notMissing >= 0.95)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[, traitsWhere(traitType == "data" & notMissing >= 0.95)]

myPheno |> filterTraits(traitType == "data", notMissing >= 0.95)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[, 1:3]

myPheno |> sliceTraits(1:3)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> sliceTraits(-1)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[, !traits("EarDia")]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[, !traitsWhere(traitType == "covariate")]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> selectTraits(-EarDia)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> selectTraits(-starts_with("Q"))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> sliceTraits(-(5:7))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTraits(traitType != "covariate")


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |> filterTaxa(location != "A")


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[
    taxaWhere(location == "A" & !is.na(EarDia)),
    traits("EarHT", "dpoll", "EarDia")
]

myPheno |>
    filterTaxa(location == "A", !is.na(EarDia)) |>
    selectTraits(EarHT, dpoll, EarDia)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
completeCases <- myPheno |>
    filterTaxa(!is.na(EarHT), !is.na(dpoll), !is.na(EarDia))

completeCases |> filterTraits(traitType == "data")


## ----eval=FALSE---------------------------------------------------------------
# myPheno[taxaWhere(location == "A"), ] |>
#     selectTraits(EarHT, dpoll)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |>
    filterTraits(notMissing >= 0.95) |>
    traitNames()

myPheno |>
    filterTaxa(!is.na(EarDia)) |>
    filterTraits(notMissing >= 0.95) |>
    traitNames()


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[, traitsWhere(notMissing >= 0.95)] |>
    traitNames()

myPheno[taxaWhere(!is.na(EarDia)), traitsWhere(notMissing >= 0.95)] |>
    traitNames()


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno[taxaWhere(location == "A"), !traits("location")]

myPheno |>
    filterTaxa(location == "A") |>
    selectTraits(-location)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myBLUE <- myPheno |>
    filterTaxa(!is.na(EarHT), !is.na(dpoll), !is.na(EarDia)) |>
    selectTraits(location, EarHT, dpoll) |>
    assocModelFitter(
        formula    = . ~ .,
        fitMarkers = FALSE
    )

myBLUE |> tableReport() |> head()


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myPheno |>
    filterTaxa(EarHT > 100) |>
    selectTraits(location, EarHT) |>
    as.data.frame()


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
genoPath <- system.file("extdata", "mdp_genotype.hmp.txt", package = "rTASSEL")

myDataset <- readGenomicDataset(genoPath, myPheno)
myDataset


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myDataset[, traits("EarHT", "Q1", "Q2", "Q3")]

myDataset |> selectTraits(EarHT, starts_with("Q"))


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# Observations: a phenotype column is named
myDataset[taxaWhere(EarHT > 100), ]

# Genotype table: only genotype metadata is named
myDataset[taxaWhere(notMissing >= 0.9), ]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myDataset |> filterTaxa(EarHT > 100)

myDataset |> filterTaxa(notMissing >= 0.9)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myDataset[taxaWhere(notMissing >= 0.9 & location == "A"), ]

myDataset |> filterTaxa(notMissing >= 0.9, location == "A")


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
tallEarsDataset <- myDataset |> filterTaxa(EarHT > 100)

tallEarsDataset |> genotype()

tallEarsDataset |> phenotype()


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myDataset[, sitesWhere(maf >= 0.05)][, traits("EarHT", "dpoll")]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myDataset |>
    filterTaxa(notMissing >= 0.9, location == "A") |>
    selectTraits(EarHT, dpoll, starts_with("Q")) |>
    filterSites(maf >= 0.05)


## ----eval=TRUE, echo=TRUE, error=TRUE-----------------------------------------
try({
myPheno |> filterTaxa(notMissing >= 0.9)
})

