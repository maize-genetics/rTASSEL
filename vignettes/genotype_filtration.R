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


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
genoPath <- system.file("extdata", "mdp_genotype.hmp.txt", package = "rTASSEL")

myRealGT <- readGenotype(genoPath)
myRealGT


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT[, sitesWhere(maf >= 0.05)]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT |> filterSites(maf >= 0.05)


## ----eval=FALSE---------------------------------------------------------------
# myGT[taxa("B73", "Mo17"), ]           # myGT |> selectTaxa("B73", "Mo17")
# myGT[taxaWhere(notMissing >= 0.8), ]  # myGT |> filterTaxa(notMissing >= 0.8)
# myGT[, sites(1:100)]                  # myGT |> sliceSites(1:100)
# myGT[, !sites(1:10)]                  # myGT |> sliceSites(-(1:10))
# myGT[, siteIds(myMarkers)]            # myGT |> selectSites(all_of(myMarkers))
# myGT[, chrom("1", "5")]               # myGT |> filterSites(chrom %in% c("1", "5"))
# myGT[, region("1", 1e6, 2e6)]         # myGT |> filterSites(chrom == "1", pos >= 1e6, pos <= 2e6)
# myGT[, region(gr)]                    # myGT |> filterSites(overlaps(gr))
# myGT[, sitesWhere(maf >= 0.05 & !isIndel)]
#                                       # myGT |> filterSites(maf >= 0.05, !isIndel)


## ----eval=FALSE---------------------------------------------------------------
# myGT |> filterSites(overlaps(gr), maf >= 0.05)


## ----eval=FALSE---------------------------------------------------------------
# myGT


## ----eval=FALSE---------------------------------------------------------------
# myGT[, sitesWhere(maf >= 0.3)]


## ----eval=FALSE---------------------------------------------------------------
# myGT[, sitesWhere(maf >= 0.2 & maf <= 0.3)]


## ----eval=FALSE---------------------------------------------------------------
# myGT[, sitesWhere(alleleCount >= 10)]


## ----eval=FALSE---------------------------------------------------------------
# myGT[, sitesWhere(het <= 0.1)]


## ----eval=FALSE---------------------------------------------------------------
# myGT[, sites(2:4)]


## ----eval=FALSE---------------------------------------------------------------
# myGT[, 2:4]


## ----eval=FALSE---------------------------------------------------------------
# myGT[, chrom("1")]


## ----eval=FALSE---------------------------------------------------------------
# myGT[, region("1", 250, 542)]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT[, region("2", 20e6, 30e6)]


## ----eval=FALSE---------------------------------------------------------------
# gr


## ----eval=FALSE---------------------------------------------------------------
# myGT[, region(gr)]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
gr <- GenomicRanges::GRanges(
    seqnames = c("1", "2"),
    ranges   = IRanges::IRanges(
        start = c(20e6, 100e6),
        end   = c(30e6, 150e6)
    )
)

myRealGT[, region(gr)]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT[, sitesWhere(overlaps(gr) & maf >= 0.05)]


## $ cat my_ranges.bed
## 
## chr1    250   500
## chr2    213   400
## chr2    500   700

## ----eval=FALSE---------------------------------------------------------------
# myGT[, region(rtracklayer::import("my_ranges.bed"))]


## $ cat my_chr_pos.tsv
## 
## Chromosome  Position
## 1   300
## 2   213
## 2   665

## ----eval=FALSE---------------------------------------------------------------
# chrPos <- read.table("my_chr_pos.tsv", sep = "\t", header = TRUE)
# 
# gr <- GenomicRanges::GRanges(
#     seqnames = chrPos$Chromosome,
#     ranges   = IRanges::IRanges(start = chrPos$Position, end = chrPos$Position)
# )
# 
# myGT[, region(gr)]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT[, siteIds("PZB00859.1", "PZA01271.1")]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myMarkers <- c("PZB00859.1", "PZA01271.1")

myRealGT[, myMarkers]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT[, !sites(1:10)]


## ----eval=FALSE---------------------------------------------------------------
# myGT[taxaWhere(notMissing >= 1.0), ]


## ----eval=FALSE---------------------------------------------------------------
# myGT[taxaWhere(het <= 0.0), ]


## ----eval=FALSE---------------------------------------------------------------
# myGT[taxa("B73", "B97", "Ky21"), ]


## ----eval=FALSE---------------------------------------------------------------
# myFavTaxa <- c("B73", "B97", "Ky21")
# 
# myGT[taxa(myFavTaxa), ]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT[c("33-16", "38-11", "4226"), ]


## ----eval=FALSE---------------------------------------------------------------
# # All taxa whose ID starts with "B" or "Ky"
# myGT[taxaWhere(grepl("^B|^Ky", taxaId)), ]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT[taxaWhere(startsWith(taxaId, "CML")), ]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT[taxaWhere(startsWith(taxaId, "CML") & notMissing >= 0.9), ]


## ----eval=FALSE---------------------------------------------------------------
# myGT[taxaWhere(grepl("^B|^Ky", taxaId)), sites(2:4)]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myFiltGT <- myRealGT[taxaWhere(notMissing >= 0.8), ]
myFiltGT <- myFiltGT[, sitesWhere(maf >= 0.05)]

myFiltGT


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT |>
    filterTaxa(notMissing >= 0.8) |>
    filterSites(maf >= 0.05)


## ----eval=FALSE---------------------------------------------------------------
# myRealGT[taxaWhere(notMissing >= 0.8), sitesWhere(maf >= 0.05)] |>
#     exportGenotypeTable(
#         file   = "my_filtered_gt.vcf",
#         format = "vcf"
#     )


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
phenoPath <- system.file("extdata", "mdp_traits.txt", package = "rTASSEL")

myDataset <- readGenomicDataset(myRealGT, phenoPath)

myDataset[taxaWhere(notMissing >= 0.9), sitesWhere(maf >= 0.05)]


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myDataset |>
    filterTaxa(notMissing >= 0.9) |>
    filterSites(maf >= 0.05)


## ----eval=TRUE, echo=TRUE-----------------------------------------------------
myRealGT |> removeMinorSNPStates()


## ----eval=FALSE---------------------------------------------------------------
# myRealGT[, sitesWhere(maf >= 0.05)] |> removeMinorSNPStates()


## ----eval=FALSE---------------------------------------------------------------
# # Before
# myGT |>
#     filterGenotypeTableSites(
#         siteMinCount      = 150,
#         siteMinAlleleFreq = 0.05
#     )
# 
# # After
# myGT[, sitesWhere(alleleCount >= 150 & maf >= 0.05)]


## ----eval=FALSE---------------------------------------------------------------
# # Before
# myGT |>
#     filterGenotypeTableTaxa(minNotMissing = 0.8) |>
#     filterGenotypeTableSites(siteMinAlleleFreq = 0.05)
# 
# # After
# myGT[taxaWhere(notMissing >= 0.8), sitesWhere(maf >= 0.05)]

