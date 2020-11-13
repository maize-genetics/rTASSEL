## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
    fig.path='figure/graphics-',
    cache.path='cache/graphics-',
    fig.align='center',
    external=TRUE,
    echo=TRUE,
    warning=FALSE
    # fig.pos="H"
)

## ---- eval=FALSE--------------------------------------------------------------
#  myGT

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableSites(
#          siteMinAlleleFreq = 0.3,
#          siteMaxAlleleFreq = 1.0
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableSites(
#          siteMinAlleleFreq = 0.2,
#          siteMaxAlleleFreq = 0.3
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableSites(
#          siteMinCount = 5
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableSites(
#          minHeterozygous = 0.0,
#          maxHeterozygous = 0.1
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableSites(
#          siteRangeFilterType = "sites",
#          startSite = 1,
#          endSite = 3
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableSites(
#          siteRangeFilterType = "position",
#          startChr = 1,
#          endChr = 1,
#          startPos = NULL,
#          endPos = NULL
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableSites(
#          siteRangeFilterType = "position",
#          startChr = 1,
#          endChr = 2,
#          startPos = 250,
#          endPos = 700
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  gr

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableSites(
#          siteRangeFilterType = "none",
#          gRangesObj = gr
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableSites(
#          siteRangeFilterType = "none",
#          bedFile = "my_ranges.bed"
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableSites(
#          siteRangeFilterType = "none",
#          bedFile = "my_chr_pos.tsv"
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableTaxa(
#          minNotMissing = 1.0
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableTaxa(
#          minHeterozygous = 0.0,
#          maxHeterozygous = 0.0
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myFavTaxa

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableTaxa(
#          taxa = myFavTaxa
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  allTaxa

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableTaxa(
#          taxa = allTaxa %>% str_subset("^B|^Ky")
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableTaxa(
#          taxa = allTaxa %>% str_subset("^B|^Ky")
#      ) %>%
#      filterGenotypeTableSites(
#          siteRangeFilterType = "sites",
#          startSite = 1,
#          endSite = 3
#      )

## ---- eval=FALSE--------------------------------------------------------------
#  myGT %>%
#      filterGenotypeTableTaxa(
#          taxa = allTaxa %>% str_subset("^B|^Ky")
#      ) %>%
#      filterGenotypeTableSites(
#          siteRangeFilterType = "sites",
#          startSite = 1,
#          endSite = 3
#      ) %>%
#      exportGenotypeTable(
#          file = "my_filtered_gt.vcf"
#          format = "vcf"
#      )

