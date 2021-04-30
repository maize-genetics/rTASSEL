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

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  if (!require("devtools")) install.packages("devtools")
#  devtools::install_bitbucket(
#      repo = "bucklerlab/rtassel",
#      ref = "master",
#      build_vignettes = TRUE
#  )

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  library(rTASSEL)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  options(java.parameters = c("-Xmx<memory>", "-Xms<memory>"))

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  rTASSEL::startLogger(fullPath = NULL, fileName = NULL)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
# Load hapmap data
genoPathHMP <- system.file(
    "extdata",
    "mdp_genotype.hmp.txt",
    package = "rTASSEL"
)
genoPathHMP

# Load VCF data
genoPathVCF <- system.file(
    "extdata",
    "maize_chr9_10thin40000.recode.vcf",
    package = "rTASSEL"
)
genoPathVCF

