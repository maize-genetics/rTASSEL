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
#  devtools::install_github(
#      repo = "maize-genetics/rTASSEL",
#      ref = "master",
#      build_vignettes = TRUE,
#      dependencies = TRUE
#  )

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  if (!require("devtools")) install.packages("devtools")
#  devtools::install_github("maize-genetics/rTASSEL")

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
library(rTASSEL)

## ---- eval=FALSE, echo = TRUE-------------------------------------------------
#  devtools::install_github(
#      repo = "maize-genetics/rTASSEL",
#      ref = "master",
#      build_vignettes = FALSE,
#      INSTALL_opts = "--no-multiarch"
#  )

