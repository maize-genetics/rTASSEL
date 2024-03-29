% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_qq.R
\name{plotQQ}
\alias{plotQQ}
\title{Create a QQ plot from rTASSEL association output}
\usage{
plotQQ(assocRes, trait = NULL, overlay = TRUE, interactive = FALSE)
}
\arguments{
\item{assocRes}{An object of type \code{AssociationResults}}

\item{trait}{Which phenotypic trait do you want to plot? If set to
\code{NULL}, this will generate a faceted plot with all mapped traits.}

\item{overlay}{Do you want trait results faceted or overlayed into one
single plot? Defaults to \code{TRUE}.}

\item{interactive}{Do you want to produce an interactive visualization?
Defaults to \code{FALSE}.}
}
\value{
Returns a \code{ggplot2} object
}
\description{
This function allows for quick generation of a QQ
   plot from rTASSEL association statistical output data.
}
