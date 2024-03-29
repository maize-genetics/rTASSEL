% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotype_table_functions.R
\name{getSumExpFromGenotypeTable}
\alias{getSumExpFromGenotypeTable}
\title{Create Summarized Experiment from a TASSEL Genotype Table}
\usage{
getSumExpFromGenotypeTable(tasObj, coerceDosageToInt = TRUE, verbose = FALSE)
}
\arguments{
\item{tasObj}{An object of class \code{TasselGenotypePenotype}.}

\item{coerceDosageToInt}{Should dosage array be returned as \code{integer}
values? If \code{FALSE}, dosage array will be returned as type
\code{raw} byte values. Returning \code{raw} byte values. Will greatly
save on memory. Defaults to \code{TRUE}.}

\item{verbose}{Should messages be displayed to console? Defaults to
\code{FALSE}.}
}
\value{
Returns a \code{SummarizedExperiment} of TASSEL genotype data.
}
\description{
This function will generate an object of
   \code{SummarizedExperiment} class for marker data derived from a
   \code{TasselGenotypePhenotype} class object.
}
