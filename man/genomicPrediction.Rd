% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomic_prediction.R
\name{genomicPrediction}
\alias{genomicPrediction}
\title{R interface for TASSEL's genomic prediction capabilities}
\usage{
genomicPrediction(tasPhenoObj, kinship, doCV = FALSE, kFolds, nIter)
}
\arguments{
\item{tasPhenoObj}{An object of class \code{TasselGenotypePenotype} that
contains a phenotype object.}

\item{kinship}{A TASSEL kinship object of class \code{TasselDistanceMatrix}.}

\item{doCV}{Do you want to perform k-fold cross-validation? Defaults to
\code{FALSE}.}

\item{kFolds}{Number of folds to be entered.}

\item{nIter}{Number of iterations to be ran.}
}
\value{
Returns a \code{DataFrame}-based data frame
}
\description{
This function acts as a front-end for TASSEL's genomic
  prediction functions. This analysis method uses gBLUP (genomic BLUP) to
  predict phenotypes from genotypes. It proceeds by fitting a mixed model
  that uses kinship to capture covariance between taxa. The mixed model can
  calculate BLUPs for taxa that do not have phenotypes based on the
  phenotypes of lines with relationship information.

  A phenotype dataset and a kinship matrix must be supplied as input to the
  method by selecting both then choosing Analysis/Genomic Selection. In
  addition to trait values, the phenotype dataset may also contain factors
  or covariates which will be used as fixed effects in the model. All taxa
  in the phenotype dataset can only appear once. No repeated values are
  allowed for a single taxon. When the analysis is run, the user is
  presented with the choice to run k-fold cross-validation. If cross-
  validation is selected, then the number of folds and the number of
  iterations can be entered. For each iteration and each fold within an
  iteration, the correlation between the observed and predicted values will
  be reported. If cross-validation is not selected, then the original
  observations, predicted values and PEVs (prediction error variance) will
  be reported for all taxa in the dataset.

  When k-fold cross-validation is performed, only taxa with phenotypes and
  rows in the kinship matrix are used. That set of taxa are divided into k
  subsets of equal size. Each subset in turn is used as the validation set.
  Phenotypes of the individuals in the validation are set to 0 then
  predicted using the remaining individuals as the training set. The
  correlation (r) of the observed values and predicted values is calculated
  for the validation set and reported. The mean and standard deviation of
  the mean of the r's are calculated for each trait and reported in the
  comments section of the "Accuracy" data set that is output by the
  analysis. In general, the results are not very sensitive to the choice of
  k. The number of iterations affects the standard error of the mean for the
  accuracy estimates. The defaults of k = 5 and iterations = 20 will be
  adequate for most users.
}
