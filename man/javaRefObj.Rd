% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_all_generics.R,
%   R/class_tassel_genotype.R, R/class_tassel_numeric_genotype.R,
%   R/class_tassel_phenotype.R
\name{javaRefObj}
\alias{javaRefObj}
\alias{javaRefObj,TasselGenotype-method}
\alias{javaRefObj,TasselNumericGenotype-method}
\alias{javaRefObj,TasselPhenotype-method}
\title{Return \code{rJava} reference object}
\usage{
javaRefObj(object, ...)

\S4method{javaRefObj}{TasselGenotype}(object)

\S4method{javaRefObj}{TasselNumericGenotype}(object)

\S4method{javaRefObj}{TasselPhenotype}(object)
}
\arguments{
\item{object}{an \code{rTASSEL} object}

\item{...}{Additional arguments, for use in specific methods}
}
\description{
Returns the \code{rJava} memory reference for a given \code{rTASSEL} object
}
