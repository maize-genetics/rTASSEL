## ----
#' @title Get taxa IDs from genotype data
#'
#' @description
#' Returns a character vector of taxa (sample) IDs from a
#' \code{\linkS4class{TasselGenotypePhenotype}} or
#' \code{\linkS4class{TasselGenotype}} object.
#'
#' @param tasObj A \code{TasselGenotypePhenotype} or \code{TasselGenotype}
#'   object containing genotype data.
#'
#' @return A character vector of taxa IDs.
#'
#' @rdname taxaList
#' @export
setGeneric("taxaList", function(tasObj) standardGeneric("taxaList"))


## ----
#' @title Get position list metadata from genotype data
#'
#' @description
#' Returns positional metadata (chromosome, position, etc.) from a
#' \code{\linkS4class{TasselGenotypePhenotype}} or
#' \code{\linkS4class{TasselGenotype}} object.
#'
#' @param tasObj A \code{TasselGenotypePhenotype} or \code{TasselGenotype}
#'   object containing genotype data.
#'
#' @return A \code{tibble} of positional metadata.
#'
#' @rdname positionList
#' @export
setGeneric("positionList", function(tasObj) standardGeneric("positionList"))


## ----
#' @title Get site summary of genotype table
#'
#' @description
#' Returns per-site summary statistics (allele frequencies, heterozygosity,
#' missingness, etc.) from genotype data stored in a
#' \code{\linkS4class{TasselGenotypePhenotype}} or
#' \code{\linkS4class{TasselGenotype}} object.
#'
#' @param tasObj A \code{TasselGenotypePhenotype} or \code{TasselGenotype}
#'   object containing genotype data.
#'
#' @return A \code{data.frame} of per-site summary statistics.
#'
#' @rdname siteSummary
#' @export
setGeneric("siteSummary", function(tasObj) standardGeneric("siteSummary"))


## ----
#' @title Get taxa summary of genotype table
#'
#' @description
#' Returns per-taxon summary statistics (missingness, heterozygosity, etc.)
#' from genotype data stored in a
#' \code{\linkS4class{TasselGenotypePhenotype}} or
#' \code{\linkS4class{TasselGenotype}} object.
#'
#' @param tasObj A \code{TasselGenotypePhenotype} or \code{TasselGenotype}
#'   object containing genotype data.
#'
#' @return A \code{data.frame} of per-taxon summary statistics.
#'
#' @rdname taxaSummary
#' @export
setGeneric("taxaSummary", function(tasObj) standardGeneric("taxaSummary"))


## ----
#' @title Return GWAS association type
#'
#' @description
#' Returns association type for given  \code{\linkS4class{AssociationResults}}
#' object.
#'
#' @param assocRes a \code{\linkS4class{AssociationResults}} object
#'
#' @rdname associationType
#' @export
setGeneric("associationType", function(assocRes) standardGeneric("associationType"))


## ----
#' @title Return TASSEL attribute data
#'
#' @description
#' Returns attribute data from an rTASSEL phenotype object
#'
#' @param object an \code{rTASSEL} object
#' @param ... Additional arguments, for use in specific methods
#'
#' @rdname attributeData
#' @export
setGeneric("attributeData", function(object, ...) standardGeneric("attributeData"))


## ----
#' @title Return \code{rJava} reference object
#'
#' @description
#' Returns the \code{rJava} memory reference for a given \code{rTASSEL} object
#'
#' @param object an \code{rTASSEL} object
#' @param ... Additional arguments, for use in specific methods
#'
#' @rdname javaRefObj
#' @export
setGeneric("javaRefObj", function(object, ...) standardGeneric("javaRefObj"))


## ----
#' @title Return report names
#'
#' @description
#' Returns a \code{character} vector of table report names
#'
#' @param object a \code{\linkS4class{AssociationResults}} object
#'
#' @rdname reportNames
#' @export
setGeneric("reportNames", function(object) standardGeneric("reportNames"))


## ----
#' @title Return selected table report
#'
#' @description
#' Returns a \code{data.frame} object of association table reports
#'
#' @param assocRes a \code{\linkS4class{AssociationResults}} object
#' @param reportName a specific table report to return
#'
#' @rdname tableReport
#' @export
setGeneric("tableReport", function(assocRes = missing(), reportName = missing()) standardGeneric("tableReport"))


## ----
#' @title Return trait names
#'
#' @description
#' Returns a \code{character} vector of trait names
#'
#' @param object a \code{\linkS4class{AssociationResults}} object
#'
#' @rdname traitNames
#' @export
setGeneric("traitNames", function(object) standardGeneric("traitNames"))




