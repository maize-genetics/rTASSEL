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


