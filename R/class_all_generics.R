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
#' @export
setGeneric("attributeData", function(object) standardGeneric("attributeData"))


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




