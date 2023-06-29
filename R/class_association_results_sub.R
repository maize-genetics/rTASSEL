## ----
#' @title AssociationResultsBLUE Class
#'
#' @description
#' Class \code{AssociationResultsBLUE} defines a \code{rTASSEL}
#' Class for storing TASSEL 5 BLUE GWAS results
#'
#' @name AssociationResultsBLUE-class
#' @rdname AssociationResultsBLUE-class
#' @exportClass AssociationResultsBLUE
setClass(
    Class = "AssociationResultsBLUE",
    contains = "AssociationResults"
)


## ----
#' @title AssociationResultsGLM Class
#'
#' @description
#' Class \code{AssociationResultsGLM} defines a \code{rTASSEL}
#' Class for storing TASSEL 5 GLM GWAS results
#'
#' @name AssociationResultsGLM-class
#' @rdname AssociationResultsGLM-class
#' @exportClass AssociationResultsGLM
setClass(
    Class = "AssociationResultsGLM",
    contains = "AssociationResults"
)


## ----
#' @title AssociationResultsMLM Class
#'
#' @description
#' Class \code{AssociationResultsMLM} defines a \code{rTASSEL}
#' Class for storing TASSEL 5 MLM GWAS results
#'
#' @name AssociationResultsMLM-class
#' @rdname AssociationResultsMLM-class
#' @exportClass AssociationResultsMLM
setClass(
    Class = "AssociationResultsMLM",
    contains = "AssociationResults"
)


## ----
#' @title AssociationResultsFast Class
#'
#' @description
#' Class \code{AssociationResultsFast} defines a \code{rTASSEL}
#' Class for storing TASSEL 5 FastAssociation GWAS results
#'
#' @name AssociationResultsFast-class
#' @rdname AssociationResultsFast-class
#' @exportClass AssociationResultsFast
setClass(
    Class = "AssociationResultsFast",
    contains = "AssociationResults"
)


