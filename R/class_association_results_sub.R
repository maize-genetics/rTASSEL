# === BLUE Reports ==================================================

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
#' @rdname tableReport
#' @export
setMethod(
    f = "tableReport",
    signature = signature(
        assocRes   = "AssociationResultsBLUE",
        reportName = "ANY"
    ),
    definition = function(assocRes, reportName) {
        if (missing(reportName)) reportName <- NULL
        returnReportElements(
            assocRes             = assocRes,
            reportName           = reportName,
            defaultReportElement = "BLUE"
        )
    }
)



# === GLM Reports ===================================================

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
#' @rdname tableReport
#' @export
setMethod(
    f = "tableReport",
    signature = signature(
        assocRes   = "AssociationResultsGLM",
        reportName = "ANY"
    ),
    definition = function(assocRes, reportName) {
        if (missing(reportName)) reportName <- NULL
        returnReportElements(
            assocRes             = assocRes,
            reportName           = reportName,
            defaultReportElement = "GLM_Stats"
        )
    }
)



# === MLM Reports ===================================================

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
#' @rdname tableReport
#' @export
setMethod(
    f = "tableReport",
    signature = signature(
        assocRes   = "AssociationResultsMLM",
        reportName = "ANY"
    ),
    definition = function(assocRes, reportName) {
        if (missing(reportName)) reportName <- NULL
        returnReportElements(
            assocRes             = assocRes,
            reportName           = reportName,
            defaultReportElement = "MLM_Stats"
        )
    }
)



# === Fast Association (Shabalin) Reports ===========================

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


## ----
#' @rdname tableReport
#' @export
setMethod(
    f = "tableReport",
    signature = signature(
        assocRes   = "AssociationResultsFast",
        reportName = "ANY"
    ),
    definition = function(assocRes, reportName) {
        if (missing(reportName)) reportName <- NULL
        returnReportElements(
            assocRes             = assocRes,
            reportName           = reportName,
            defaultReportElement = "FastAssociation"
        )
    }
)



# === Stepwise results ==============================================

## ----
#' @title AssociationResultsStepwise Class
#'
#' @description
#' Class \code{AssociationResultsStepwise} defines a \code{rTASSEL}
#' Class for storing TASSEL 5 Stepwise results
#'
#' @name AssociationResultsStepwise-class
#' @rdname AssociationResultsStepwise-class
#' @exportClass AssociationResultsStepwise
setClass(
    Class = "AssociationResultsStepwise",
    contains = "AssociationResults"
)


## ----
#' @rdname tableReport
#' @export
setMethod(
    f = "tableReport",
    signature = signature(
        assocRes   = "AssociationResultsStepwise",
        reportName = "ANY"
    ),
    definition = function(assocRes, reportName) {
        if (missing(reportName)) reportName <- NULL
        returnReportElements(
            assocRes             = assocRes,
            reportName           = reportName,
            defaultReportElement = "ANOVA_report"
        )
    }
)


