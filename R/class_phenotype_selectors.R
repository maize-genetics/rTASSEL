# /// Selector Classes ///////////////////////////////////////////////

## ----
#' @title TraitSelector Class
#'
#' @description
#' S4 class representing trait selection criteria for bracket-based
#' filtering of \code{TasselPhenotype} and
#' \code{TasselGenomicDataset} objects.
#'
#' @slot type Character indicating selector type: \code{"names"} or
#'   \code{"predicate"}.
#' @slot ids Character vector of trait names (used when
#'   \code{type = "names"}).
#' @slot quo A quosure, or a list of quosures combined with \code{&},
#'   for predicate evaluation (used when \code{type = "predicate"}).
#' @slot negate Logical indicating whether to negate the selection.
#'
#' @name TraitSelector-class
#' @rdname TraitSelector-class
#' @exportClass TraitSelector
setClass("TraitSelector", slots = c(
    type   = "character",
    ids    = "character",
    quo    = "ANY",
    negate = "logical"
))


# /// Selector Constructors //////////////////////////////////////////

## ----
#' @title Select Traits by Name
#'
#' @description
#' Creates a \code{\linkS4class{TraitSelector}} for filtering a
#' \code{TasselPhenotype} by trait names.
#'
#' @details
#' The taxa column is the axis the observations sit on rather than a
#' trait, so it is never selectable and never dropped.
#'
#' @param ... Character trait names to select.
#'
#' @return A \code{\linkS4class{TraitSelector}} object.
#'
#' @examples
#' \dontrun{
#' ph[, traits("EarHT", "dpoll")]
#' }
#'
#' @export
traits <- function(...) {
    ids <- as.character(c(...))
    if (length(ids) == 0) rlang::abort("At least one trait name must be provided")
    methods::new("TraitSelector",
        type = "names", ids = ids, quo = NULL, negate = FALSE
    )
}


## ----
#' @title Select Traits by Predicate
#'
#' @description
#' Creates a \code{\linkS4class{TraitSelector}} using a predicate
#' expression evaluated against trait metadata. Available columns in
#' the data mask are \code{traitIndex} (1-based position of the
#' trait), \code{traitId} (trait name), \code{traitType} (TASSEL
#' attribute type), and \code{notMissing} (proportion of observations
#' that are not missing).
#'
#' @param expr An unquoted expression evaluated against trait
#'   metadata.
#'
#' @return A \code{\linkS4class{TraitSelector}} object.
#'
#' @examples
#' \dontrun{
#' ph[, traitsWhere(traitType == "covariate")]
#' ph[, traitsWhere(notMissing >= 0.95)]
#' }
#'
#' @export
traitsWhere <- function(expr) {
    predicateTraitSelector(rlang::enquo(expr))
}

## ----
# Wrap a predicate in a TraitSelector.  Shared by traitsWhere(), which
# passes a single quosure, and filterTraits(), which passes a list of them.
predicateTraitSelector <- function(quo) {
    methods::new("TraitSelector",
        type = "predicate", ids = character(0), quo = quo, negate = FALSE
    )
}


# /// Negation Methods ///////////////////////////////////////////////

## ----
#' @param x A \code{TraitSelector} object.
#'
#' @rdname TraitSelector-class
#' @aliases !,TraitSelector-method
setMethod("!", "TraitSelector", function(x) {
    x@negate <- !x@negate
    x
})
