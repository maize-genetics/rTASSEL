## ----
#' @title Remove minor SNP states from genotype data
#'
#' @description
#' Collapses every site down to its two most common allelic states,
#' recoding any remaining (minor) states as missing. This is the
#' bracket-API replacement for the \code{removeMinorSNPStates} argument of
#' the deprecated \code{\link{filterGenotypeTableSites}()}: it transforms
#' genotype calls rather than selecting a subset of taxa or sites, so it
#' does not fit the \code{\link{sitesWhere}()} selector model.
#'
#' @param x An object of class \code{\linkS4class{TasselGenotype}} or
#'    \code{\linkS4class{TasselGenomicDataset}}. Objects of the deprecated
#'    \code{TasselGenotypePhenotype} class are still accepted.
#'
#' @return An object of the same class as \code{x} holding the recoded
#'    genotype data.
#'
#' @examples
#' \dontrun{
#' gt <- readGenotype("path/to/genotype.hmp.txt")
#' gt |> removeMinorSNPStates()
#'
#' # Chains with bracket filtering
#' gt[, sitesWhere(maf >= 0.05)] |> removeMinorSNPStates()
#' }
#'
#' @importFrom rJava new
#' @importFrom rJava J
#' @importFrom rJava .jnull
#'
#' @export
removeMinorSNPStates <- function(x) {
    tasIn <- .resolveTasselInput(x, "genotype", "removeMinorSNPStates")

    plugin <- rJava::new(
        rJava::J("net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin"),
        rJava::.jnull(),
        FALSE
    )
    plugin$setParameter("removeMinorSNPStates", "true")

    jGt <- tryCatch(
        plugin$runPlugin(tasIn$jGt),
        error = function(e) NULL
    )

    if (!inherits(jGt, "jobjRef") || rJava::is.jnull(jGt)) {
        rlang::abort("No sites remain after removing minor SNP states")
    }

    .wrapGenotypeResult(jGt, tasIn)
}
