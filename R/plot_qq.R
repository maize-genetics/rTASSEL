## ----
#' @title Create a QQ plot from rTASSEL association output
#'
#' @description This function allows for quick generation of a QQ
#'    plot from rTASSEL association statistical output data.
#'
#' @param assocRes An object of type \code{AssociationResults}
#' @param trait Which phenotypic trait do you want to plot? If set to
#'    \code{NULL}, this will generate a faceted plot with all mapped traits.
#' @param overlay Do you want trait results faceted or overlayed into one
#'    single plot? Defaults to \code{TRUE}.
#'
#' @return Returns a \code{ggplot2} object
#'
#' @export
plotQQ <- function(
    assocRes,
    trait = NULL,
    overlay = TRUE
) {
    if (!is(assocRes, "AssociationResults")) {
        stop(
            "The object '", deparse(substitute(assocRes)),
            "' is not an 'AssociationResults' object"
        )
    }

    traitValidityChecker(trait, assocRes)

    assocStats <- switch (
        associationType(assocRes),
        "GLM"       = tableReport(assocRes, "GLM_Stats"),
        "MLM"       = tableReport(assocRes, "MLM_Stats"),
        "FastAssoc" = tableReport(assocRes, "FastAssociation"),
        "default"   = NULL
    )

    if (is.null(assocStats)) {
        stop("Association Type not defined")
    }

    # Make a more robust/'OOP' object to pass
    pQQCoreParams <- list(
        "assocStats" = assocStats,
        "trait"      = trait,
        "overlay"    = overlay
    )

    plotQQCore(pQQCoreParams)
}


## ----
#' @title Core visual engine for QQ plotting
#' @param params A list of parameter variables
#' @importFrom rlang .data
plotQQCore <- function(params) {

    ## Parse parameters
    assocStats <- params$assocStats
    trait      <- params$trait
    overlay    <- params$overlay

    ## "Prime" and filter data for plotting
    qqStats <- generateQQData(params)

    ## Check for overlay
    if (overlay) {
        aesData <- list(
            "point_aes" = ggplot2::geom_point(
                size = 0.8,
                na.rm = TRUE
            ),
            "primaryAes" = ggplot2::aes(
                x = -log10(.data$expected),
                y = -log10(.data$observed),
                color = .data$trait_id
            ),
            "trait_labs" = ggplot2::labs(
                color = "Trait"
            )
        )
    } else {
        aesData <- list(
            "facets" = ggplot2::facet_wrap(
                facets = . ~ trait_id,
                scales = "free",
                ncol = 3
            ),
            "point_aes" = ggplot2::geom_point(
                size = 0.8,
                na.rm = TRUE,
                color = "#3e619b"
            ),
            "primaryAes" = ggplot2::aes(
                x = -log10(.data$expected),
                y = -log10(.data$observed),
            ),
            "theme" = ggplot2::theme(legend.position = "none")
        )
    }

    ## Plot components
    p <- ggplot2::ggplot(data = qqStats) +
        aesData$primaryAes +
        ggplot2::geom_abline(
            intercept = 0,
            slope = 1
        ) +
        aesData$point_aes +
        ggplot2::xlab(bquote(~-log[10]~ '('*italic(p)*'-value) (expected)')) +
        ggplot2::ylab(bquote(~-log[10]~ '('*italic(p)*'-value) (observed)')) +
        aesData$facets +
        aesData$trait_labs +
        ggplot2::theme_bw() +
        aesData$theme

    return(p)
}


## ----
#' @title Generate QQ plot data
#' @param params A list of parameter variables
generateQQData <- function(params) {
    ## Parse parameters
    assocStats <- params$assocStats
    trait      <- params$trait

    ## Sanity check coerce data frame from assocStats object
    # NOTE - probably don't need to do this right now but a security
    #        feature if I change TableReport classes in the future...
    assocStats <- as.data.frame(assocStats)

    ## Sanity check for columns
    neededCols <- c("Trait", "p")
    checkForValidColumns(assocStats, neededCols)

    ## Coerce association stats trait to character to stop ghost levels
    assocStats       <- assocStats[, c("Trait", "p")]
    assocStats$Trait <- as.character(assocStats$Trait)
    assocStats$p     <- as.numeric(assocStats$p)

    ## Check trait levels
    if (!is.null(trait)) {
        assocStats <- assocStats[assocStats$Trait %in% trait, ]
    }


    ## Generate QQ results
    qqDf <- do.call(
        what = "rbind",
        args = lapply(
            X = split(assocStats, assocStats$Trait),
            FUN = function(i) {
                expectedQuantiles(i$p, unique(i$Trait))
            }
        )
    )
    rownames(qqDf) <- NULL

    return(qqDf)
}


## ----
#' @title Get expected quantiles from p-val vectors
#' @param pVal A \code{numeric} vector containing association p-values
#' @param trait A trait ID that correlates with p-value data
expectedQuantiles <- function(pVal, trait) {
    ## Remove NaN values
    filtPVal <- pVal[!is.na(pVal)]

    ## CDF(p) = P[measurement with p-value of ≤ p] = p
    n        <- length(filtPVal)
    obsQuant <- sort(filtPVal)
    expQuant <- seq_along(filtPVal) / n # {1/n, 2/n, ... n/n}

    ## Return
    eDf <- data.frame(
        trait_id = trait,
        expected = expQuant,
        observed = obsQuant
    )
    return(eDf)
}


