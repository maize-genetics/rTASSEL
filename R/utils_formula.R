## ----
# Check recursively whether an expression contains a dot (".")
hasDot <- function(expr) {
    if (is.symbol(expr)) {
        return(as.character(expr) == ".")
    } else if (is.call(expr)) {
        return(any(sapply(as.list(expr[-1]), hasDot)))
    } else {
        return(FALSE)
    }
}


## ----
# Recursively evaluate an expression in a left-to-right “chain” manner.
# Returns a list with two elements:
#   - result: the character vector of variable names from evaluating the expression.
#   - explicit: the set of names that were explicitly mentioned (ignoring injected default via ".").
# The parameter 'defaultSet' is used if a dot is encountered as the very first token.
# 'isFirst' indicates whether the token is the first element in a chain.
evalExpr <- function(expr, defaultSet, isFirst = TRUE) {
    if (is.symbol(expr)) {
        varName <- as.character(expr)
        if (varName == ".") {
            if (isFirst) {
                return(list(result = defaultSet, explicit = character(0)))
            } else {
                return(list(result = character(0), explicit = character(0)))
            }
        } else {
            return(list(result = varName, explicit = varName))
        }
    } else if (is.call(expr)) {
        op <- as.character(expr[[1]])
        if (op == "c") {
            evalList    <- lapply(as.list(expr[-1]), function(e) evalExpr(e, defaultSet, isFirst = TRUE))
            resultVec   <- unique(unlist(lapply(evalList, function(x) x$result)))
            explicitVec <- unique(unlist(lapply(evalList, function(x) x$explicit)))
            return(list(result = resultVec, explicit = explicitVec))
        } else if (op %in% c("+", "-") && length(expr) == 3) {
            leftEval  <- evalExpr(expr[[2]], defaultSet, isFirst)
            rightEval <- evalExpr(expr[[3]], defaultSet, isFirst = FALSE)
            if (op == "+") {
                resultVec <- union(leftEval$result, rightEval$result)
            } else {
                resultVec <- setdiff(leftEval$result, rightEval$result)
            }
            explicitVec <- unique(c(leftEval$explicit, rightEval$explicit))
            return(list(result = resultVec, explicit = explicitVec))
        } else if (op == "-" && length(expr) == 2) {
            if (isFirst) {
                operandEval <- evalExpr(expr[[2]], defaultSet, isFirst = TRUE)
                resultVec   <- setdiff(defaultSet, operandEval$result)
                return(list(result = resultVec, explicit = operandEval$explicit))
            } else {
                rlang::abort("Unary minus encountered in non-first position is not supported.")
            }
        } else {
            rlang::abort("Unknown expression type encountered during formula parsing.")
        }
    } else {
        rlang::abort("Unsupported expression type encountered.")
    }
}


## ----
# Main function to parse a formula into response and predictor variable sets.
#
# Parameters:
#   - formulaStr: a character string specifying the formula (e.g. "c(dpoll, d3) ~ . + d5")
#   - df: a data frame of attributes that must contain at least:
#           * trait_id: variable identifiers
#           * trait_type: classification ("data", "covariate", "factor", etc.)
#
# Behavior:
#  - For the response (left-hand side), if a dot is present, then the default set is
#    all trait_id with trait_type "data"; otherwise, only the explicit tokens are used.
#  - For the predictors (right-hand side), if a dot is present, then the default set is
#    all trait_id with trait_type in c("covariate", "factor"); otherwise, only the explicit
#    tokens are used.
#
# After evaluation, explicit tokens are checked:
#   * If some explicit variable names are not in df$trait_id, they are removed.
#   * If all explicit names on a side are missing, an error is raised.
#   * Otherwise a warning is issued listing the missing variables.
#
# Finally, type checks are enforced:
#   - The response may contain only "data" type variables.
#   - The predictors may not contain any "data" type variables.
parseFormula <- function(formulaStr, df) {
    validVars         <- df$trait_id
    defaultResponse   <- df$trait_id[df$trait_type == "data"]
    defaultPredictors <- df$trait_id[df$trait_type %in% c("covariate", "factor")]

    f       <- as.formula(formulaStr)
    lhsExpr <- f[[2]]
    rhsExpr <- f[[3]]

    responseHasDefault  <- hasDot(lhsExpr)
    predictorHasDefault <- hasDot(rhsExpr)

    responseDefaultForEval  <- if (responseHasDefault) defaultResponse else character(0)
    predictorDefaultForEval <- if (predictorHasDefault) defaultPredictors else character(0)

    lhsEval <- evalExpr(lhsExpr, responseDefaultForEval, isFirst = TRUE)
    rhsEval <- evalExpr(rhsExpr, predictorDefaultForEval, isFirst = TRUE)

    responseVars  <- unique(lhsEval$result)
    predictorVars <- unique(rhsEval$result)

    checkMissingExplicit <- function(explicit, sideName, sideHasDefault) {
        missingVars <- setdiff(explicit, validVars)
        if (length(explicit) > 0 && length(missingVars) > 0) {
            if (!sideHasDefault && length(explicit) == length(missingVars)) {
                rlang::abort(
                    message = sprintf(
                        "Error in %s: All explicit variable(s) not found in attributes: %s",
                        sideName, paste(missingVars, collapse = ", ")
                    )
                )
            } else {
                rlang::warn(
                    message = sprintf(
                        "Warning in %s: The following explicit variable(s) were not found in attributes and will be ignored: %s",
                        sideName, paste(missingVars, collapse = ", ")
                    )
                )
            }
        }
    }

    checkMissingExplicit(lhsEval$explicit, "response", responseHasDefault)
    checkMissingExplicit(rhsEval$explicit, "predictors", predictorHasDefault)

    responseFinal  <- intersect(responseVars, validVars)
    predictorFinal <- intersect(predictorVars, validVars)

    if (length(responseFinal) == 0) {
        rlang::abort("Error in response: No valid response variables selected.")
    }
    if (length(predictorFinal) == 0) {
        rlang::abort("Error in predictors: No valid predictor variables selected.")
    }

    wrongTypeResponse <- df$trait_type[df$trait_id %in% responseFinal & df$trait_type != "data"]
    if (length(wrongTypeResponse) > 0) {
        rlang::abort("Error in response: Only 'data' type variables are allowed in the response. Found non-'data' types.")
    }

    wrongTypePredictor <- df$trait_type[df$trait_id %in% predictorFinal & df$trait_type == "data"]
    if (length(wrongTypePredictor) > 0) {
        rlang::abort("Error in predictors: 'data' type variables are not allowed in predictors. Only 'covariate' and 'factor' types are permitted.")
    }

    return(list(response = responseFinal, predictors = predictorFinal))
}


