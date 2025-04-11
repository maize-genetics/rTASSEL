# Helper: recursively check whether an expression contains a dot (".")
hasDot <- function(expr) {
    if (is.symbol(expr)) {
        return(as.character(expr) == ".")
    } else if (is.call(expr)) {
        return(any(sapply(as.list(expr[-1]), hasDot)))
    } else {
        return(FALSE)
    }
}

# Recursive evaluator for left-to-right “chain” evaluation.
#
# This version returns a list with two elements:
#   - result: a character vector of variable names produced by evaluating the expression.
#   - explicit: the set of variable names that were explicitly mentioned (i.e. not injected via a default dot).
#
# Parameters:
#   - expr: the expression to evaluate.
#   - defaultSet: the set to inject if a dot is encountered as the first token.
#   - isFirst: a logical flag indicating whether the current term is first in a chain.
#   - side: a string ("response" or "predictor"). For side == "predictor", special tokens I(cov) and I(fct) are supported.
#   - df: if side is "predictor", the attribute data frame must be provided.
evalExpr <- function(expr, defaultSet, isFirst = TRUE, side = "response", df = NULL) {
    if (is.symbol(expr)) {
        varName <- as.character(expr)
        if (varName == ".") {
            if (isFirst) {
                return(list(result = defaultSet, explicit = character(0)))
            } else {
                # Subsequent occurrences of dot are treated as identity (no effect).
                return(list(result = character(0), explicit = character(0)))
            }
        } else {
            return(list(result = varName, explicit = varName))
        }
    } else if (is.call(expr)) {
        op <- as.character(expr[[1]])

        # --- Special case for predictor side: I(cov) and I(fct) ---
        if (side == "predictor" && op == "I") {
            if (length(expr) != 2) {
                stop("I() call in predictor side must have exactly one argument")
            }
            arg <- expr[[2]]
            if (!is.symbol(arg)) {
                stop("I() call in predictor side must have a symbol argument (cov or fct)")
            }
            token <- as.character(arg)
            if (token == "cov") {
                specialSet <- df$trait_id[df$trait_type == "covariate"]
            } else if (token == "fct") {
                specialSet <- df$trait_id[df$trait_type == "factor"]
            } else {
                stop("I() call in predictor side must be either I(cov) or I(fct)")
            }
            # The special token returns its corresponding set. Its explicit token is recorded as the token value.
            return(list(result = specialSet, explicit = token))
        }

        # --- Process other calls ---
        if (op == "c") {
            # Evaluate each argument of c(...), each treated as first.
            evalList <- lapply(as.list(expr[-1]),
                               function(e) evalExpr(e, defaultSet, isFirst = TRUE, side = side, df = df))
            resultVec <- unique(unlist(lapply(evalList, function(x) x$result)))
            explicitVec <- unique(unlist(lapply(evalList, function(x) x$explicit)))
            return(list(result = resultVec, explicit = explicitVec))
        } else if (op %in% c("+", "-") && length(expr) == 3) {
            # For binary operators, evaluate the left operand (preserving isFirst flag) and the right operand with isFirst = FALSE.
            leftEval <- evalExpr(expr[[2]], defaultSet, isFirst, side = side, df = df)
            rightEval <- evalExpr(expr[[3]], defaultSet, isFirst = FALSE, side = side, df = df)
            if (op == "+") {
                resultVec <- union(leftEval$result, rightEval$result)
            } else {  # op == "-"
                resultVec <- setdiff(leftEval$result, rightEval$result)
            }
            explicitVec <- unique(c(leftEval$explicit, rightEval$explicit))
            return(list(result = resultVec, explicit = explicitVec))
        } else if (op == "-" && length(expr) == 2) {
            # A unary minus at the beginning is interpreted as: defaultSet - operand.
            if (isFirst) {
                operandEval <- evalExpr(expr[[2]], defaultSet, isFirst = TRUE, side = side, df = df)
                resultVec <- setdiff(defaultSet, operandEval$result)
                return(list(result = resultVec, explicit = operandEval$explicit))
            } else {
                stop("Unary minus encountered in non-first position is not supported.")
            }
        } else {
            stop("Unknown expression type encountered during formula parsing.")
        }
    } else {
        stop("Unsupported expression type encountered.")
    }
}

# Main function: parseFormula
#
# Parameters:
#   - formulaStr: a character string specifying the formula.
#       e.g. "c(dpoll, d3) ~ . - I(cov)" or ". ~ Q1 + I(fct)"
#   - df: an attribute data frame with at least the following columns:
#           * trait_id: variable identifiers,
#           * trait_type: classification (e.g., "data", "covariate", "factor", etc.)
#
# Behavior:
#   - Response (LHS):
#       * If a dot is present then the default is all variables with trait_type "data".
#       * Otherwise, only the explicit tokens (e.g. via c(...) or plain symbols) are used.
#
#   - Predictors (RHS):
#       * If a dot is present then the default is all variables with trait_type in c("covariate", "factor").
#       * Otherwise, only the explicit tokens are used.
#       * In addition, the special tokens I(cov) and I(fct) (allowed only on the predictor side)
#         return the sets of variables with trait_type "covariate" and "factor" respectively.
#
# After evaluating each side, the function:
#   1. Checks the explicit tokens against the attribute data:
#         - If some explicit variables are missing (i.e. not found in df$trait_id) and at least one valid explicit exists
#           or a default is present, a warning is issued and those missing names are dropped.
#         - If all explicit tokens on a side are missing, an error is raised.
#   2. Enforces type consistency:
#         - The response must contain only variables with trait_type "data".
#         - The predictors must not include any variables with trait_type "data".
parseFormula <- function(formulaStr, df) {
    # Valid variable names from the attribute data.
    validVars <- df$trait_id

    # Default sets for response and predictors.
    defaultResponse <- df$trait_id[df$trait_type == "data"]
    defaultPredictors <- df$trait_id[df$trait_type %in% c("covariate", "factor")]

    # Convert the formula string to a formula object.
    f <- as.formula(formulaStr)
    lhsExpr <- f[[2]]  # left-hand side (response)
    rhsExpr <- f[[3]]  # right-hand side (predictors)

    # Check if a dot exists on each side.
    responseHasDefault <- hasDot(lhsExpr)
    predictorHasDefault <- hasDot(rhsExpr)

    # For evaluation, if a dot is present then inject the corresponding default; otherwise use an empty default.
    responseDefaultForEval <- if (responseHasDefault) defaultResponse else character(0)
    predictorDefaultForEval <- if (predictorHasDefault) defaultPredictors else character(0)

    # Evaluate each side.
    lhsEval <- evalExpr(lhsExpr, responseDefaultForEval, isFirst = TRUE, side = "response")
    rhsEval <- evalExpr(rhsExpr, predictorDefaultForEval, isFirst = TRUE, side = "predictor", df = df)

    responseVars <- unique(lhsEval$result)
    predictorVars <- unique(rhsEval$result)

    # Helper function to check explicit tokens on a given side.
    # 'sideName' is "response" or "predictors", 'sideHasDefault' indicates if a dot was present.
    checkMissingExplicit <- function(explicit, sideName, sideHasDefault) {
        # For predictors, filter out the special tokens "cov" and "fct".
        if (sideName == "predictors") {
            explicit <- explicit[!(explicit %in% c("cov", "fct"))]
        }
        missingVars <- setdiff(explicit, validVars)
        if (length(explicit) > 0 && length(missingVars) > 0) {
            # If no valid explicit variable remains and no default is injected, raise an error.
            if (!sideHasDefault && length(explicit) == length(missingVars)) {
                stop(sprintf("Error in %s: All explicit variable(s) not found in attributes: %s",
                             sideName, paste(missingVars, collapse = ", ")))
            } else {
                warning(sprintf("Warning in %s: The following explicit variable(s) were not found in attributes and will be ignored: %s",
                                sideName, paste(missingVars, collapse = ", ")))
            }
        }
    }

    # Check explicit tokens on each side.
    checkMissingExplicit(lhsEval$explicit, "response", responseHasDefault)
    checkMissingExplicit(rhsEval$explicit, "predictors", predictorHasDefault)

    # Remove any names not found in the attribute data.
    responseFinal <- intersect(responseVars, validVars)
    predictorFinal <- intersect(predictorVars, validVars)

    if (length(responseFinal) == 0) {
        stop("Error in response: No valid response variables selected.")
    }
    # Note: The check for predictorFinal having zero length has been removed.

    # Type consistency checks.
    # For response, only variables with trait_type "data" are allowed.
    wrongTypeResponse <- df$trait_type[df$trait_id %in% responseFinal & df$trait_type != "data"]
    if (length(wrongTypeResponse) > 0) {
        stop("Error in response: Only 'data' type variables are allowed in the response. Found non-'data' types.")
    }

    # For predictors, variables of type "data" are not permitted.
    wrongTypePredictor <- df$trait_type[df$trait_id %in% predictorFinal & df$trait_type == "data"]
    if (length(wrongTypePredictor) > 0) {
        stop("Error in predictors: 'data' type variables are not allowed in predictors. Only 'covariate' and 'factor' types are permitted.")
    }

    return(list(response = responseFinal, predictors = predictorFinal))
}
# Helper: recursively check whether an expression contains a dot (".")
hasDot <- function(expr) {
    if (is.symbol(expr)) {
        return(as.character(expr) == ".")
    } else if (is.call(expr)) {
        return(any(sapply(as.list(expr[-1]), hasDot)))
    } else {
        return(FALSE)
    }
}

# Recursive evaluator for left-to-right “chain” evaluation.
#
# This version returns a list with two elements:
#   - result: a character vector of variable names produced by evaluating the expression.
#   - explicit: the set of variable names that were explicitly mentioned (i.e. not injected via a default dot).
#
# Parameters:
#   - expr: the expression to evaluate.
#   - defaultSet: the set to inject if a dot is encountered as the first token.
#   - isFirst: a logical flag indicating whether the current term is first in a chain.
#   - side: a string ("response" or "predictor"). For side == "predictor", special tokens I(cov) and I(fct) are supported.
#   - df: if side is "predictor", the attribute data frame must be provided.
evalExpr <- function(expr, defaultSet, isFirst = TRUE, side = "response", df = NULL) {
    if (is.symbol(expr)) {
        varName <- as.character(expr)
        if (varName == ".") {
            if (isFirst) {
                return(list(result = defaultSet, explicit = character(0)))
            } else {
                # Subsequent occurrences of dot are treated as identity (no effect).
                return(list(result = character(0), explicit = character(0)))
            }
        } else {
            return(list(result = varName, explicit = varName))
        }
    } else if (is.call(expr)) {
        op <- as.character(expr[[1]])

        # --- Special case for predictor side: I(cov) and I(fct) ---
        if (side == "predictor" && op == "I") {
            if (length(expr) != 2) {
                stop("I() call in predictor side must have exactly one argument")
            }
            arg <- expr[[2]]
            if (!is.symbol(arg)) {
                stop("I() call in predictor side must have a symbol argument (cov or fct)")
            }
            token <- as.character(arg)
            if (token == "cov") {
                specialSet <- df$trait_id[df$trait_type == "covariate"]
            } else if (token == "fct") {
                specialSet <- df$trait_id[df$trait_type == "factor"]
            } else {
                stop("I() call in predictor side must be either I(cov) or I(fct)")
            }
            # The special token returns its corresponding set. Its explicit token is recorded as the token value.
            return(list(result = specialSet, explicit = token))
        }

        # --- Process other calls ---
        if (op == "c") {
            # Evaluate each argument of c(...), each treated as first.
            evalList <- lapply(as.list(expr[-1]),
                               function(e) evalExpr(e, defaultSet, isFirst = TRUE, side = side, df = df))
            resultVec <- unique(unlist(lapply(evalList, function(x) x$result)))
            explicitVec <- unique(unlist(lapply(evalList, function(x) x$explicit)))
            return(list(result = resultVec, explicit = explicitVec))
        } else if (op %in% c("+", "-") && length(expr) == 3) {
            # For binary operators, evaluate the left operand (preserving isFirst flag) and the right operand with isFirst = FALSE.
            leftEval <- evalExpr(expr[[2]], defaultSet, isFirst, side = side, df = df)
            rightEval <- evalExpr(expr[[3]], defaultSet, isFirst = FALSE, side = side, df = df)
            if (op == "+") {
                resultVec <- union(leftEval$result, rightEval$result)
            } else {  # op == "-"
                resultVec <- setdiff(leftEval$result, rightEval$result)
            }
            explicitVec <- unique(c(leftEval$explicit, rightEval$explicit))
            return(list(result = resultVec, explicit = explicitVec))
        } else if (op == "-" && length(expr) == 2) {
            # A unary minus at the beginning is interpreted as: defaultSet - operand.
            if (isFirst) {
                operandEval <- evalExpr(expr[[2]], defaultSet, isFirst = TRUE, side = side, df = df)
                resultVec <- setdiff(defaultSet, operandEval$result)
                return(list(result = resultVec, explicit = operandEval$explicit))
            } else {
                stop("Unary minus encountered in non-first position is not supported.")
            }
        } else {
            stop("Unknown expression type encountered during formula parsing.")
        }
    } else {
        stop("Unsupported expression type encountered.")
    }
}

# Main function: parseFormula
#
# Parameters:
#   - formulaStr: a character string specifying the formula.
#       e.g. "c(dpoll, d3) ~ . - I(cov)" or ". ~ Q1 + I(fct)"
#   - df: an attribute data frame with at least the following columns:
#           * trait_id: variable identifiers,
#           * trait_type: classification (e.g., "data", "covariate", "factor", etc.)
#
# Behavior:
#   - Response (LHS):
#       * If a dot is present then the default is all variables with trait_type "data".
#       * Otherwise, only the explicit tokens (e.g. via c(...) or plain symbols) are used.
#
#   - Predictors (RHS):
#       * If a dot is present then the default is all variables with trait_type in c("covariate", "factor").
#       * Otherwise, only the explicit tokens are used.
#       * In addition, the special tokens I(cov) and I(fct) (allowed only on the predictor side)
#         return the sets of variables with trait_type "covariate" and "factor" respectively.
#
# After evaluating each side, the function:
#   1. Checks the explicit tokens against the attribute data:
#         - If some explicit variables are missing (i.e. not found in df$trait_id) and at least one valid explicit exists
#           or a default is present, a warning is issued and those missing names are dropped.
#         - If all explicit tokens on a side are missing, an error is raised.
#   2. Enforces type consistency:
#         - The response must contain only variables with trait_type "data".
#         - The predictors must not include any variables with trait_type "data".
parseFormula <- function(formulaStr, df) {
    # Valid variable names from the attribute data.
    validVars <- df$trait_id

    # Default sets for response and predictors.
    defaultResponse <- df$trait_id[df$trait_type == "data"]
    defaultPredictors <- df$trait_id[df$trait_type %in% c("covariate", "factor")]

    # Convert the formula string to a formula object.
    f <- as.formula(formulaStr)
    lhsExpr <- f[[2]]  # left-hand side (response)
    rhsExpr <- f[[3]]  # right-hand side (predictors)

    # Check if a dot exists on each side.
    responseHasDefault <- hasDot(lhsExpr)
    predictorHasDefault <- hasDot(rhsExpr)

    # For evaluation, if a dot is present then inject the corresponding default; otherwise use an empty default.
    responseDefaultForEval <- if (responseHasDefault) defaultResponse else character(0)
    predictorDefaultForEval <- if (predictorHasDefault) defaultPredictors else character(0)

    # Evaluate each side.
    lhsEval <- evalExpr(lhsExpr, responseDefaultForEval, isFirst = TRUE, side = "response")
    rhsEval <- evalExpr(rhsExpr, predictorDefaultForEval, isFirst = TRUE, side = "predictor", df = df)

    responseVars <- unique(lhsEval$result)
    predictorVars <- unique(rhsEval$result)

    # Helper function to check explicit tokens on a given side.
    # 'sideName' is "response" or "predictors", 'sideHasDefault' indicates if a dot was present.
    checkMissingExplicit <- function(explicit, sideName, sideHasDefault) {
        # For predictors, filter out the special tokens "cov" and "fct".
        if (sideName == "predictors") {
            explicit <- explicit[!(explicit %in% c("cov", "fct"))]
        }
        missingVars <- setdiff(explicit, validVars)
        if (length(explicit) > 0 && length(missingVars) > 0) {
            # If no valid explicit variable remains and no default is injected, raise an error.
            if (!sideHasDefault && length(explicit) == length(missingVars)) {
                stop(sprintf("Error in %s: All explicit variable(s) not found in attributes: %s",
                             sideName, paste(missingVars, collapse = ", ")))
            } else {
                warning(sprintf("Warning in %s: The following explicit variable(s) were not found in attributes and will be ignored: %s",
                                sideName, paste(missingVars, collapse = ", ")))
            }
        }
    }

    # Check explicit tokens on each side.
    checkMissingExplicit(lhsEval$explicit, "response", responseHasDefault)
    checkMissingExplicit(rhsEval$explicit, "predictors", predictorHasDefault)

    # Remove any names not found in the attribute data.
    responseFinal <- intersect(responseVars, validVars)
    predictorFinal <- intersect(predictorVars, validVars)

    if (length(responseFinal) == 0) {
        stop("Error in response: No valid response variables selected.")
    }
    # Note: The check for predictorFinal having zero length has been removed.

    # Type consistency checks.
    # For response, only variables with trait_type "data" are allowed.
    wrongTypeResponse <- df$trait_type[df$trait_id %in% responseFinal & df$trait_type != "data"]
    if (length(wrongTypeResponse) > 0) {
        stop("Error in response: Only 'data' type variables are allowed in the response. Found non-'data' types.")
    }

    # For predictors, variables of type "data" are not permitted.
    wrongTypePredictor <- df$trait_type[df$trait_id %in% predictorFinal & df$trait_type == "data"]
    if (length(wrongTypePredictor) > 0) {
        stop("Error in predictors: 'data' type variables are not allowed in predictors. Only 'covariate' and 'factor' types are permitted.")
    }

    return(list(response = responseFinal, predictors = predictorFinal))
}
