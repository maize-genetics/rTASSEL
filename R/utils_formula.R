## ----
# Recursively check whether an expression contains a dot (".")
#
# @description
# This function recursively checks whether a given R expression 
# contains the dot symbol (`.`). It can handle symbols, calls, and 
# other types of expressions.
#
# @param expr
# An R expression to be checked. This can be a symbol, a call, 
# or other types of expressions.
#
# @return
# A logical value: `TRUE` if the dot symbol (`.`) is present in the 
# expression, and `FALSE` otherwise.
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
# Evaluate an Expression in the Context of a Formula
#
# @description
# This function evaluates an expression within the context of a 
# formula, supporting operations such as symbol resolution, special 
# handling for predictors, and binary/unary operators. It is designed 
# to process expressions in formula parsing, particularly for 
# response and predictor sides.
#
# @details
# - Symbols are resolved to their character representation. A dot 
#   (`.`) is replaced with the `defaultSet` on the first occurrence, 
#   and treated as identity (no effect) on subsequent occurrences.
# - Special handling is provided for `I(cov)` and `I(fct)` on the 
#   predictor side, which resolve to covariate and factor sets, 
#   respectively.
# - Binary operators (`+`, `-`) are supported for combining or 
#   subtracting sets of variables. Unary minus (`-`) is interpreted 
#   as subtracting from the `defaultSet` when it appears at the 
#   beginning.
# - Calls to `c(...)` are evaluated recursively, combining results 
#   from each argument.
# - Unsupported expression types or invalid usage will result in an 
#   error.
# 
# @param expr
# An R expression to evaluate. This can be a symbol, a call, or 
# other supported expression types.
# @param defaultSet
# A character vector representing the default set of variables to 
# use when the expression contains a dot (`.`).
# @param isFirst
# Logical, indicating whether this is the first occurrence of the 
# expression in the formula. Defaults to `TRUE`.
# @param side
# A string indicating the side of the formula being evaluated. Can 
# be `"response"` or `"predictor"`.
# @param df
# A data frame containing metadata about traits. Required when 
# evaluating special cases like `I(cov)` or `I(fct)` on the predictor 
# side.
#
# @return
# A list with two elements:
#   - result.....: A character vector of resolved variable names or 
#                  sets.
#   - explicit...: A character vector of explicitly mentioned 
#                  variable names in the expression.
#
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
                rlang::abort("I() call in predictor side must have exactly one argument")
            }
            arg <- expr[[2]]
            if (!is.symbol(arg)) {
                rlang::abort("I() call in predictor side must have a symbol argument (cov or fct)")
            }
            token <- as.character(arg)
            if (token == "cov") {
                specialSet <- df$trait_id[df$trait_type == "covariate"]
            } else if (token == "fct") {
                specialSet <- df$trait_id[df$trait_type == "factor"]
            } else {
                rlang::abort("I() call in predictor side must be either I(cov) or I(fct)")
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
# Parse a Formula String and Validate Variables
#
# @description
# This function parses a formula string and validates the variables 
# used in the formula against a provided data frame of valid 
# variables and their types. It supports the use of a dot (`.`) in 
# the formula to represent default sets of response and predictor 
# variables.
#
# @details
# The function performs the following steps:
#   1. Parses the formula string into a formula object.
#   2. Identifies the left-hand side (response) and right-hand side 
#      (predictors) of the formula.
#   3. Evaluates the formula, injecting default variables if a dot 
#      (`.`) is present.
#   4. Validates the explicit variables in the formula against the 
#      valid variables in `df`.
#   5. Ensures type consistency:
#      - Response variables must be of type `"data"`.
#      - Predictor variables cannot be of type `"data"`.
#   6. Returns the validated response and predictor variables.
#
# @param formulaStr
# A character string representing the formula to be parsed. For 
# example, `"response ~ predictors"`.
# @param df
# A data frame containing valid variable information. It must have 
# the following columns:
#   - `trait_id`: Character vector of variable names.
#   - `trait_type`: Character vector indicating the type of each 
#     variable. Expected values are:
#       - `"data"`: Variables allowed in the response.
#       - `"covariate"` or `"factor"`: Variables allowed in the 
#         predictors.
#
# @note
# - If no valid response variables are selected, an error is raised.
# - If no valid predictor variables are selected, the function does 
#   not raise an error but returns an empty vector for predictors.
# - Warnings are issued for any explicit variables in the formula 
#   that are not found in `df`.
# 
# @return
# A list with two elements:
#   - `response`: A character vector of valid response variables.
#   - `predictors`: A character vector of valid predictor variables.
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
                rlang::abort(sprintf("Error in %s: All explicit variable(s) not found in attributes: %s",
                             sideName, paste(missingVars, collapse = ", ")))
            } else {
                rlang::warn(sprintf("Warning in %s: The following explicit variable(s) were not found in attributes and will be ignored: %s",
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
        rlang::abort("Error in response: No valid response variables selected.")
    }
    # Note: The check for predictorFinal having zero length has been removed.

    # Type consistency checks.
    # For response, only variables with trait_type "data" are allowed.
    wrongTypeResponse <- df$trait_type[df$trait_id %in% responseFinal & df$trait_type != "data"]
    if (length(wrongTypeResponse) > 0) {
        rlang::abort("Error in response: Only 'data' type variables are allowed in the response. Found non-'data' types.")
    }

    # For predictors, variables of type "data" are not permitted.
    wrongTypePredictor <- df$trait_type[df$trait_id %in% predictorFinal & df$trait_type == "data"]
    if (length(wrongTypePredictor) > 0) {
        rlang::abort("Error in predictors: 'data' type variables are not allowed in predictors. Only 'covariate' and 'factor' types are permitted.")
    }

    return(list(response = responseFinal, predictors = predictorFinal))
}


