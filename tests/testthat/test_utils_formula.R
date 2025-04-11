context("Unit tests for utils_formula.R functions")

## Tests for the hasDot() function ----

test_that("hasDot correctly identifies a dot symbol", {
    expect_true(hasDot(quote(.)))           # Single dot symbol should return TRUE
})

test_that("hasDot returns FALSE for non-dot symbols", {
    expect_false(hasDot(quote(a)))
})

test_that("hasDot detects dot recursively in call expressions", {
    # In c(., a), the dot is part of the expression
    expect_true(hasDot(quote(c(., a))))

    # In c(a, b) there is no dot in any element
    expect_false(hasDot(quote(c(a, b))))
})


## Tests for the evalExpr() function ----

test_that("evalExpr handles symbols and dot correctly", {
    default_set <- c("x", "y")

    # For a regular symbol that is not a dot:
    res_symbol <- evalExpr(quote(a), default_set)
    expect_equal(res_symbol$result, "a")
    expect_equal(res_symbol$explicit, "a")

    # Dot encountered in the first position should return the default set:
    res_dot_first <- evalExpr(quote(.), default_set, isFirst = TRUE)
    expect_equal(res_dot_first$result, default_set)
    expect_equal(res_dot_first$explicit, character(0))

    # Dot encountered when not in the first position should return an empty result:
    res_dot_nonfirst <- evalExpr(quote(.), default_set, isFirst = FALSE)
    expect_equal(res_dot_nonfirst$result, character(0))
    expect_equal(res_dot_nonfirst$explicit, character(0))
})

test_that("evalExpr handles the 'c' call", {
    default_set <- c("x", "y")

    # c(a, .) should combine the explicit symbol "a" and the default set for "."
    res_c_call <- evalExpr(quote(c(a, .)), default_set)
    expected_result <- unique(c("a", default_set))
    expect_equal(sort(res_c_call$result), sort(expected_result))
    expect_equal(res_c_call$explicit, "a")
})

test_that("evalExpr handles binary operators '+' and '-'", {
    default_set <- c("x", "y")

    # Binary plus: a + b should return both symbols.
    res_plus <- evalExpr(quote(a + b), default_set)
    expect_equal(sort(res_plus$result), sort(c("a", "b")))
    expect_equal(sort(res_plus$explicit), sort(c("a", "b")))

    # Binary minus: a - b should return the set difference (i.e. "a" if they are distinct)
    res_minus <- evalExpr(quote(a - b), default_set)
    expect_equal(res_minus$result, "a")
    expect_equal(sort(res_minus$explicit), sort(c("a", "b")))
})

test_that("evalExpr handles unary minus correctly", {
    default_set <- c("a", "b", "c")

    # Unary minus in first position: -a should remove "a" from the default set.
    res_unary <- evalExpr(quote(-a), default_set, isFirst = TRUE)
    expect_equal(sort(res_unary$result), sort(setdiff(default_set, "a")))
    expect_equal(res_unary$explicit, "a")

    # Unary minus in a non-first position should trigger an error.
    expect_error(evalExpr(quote(-a), default_set, isFirst = FALSE),
                 "Unary minus encountered in non-first position is not supported.")
})

test_that("evalExpr aborts on unknown operator", {
    default_set <- c("x", "y")

    # An expression with an unsupported operator (like *) should abort.
    expect_error(evalExpr(quote(a * b), default_set),
                 "Unknown expression type encountered during formula parsing.")
})


## Tests for the parseFormula() function ----

df <- data.frame(
    trait_id = c("a", "b", "c", "d", "e", "e.2"),
    trait_type = c("data", "covariate", "factor", "data", "covariate", "data"),
    stringsAsFactors = FALSE
)

test_that("parseFormula works with explicit variable names", {
    # Both response and predictors are provided explicitly.
    res <- parseFormula("a ~ b + c", df)
    expect_equal(res$response, "a")
    expect_equal(sort(res$predictors), sort(c("b", "c")))
})

test_that("parseFormula uses default response when dot is present on LHS", {
    # Dot on the left-hand side should default to all "data" variables.
    res <- parseFormula(". ~ b", df)
    # Default response should be trait_id with trait_type "data": a and d.
    expect_equal(sort(res$response), sort(c("a", "d")))
    expect_equal(res$predictors, "b")
})

test_that("parseFormula uses default predictors when dot is present on RHS", {
    # Dot on the right-hand side should default to all covariates and factors.
    res <- parseFormula("a ~ .", df)
    expect_equal(res$response, "a")
    # Default predictors are trait_id where trait_type is either covariate or factor: b, c, and e.
    expect_equal(sort(res$predictors), sort(c("b", "c", "e")))
})

test_that("parseFormula issues a warning when some explicit predictors are missing", {
    # Construct a smaller data frame.
    df_small <- data.frame(
        trait_id = c("a", "b", "c"),
        trait_type = c("data", "covariate", "factor"),
        stringsAsFactors = FALSE
    )
    # Predictor "nonexistent" is not in df_small; since not all predictors are missing, a warning should be issued.
    expect_warning(
        res <- parseFormula("a ~ b + nonexistent", df_small),
        "Warning in predictors:"
    )
    # The missing predictor should be omitted.
    expect_equal(res$response, "a")
    expect_equal(res$predictors, "b")
})

test_that("parseFormula aborts when explicit response variable is missing", {
    # "nonexistent" is not found in the data frame, so expect an error.
    expect_error(
        parseFormula("nonexistent ~ b", df),
        "Error in response: All explicit variable\\(s\\) not found in attributes:"
    )
})

test_that("parseFormula enforces type checking for response variables", {
    # Create a data frame where the response variable is not of type 'data'
    df_invalid_response <- data.frame(
        trait_id = c("a", "b"),
        trait_type = c("covariate", "factor"),
        stringsAsFactors = FALSE
    )
    expect_error(
        parseFormula("a ~ b", df_invalid_response),
        "Error in response: Only 'data' type variables are allowed in the response."
    )
})

test_that("parseFormula enforces type checking for predictor variables", {
    # Create a data frame where the predictor variable is of type 'data', which is not allowed.
    df_invalid_predictor <- data.frame(
        trait_id = c("a", "b"),
        trait_type = c("data", "data"),
        stringsAsFactors = FALSE
    )
    expect_error(
        parseFormula("a ~ b", df_invalid_predictor),
        "Error in predictors: 'data' type variables are not allowed in predictors."
    )
})


