test_that("hasDot function works correctly", {
    # Test with symbol
    expect_true(hasDot(as.symbol(".")))
    expect_false(hasDot(as.symbol("x")))
    
    # Test with calls
    expect_true(hasDot(as.call(list(as.symbol("+"), as.symbol("."), as.symbol("x")))))
    expect_false(hasDot(as.call(list(as.symbol("+"), as.symbol("x"), as.symbol("y")))))
    
    # Test with other types
    expect_false(hasDot(1))
    expect_false(hasDot("x"))
})

test_that("evalExpr handles basic symbols correctly", {
    defaultSet <- c("var1", "var2", "var3")
    
    # Test dot symbol
    result <- evalExpr(as.symbol("."), defaultSet)
    expect_equal(result$result, defaultSet)
    expect_equal(result$explicit, character(0))
    
    # Test regular symbol
    result <- evalExpr(as.symbol("x"), defaultSet)
    expect_equal(result$result, "x")
    expect_equal(result$explicit, "x")
})

test_that("evalExpr handles special predictor cases correctly", {
    df <- data.frame(
        trait_id = c("cov1", "cov2", "fct1", "fct2"),
        trait_type = c("covariate", "covariate", "factor", "factor"),
        stringsAsFactors = FALSE
    )
    
    # Test I(cov)
    result <- evalExpr(
        as.call(list(as.symbol("I"), as.symbol("cov"))),
        defaultSet = character(0),
        side = "predictor",
        df = df
    )
    expect_equal(result$result, c("cov1", "cov2"))
    expect_equal(result$explicit, "cov")
    
    # Test I(fct)
    result <- evalExpr(
        as.call(list(as.symbol("I"), as.symbol("fct"))),
        defaultSet = character(0),
        side = "predictor",
        df = df
    )
    expect_equal(result$result, c("fct1", "fct2"))
    expect_equal(result$explicit, "fct")
})

test_that("evalExpr handles operators correctly", {
    defaultSet <- c("var1", "var2", "var3")
    
    # Test addition
    result <- evalExpr(
        as.call(list(as.symbol("+"), as.symbol("x"), as.symbol("y"))),
        defaultSet
    )
    expect_equal(result$result, c("x", "y"))
    expect_equal(result$explicit, c("x", "y"))
    
    # Test subtraction
    result <- evalExpr(
        as.call(list(as.symbol("-"), as.symbol("."), as.symbol("x"))),
        defaultSet,
        isFirst = TRUE
    )
    expect_equal(result$result, c("var1", "var2", "var3")[!c("var1", "var2", "var3") %in% "x"])
    expect_equal(result$explicit, "x")
})

test_that("evalExpr handles complex nested expressions", {
    defaultSet <- c("var1", "var2", "var3")
    
    # Test nested operators
    result <- evalExpr(
        as.call(list(
            as.symbol("+"),
            as.call(list(as.symbol("-"), as.symbol("."), as.symbol("x"))),
            as.symbol("y")
        )),
        defaultSet
    )
    expect_equal(result$result, c("var1", "var2", "var3", "y"))
    expect_equal(result$explicit, c("x", "y"))
    
    # Test nested I() calls
    df <- data.frame(
        trait_id = c("cov1", "cov2", "fct1", "fct2", "data1"),
        trait_type = c("covariate", "covariate", "factor", "factor", "data"),
        stringsAsFactors = FALSE
    )
    
    result <- evalExpr(
        as.call(list(
            as.symbol("+"),
            as.call(list(as.symbol("I"), as.symbol("cov"))),
            as.call(list(as.symbol("I"), as.symbol("fct")))
        )),
        defaultSet = character(0),
        side = "predictor",
        df = df
    )
    expect_equal(result$result, c("cov1", "cov2", "fct1", "fct2"))
    expect_equal(result$explicit, c("cov", "fct"))
})

test_that("parseFormula validates formula correctly", {
    df <- data.frame(
        trait_id = c("resp1", "resp2", "cov1", "cov2", "fct1"),
        trait_type = c("data", "data", "covariate", "covariate", "factor"),
        stringsAsFactors = FALSE
    )
    
    # Test valid formula
    result <- parseFormula("resp1 ~ cov1 + fct1", df)
    expect_equal(result$response, "resp1")
    expect_equal(result$predictors, c("cov1", "fct1"))
    
    # Test formula with dot notation
    result <- parseFormula(". ~ cov1", df)
    expect_equal(result$response, c("resp1", "resp2"))
    expect_equal(result$predictors, "cov1")
    
    # Test invalid response type
    expect_error(
        parseFormula("cov1 ~ fct1", df),
        "Error in response: Only 'data' type variables are allowed"
    )
    
    # Test invalid predictor type
    expect_error(
        parseFormula("resp1 ~ resp2", df),
        "Error in predictors: 'data' type variables are not allowed"
    )
    
    # Test missing variables
    expect_warning(
        parseFormula("resp1 ~ nonexistent", df),
        "Warning in predictors: The following explicit variable\\(s\\) were not found"
    )
    
    # Test empty response
    expect_error(
        parseFormula("nonexistent ~ cov1", df),
        "Error in response: No valid response variables selected"
    )
})

test_that("parseFormula handles all special predictor types", {
    df <- data.frame(
        trait_id = c("resp1", "resp2", "cov1", "cov2", "fct1", "fct2", "data1"),
        trait_type = c("data", "data", "covariate", "covariate", 
                      "factor", "factor", "data"),
        stringsAsFactors = FALSE
    )
    
    # Test mixing special predictor types
    result <- parseFormula("resp1 ~ I(cov) + I(fct)", df)
    expect_equal(result$response, "resp1")
    expect_equal(result$predictors, c("cov1", "cov2", "fct1", "fct2"))
    
    # Test response-only formula
    result <- parseFormula("resp1 ~ 1", df)
    expect_equal(result$response, "resp1")
    expect_equal(length(result$predictors), 0)
    
    # Test dot in predictors with filtering
    result <- parseFormula("resp1 ~ . - cov1", df)
    expect_false("cov1" %in% result$predictors)
    
    # Test multiple responses with dot
    result <- parseFormula(". ~ cov1 + fct1", df)
    expect_equal(result$response, c("resp1", "resp2", "data1"))
})

test_that("hasDot handles complex nested expressions", {
    # Test dot in nested calls
    expr <- as.call(list(
        as.symbol("+"),
        as.call(list(as.symbol("-"), as.symbol("x"), as.symbol("."))),
        as.symbol("y")
    ))
    expect_true(hasDot(expr))
    
    # Test no dot in nested calls
    expr <- as.call(list(
        as.symbol("+"),
        as.call(list(as.symbol("-"), as.symbol("x"), as.symbol("z"))),
        as.symbol("y")
    ))
    expect_false(hasDot(expr))
})