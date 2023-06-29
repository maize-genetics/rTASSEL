# === Tests for AssociationClass objects ============================

test_that("AssociationResults methods return correct data and exceptions", {
    # Test instantiation
    testClass <- methods::new(
        "AssociationResults",
        results = list(),
        traits = "TEST"
    )
    expect_true(is(testClass, "AssociationResults"))
    expect_error(
        object = methods::new(
            "AssociationResults",
            results = list(mtcars, letters),
            traits = "TEST"
        ),
        regexp = "Invalid list: Element"
    )
    expect_error(
        object = methods::new(
            "AssociationResults",
            results = list(mtcars, letters, 1:5),
            traits = "TEST"
        ),
        regexp = "Invalid list: Elements"
    )
    expect_error(
        object = methods::new(
            "AssociationResults",
            results = list(mtcars, iris),
            traits = 12
        ),
        regexp = "should be or extend class \"character\""
    )

    # Test show method
    classShowFew <- capture_output(
        show(
            methods::new(
                "AssociationResults",
                results = list("td_1" = iris, "td_2" = mtcars),
                traits = paste0("trait_", 1:3)
            )
        )
    )
    classShowMany <- capture_output(
        show(
            methods::new(
                "AssociationResults",
                results = list("td_1" = iris, "td_2" = mtcars),
                traits = paste0("trait_", 1:10)
            )
        )
    )
    expect_true(grepl("AssociationResults object with 2", classShowFew))
    expect_true(grepl("(with 5 more values)", classShowMany))

    # Test getters
    testClass2 <- methods::new(
        "AssociationResults",
        results = list("td_1" = iris, "td_2" = mtcars),
        traits = paste0("trait_", 1:8)
    )

    expect_equal(reportNames(testClass2), c("td_1", "td_2"))
    expect_equal(traitNames(testClass2), paste0("trait_", 1:8))

    expect_equal(tableReport(testClass2, "td_1"), iris)
    expect_equal(tableReport(testClass2, "td_2"), mtcars)
    expect_true(is(tableReport(testClass2), "list"))
    expect_equal(length(tableReport(testClass2)), 2)
    expect_error(
        tableReport(testClass2, 123),
        regexp = "'reportName' must be of type 'character'"
    )
    expect_error(
        tableReport(testClass2, "td_3"),
        regexp = "Report ID not found in object"
    )
})


test_that("AssociationResults sub classes return correct data.", {
    testClassBLUE <- methods::new(
        "AssociationResultsBLUE",
        results = list("td_1" = iris, "td_2" = mtcars),
        traits = paste0("trait_", 1:8),
        assocType = "BLUE"
    )
    testClassGLM <- methods::new(
        "AssociationResultsGLM",
        results = list("td_1" = iris, "td_2" = mtcars),
        traits = paste0("trait_", 1:8),
        assocType = "GLM"
    )
    testClassMLM <- methods::new(
        "AssociationResultsMLM",
        results = list("td_1" = iris, "td_2" = mtcars),
        traits = paste0("trait_", 1:8),
        assocType = "MLM"
    )
    testClassFast <- methods::new(
        "AssociationResultsFast",
        results = list("td_1" = iris, "td_2" = mtcars),
        traits = paste0("trait_", 1:8),
        assocType = "FastAssociation"
    )

    expect_true(is(testClassBLUE, "AssociationResultsBLUE"))
    expect_true(is(testClassGLM, "AssociationResultsGLM"))
    expect_true(is(testClassMLM, "AssociationResultsMLM"))
    expect_true(is(testClassFast, "AssociationResultsFast"))

    expect_equal(associationType(testClassBLUE), "BLUE")
    expect_equal(associationType(testClassGLM), "GLM")
    expect_equal(associationType(testClassMLM), "MLM")
    expect_equal(associationType(testClassFast), "FastAssociation")
})




























