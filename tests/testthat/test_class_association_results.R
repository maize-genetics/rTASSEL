# === Tests for association class ===================================

test_that("AssociationResults methods return correct data", {
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
})

