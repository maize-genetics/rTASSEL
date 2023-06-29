# === Tests for association class ===================================

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
})

