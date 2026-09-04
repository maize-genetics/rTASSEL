# === Tests for genotype table coercion =============================

test_that("Genotype table coercion to matrix returns correct data", {
    m <- as.matrix(rtObjs$gt_hmp)

    expect_equal(dim(m), c(281, 3093))
})

test_that("Coercion is unaffected by the class holding the genotype table", {
    expect_equal(
        as.matrix(rtObjs$gt_hmp),
        as.matrix(rtObjsLegacy$gt_hmp)
    )
})
