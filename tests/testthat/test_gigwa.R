
exampleGigwa <- data.frame(
    `rs#` = paste0("snp_id_", 1:5),
    alleles = c("A/G", "A/G", "G/C", "G/C", "G/T"),
    chrom = c("Sb01", "Sb01", "Sb02", "Sb02", "Sb03"),
    pos = c(123, 5324, 235, 5678, 412),
    ind1 = rep(0, 5),
    ind2 = rep(1, 5),
    ind3 = rep(2, 5),
    ind4 = c(0, 1, 0, 2, 2),
    check.names = FALSE
)

test_that("`readGenotypeTableFromGigwa()` returns correct data", {
    myGt <- readGenotypeTableFromGigwa(exampleGigwa)

    expect_equal(getTaxaIDs(myGt), c("ind1", "ind2", "ind3", "ind4"))
})



