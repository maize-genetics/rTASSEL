# === Tests for association results sub-classes =====================

test_that("Method dispatch for sub classes works.", {
    ## Preamble - load data ----

    ### Start logging info
    startLogger()

    ### Load hapmap data
    genoPathHMP <- system.file(
        "extdata",
        "mdp_genotype.hmp.txt",
        package = "rTASSEL"
    )

    ### Read data - need only non missing data!
    phenoPathFast <- system.file(
        "extdata",
        "mdp_traits_nomissing.txt",
        package = "rTASSEL"
    )

    ### Create rTASSEL phenotype only object
    tasPheno <- readPhenotypeFromPath(
        path = phenoPathFast
    )

    ### Create rTASSEL genotype only object
    tasGeno <- readGenotypeTableFromPath(
        path = genoPathHMP
    )

    ### Create rTASSEL object - use prior TASSEL genotype object
    tasGenoPhenoFast <- readGenotypePhenotype(
        genoPathOrObj = genoPathHMP,
        phenoPathDFOrObj = phenoPathFast
    )

    ### Create kinship object
    tasKin <- rTASSEL::kinshipMatrix(tasGenoPhenoFast)

    ### Create error-prone kinship object
    set.seed(123)
    m <- 10
    s <- matrix(rnorm(100), m)
    s[lower.tri(s)] <- t(s)[lower.tri(s)]
    diag(s) <- 2
    colnames(s) <- rownames(s) <- paste0("s_", seq_len(m))
    tasKinError <- asTasselDistanceMatrix(s)

    ### Create multi data type pheno object
    tasGenoPhenoCov <- readGenotypePhenotype(
        genoPathOrObj = genoPathHMP,
        phenoPathDFOrObj = system.file(
            "extdata",
            "mdp_phenotype.txt",
            package = "rTASSEL"
        )
    )

    ### Create no missing data set with factors/cov (fast association)
    fullDf <- getPhenotypeDF(
        readPhenotypeFromPath(
            system.file(
                "extdata",
                "mdp_phenotype.txt",
                package = "rTASSEL"
            )
        )
    )
    noMissingDf <- fullDf[which(!is.na(fullDf$EarHT)), c("Taxa", "location", "EarHT", "Q1", "Q2", "Q3")]
    tasGenoPhenoCovNoMiss <- readGenotypePhenotype(
        genoPathOrObj = tasGeno,
        phenoPathDFOrObj = noMissingDf,
        taxaID = "Taxa",
        attributeTypes = c("factor", "data", rep("covariate", 3))
    )

    ### Association objects
    tasBLUE <- assocModelFitter(
        tasObj  = tasPheno,
        formula = . ~ .
    )
    tasGLM <- assocModelFitter(
        tasObj     = tasGenoPhenoFast,
        formula    = . ~ .,
        fitMarkers = TRUE
    )
    tasMLM <- assocModelFitter(
        tasObj     = tasGenoPhenoFast,
        formula    = . ~ .,
        fitMarkers = TRUE,
        kinship    = tasKin
    )
    tasFast <- assocModelFitter(
        tasObj          = tasGenoPhenoFast,
        formula         = . ~ .,
        fitMarkers      = TRUE,
        fastAssociation = TRUE
    )

    dimBLUE <- dim(tableReport(tasBLUE, "BLUE"))
    dimGLM  <- dim(tableReport(tasGLM, "GLM_Stats"))
    dimMLM  <- dim(tableReport(tasMLM, "MLM_Stats"))
    dimFast <- dim(tableReport(tasFast, "FastAssociation"))

    expect_true(is(tableReport(tasBLUE), "data.frame"))
    expect_true(is(tableReport(tasGLM), "data.frame"))
    expect_true(is(tableReport(tasMLM), "data.frame"))
    expect_true(is(tableReport(tasFast), "data.frame"))

    expect_equal(dim(tableReport(tasBLUE)), dimBLUE)
    expect_equal(dim(tableReport(tasGLM)), dimGLM)
    expect_equal(dim(tableReport(tasMLM)), dimMLM)
    expect_equal(dim(tableReport(tasFast)), dimFast)

})
