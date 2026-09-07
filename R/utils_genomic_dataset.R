# /// S3 - helper functions /////////////////////////////////////////

## ----
# Create a vector class for genotype display data
#
# @description
# Wraps per-observation genotype calls in a `vctrs` vector so that a
# genotype column can travel alongside the phenotype columns of a
# `tibble`. Values are stored raw; the styling is applied by
# `pillar_shaft.geno()` at print time so that `cli` can detect the
# terminal's color capabilities.
#
# @param x
# A list with one element per displayed observation. Each element is
# either a `character` vector of allele calls or a `numeric` vector
# of reference probabilities.
# @param nSites
# The total number of sites in the genotype table.
# @param nSitesShown
# The number of sites held in each element of `x`.
# @param minorAlleles
# A `character` vector of minor alleles, one per shown site. Only
# used for allele-based genotypes.
# @param numeric
# A `logical` marking `x` as reference probabilities rather than
# allele calls.
#
# @return
# A `geno` vector of the same length as `x`.
genoVctr <- function(
    x,
    nSites,
    nSitesShown,
    minorAlleles = character(),
    numeric = FALSE
) {
    vctrs::new_vctr(
        x,
        nSites       = nSites,
        nSitesShown  = nSitesShown,
        minorAlleles = minorAlleles,
        numeric      = numeric,
        class        = "geno"
    )
}


## ----
# Build the genotype column of a genomic dataset display table
#
# @description
# Collects the leading sites of a joined genotype table for the
# observations that will be displayed. Reads through the Java
# `GenotypePhenotype` object rather than the genotype table so that
# replicated and phenotype-only taxa line up with the phenotype rows.
#
# @param jGp
# A Java `GenotypePhenotype` object reference.
# @param jGt
# A Java `GenotypeTable` object reference, used for minor alleles.
# @param nRows
# The number of leading observations to collect.
# @param nSites
# The maximum number of sites to collect per observation.
#
# @return
# A `geno` vector of length `nRows`.
buildGenotypeColumn <- function(jGp, jGt, nRows, nSites = 5) {
    isDiscrete  <- jGp$areGenotypeValuesDiscrete()
    nSitesTotal <- jGt$numberOfSites()
    nShow       <- min(nSites, nSitesTotal)

    if (nRows == 0 || nShow == 0) {
        return(
            genoVctr(
                x           = vector("list", nRows),
                nSites      = nSitesTotal,
                nSitesShown = 0L,
                numeric     = !isDiscrete
            )
        )
    }

    siteIdx <- seq_len(nShow) - 1L
    rowIdx  <- seq_len(nRows)

    siteVals <- lapply(siteIdx, function(site) {
        if (isDiscrete) {
            jGp$getStringGenotype(as.integer(site))[rowIdx]
        } else {
            jGp$referenceProb(as.integer(site))[rowIdx]
        }
    })

    # A matrix keeps the single-observation case from collapsing to a
    # vector, which would scramble the per-row split below
    siteMat <- matrix(unlist(siteVals), nrow = nRows, ncol = nShow)

    minorAlleles <- if (isDiscrete) {
        vapply(siteIdx, function(site) {
            jGt$minorAlleleAsString(as.integer(site))
        }, FUN.VALUE = character(1))
    } else {
        character()
    }

    genoVctr(
        x            = lapply(rowIdx, function(i) siteMat[i, ]),
        nSites       = nSitesTotal,
        nSitesShown  = nShow,
        minorAlleles = minorAlleles,
        numeric      = !isDiscrete
    )
}


## ----
# Create a Java-Compatible Genomic Dataset Table
#
# @description
# Creates a tibble from the phenotype columns and genotype column of
# a genomic dataset and assigns the attributes needed by the header
# and footer methods. The resulting object is assigned a custom class
# `"java_geno_pheno_tbl"`.
#
# @details
# The `Genotype` column is registered as pillar's "focus" column so
# that it keeps its place at the end of the table when the console is
# too narrow to hold every phenotype column. Phenotype columns are
# squeezed into the footer instead.
#
# @param data
# A data frame or object that can be converted to a tibble.
# @param nTaxa
# An integer specifying the number of taxa.
# @param nSites
# An integer specifying the number of sites.
# @param nSitesShown
# An integer specifying the number of sites held in the genotype
# column.
# @param nCap
# An integer specifying the maximum number of rows to display.
# @param nDfRow
# An integer specifying the number of observations in the dataset.
# @param colTypes
# A single `character` value summarizing the column types.
# @param jMem
# A character value specifying the Java memory address.
#
# @return
# A tibble with additional attributes (`nTaxa`, `nSites`,
# `nSitesShown`, `nCap`, `nDfRow`, `colTypes`, `jMem`) and a custom
# class `"java_geno_pheno_tbl"`.
javaGenoPhenoTbl <- function(
    data,
    nTaxa,
    nSites,
    nSitesShown,
    nCap,
    nDfRow,
    colTypes,
    jMem
) {
    df <- tibble::as_tibble(data)
    attr(df, "nTaxa")       <- nTaxa
    attr(df, "nSites")      <- nSites
    attr(df, "nSitesShown") <- nSitesShown
    attr(df, "nCap")        <- nCap
    attr(df, "nDfRow")      <- nDfRow
    attr(df, "colTypes")    <- colTypes
    attr(df, "jMem")        <- jMem

    if ("Genotype" %in% colnames(df)) {
        attr(df, "pillar_focus") <- "Genotype"
    }

    class(df) <- c("java_geno_pheno_tbl", class(df))

    df
}


## ----
# Collapse a column type summary into a single display string
#
# @description
# Counts the columns of a genomic dataset display table by type and
# renders them in a fixed order so that the footer reads the same way
# for every dataset.
#
# @param attrSummary
# The `attrSummary` slot of a `TasselPhenotype`: a named list of
# trait counts keyed by TASSEL attribute type.
# @param hasGenotype
# A `logical` marking whether a genotype column is displayed.
#
# @return
# A single `character` value, e.g.
# `"taxa: 1, data: 3, genotype: 1"`.
formatColumnTypeSummary <- function(attrSummary, hasGenotype = TRUE) {
    counts <- unlist(attrSummary)

    if (hasGenotype) {
        counts <- c(counts, "genotype" = 1L)
    }

    if (length(counts) == 0) {
        return("no columns")
    }

    # A fixed order keeps the footer stable; anything TASSEL adds later
    # is appended rather than dropped
    canonical <- c("taxa", "factor", "data", "covariate", "genotype")
    ordered   <- c(
        intersect(canonical, names(counts)),
        setdiff(names(counts), canonical)
    )

    paste0(ordered, ": ", counts[ordered], collapse = ", ")
}


## ----
# Format Genomic Dataset Display
#
# @description
# Builds the table printed by `show()` for a `TasselGenomicDataset`
# object by pairing the phenotype's already-formatted display columns
# with a genotype column.
#
# @details
# The phenotype's `dispData` slot is reused directly, so the trait
# columns carry the same vector classes - and therefore the same
# formatting - as they do when a `TasselPhenotype` is printed on its
# own. The table is built on demand rather than cached in a slot so
# that subsetting and joining a dataset stay cheap.
#
# @param object
# A `TasselGenomicDataset` object.
# @param nSites
# The number of leading sites to show in the genotype column.
#
# @return
# A `java_geno_pheno_tbl` object.
formatGenomicDatasetDisplay <- function(object, nSites = 5) {
    phDisp <- object@phenotype@dispData
    jGt    <- object@genotype@jRefObj

    gtCol <- buildGenotypeColumn(
        jGp    = object@jRefObj,
        jGt    = jGt,
        nRows  = nrow(phDisp),
        nSites = nSites
    )

    tblData <- c(as.list(phDisp), list("Genotype" = gtCol))

    javaGenoPhenoTbl(
        data        = tibble::as_tibble(tblData),
        nTaxa       = jGt$numberOfTaxa(),
        nSites      = jGt$numberOfSites(),
        nSitesShown = attr(gtCol, "nSitesShown"),
        nCap        = attr(phDisp, "nCap"),
        nDfRow      = attr(phDisp, "nDfRow"),
        colTypes    = formatColumnTypeSummary(object@phenotype@attrSummary),
        jMem        = object@jMemAddress
    )
}
