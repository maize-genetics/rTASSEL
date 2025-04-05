## ----
#' @name TasselGenotype-class
#' @rdname TasselGenotype-class
#' @exportClass TasselGenotype
setClass(
    Class = "TasselGenotype",
    slots = c(
        dispData    = "list",
        jGeno       = "jobjRef",
        jMemAddress = "character",
        jClass      = "character"
    )
)


## ----
# targets:
#   * path
#   * matrix
#   * data frame
#' @title Helper function to build a TasselGenotype object
#' @export
readGenotype <- function(x, sortPositions = FALSE, keepDepth = FALSE) {
    if (is.character(x)) {
        xNorm <- normalizePath(x)
        if (!file.exists(xNorm)) {
            rlang::abort("The input path is not a valid file")
        }

        rJc         <- rJava::.jnew(TASSEL_JVM$R_METHODS)
        javaGt      <- rJc$read(xNorm, keepDepth, sortPositions)
        jClass      <- rJava::.jclass(javaGt)
        jMemAddress <- gsub(".*@", "", rJava::.jstrVal(javaGt))

        methods::new(
            Class = "TasselGenotype",
            dispData    = formatGtStrings(javaGt),
            jGeno       = javaGt,
            jMemAddress = jMemAddress,
            jClass      = jClass
        )
    }
}



# /// Methods (show) ////////////////////////////////////////////////

#' @export
setMethod("show", "TasselGenotype", function(object) {
    printGtDisp(
        fgs    = object@dispData,
        nTaxa  = object@jGeno$numberOfTaxa(),
        nSites = object@jGeno$numberOfSites(),
        jMem   = object@jMemAddress
    )
})


