## ----
#' @name TasselGenotype-class
#' @rdname TasselGenotype-class
#' @exportClass TasselGenotype
setClass(
    Class = "TasselGenotype",
    slots = c(
        dispData    = "matrix",
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

        rJc <- rJava::.jnew(TASSEL_JVM$R_METHODS)

        return(rJc$read(xNorm, keepDepth, sortPositions))
    }
}


