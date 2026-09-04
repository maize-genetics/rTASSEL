# /// Constants /////////////////////////////////////////////////////

## ----
# Input classes accepted by rTASSEL's data-bearing functions
#
# @description
# 'MODERN' lists the wrapper classes that are the primary path as of
# rTASSEL 0.14.0. 'LEGACY' is kept working behind a deprecation
# warning and is scheduled for removal in the next major release.
TASSEL_INPUT <- list(
    "MODERN" = c("TasselGenotype", "TasselPhenotype", "TasselGenomicDataset"),
    "LEGACY" = "TasselGenotypePhenotype"
)



# /// Class predicates //////////////////////////////////////////////

## ----
# Is an object an instance of any of the given S4 classes?
#
# @description
# Unlike 'inherits()', 'methods::is()' always honors S4 inheritance,
# so 'TasselNumericGenotype' objects are matched by 'TasselGenotype'.
# Class names that are not currently defined simply do not match.
#
# @param x
# Any R object.
# @param classes
# A 'character' vector of S4 class names.
#
# @return
# A single 'logical' value.
.isAnyClass <- function(x, classes) {
    any(vapply(classes, function(cl) methods::is(x, cl), logical(1)))
}


## ----
# Unwrap an rTASSEL wrapper class to its backing Java reference
#
# @description
# The modern wrapper classes all keep their Java object in '@jRefObj'
# ('GenotypeTable', 'Phenotype', and 'GenotypePhenotype',
# respectively), so unwrapping funnels them into the '%instanceof%'
# logic that the house-keeping getters already apply to raw 'jobjRef'
# input. Any other object is returned untouched.
#
# @param x
# Any R object.
#
# @return
# A 'jobjRef' if 'x' is a modern wrapper class, otherwise 'x'.
.unwrapTasselObject <- function(x) {
    if (.isAnyClass(x, TASSEL_INPUT$MODERN)) {
        return(x@jRefObj)
    }

    return(x)
}



# /// Shared input resolution ///////////////////////////////////////

## ----
# Resolve any accepted rTASSEL object into its Java components
#
# @description
# Single entry point replacing the copy-pasted
# 'inherits(tasObj, "TasselGenotypePhenotype")' guards and direct
# '@jGenotypeTable' / '@jPhenotypeTable' slot reads scattered across
# the package. It validates that the input carries the data the
# calling function needs and warns when the deprecated
# 'TasselGenotypePhenotype' class is used.
#
# @param x
# A 'TasselGenotype', 'TasselPhenotype', 'TasselGenomicDataset', or
# (deprecated) 'TasselGenotypePhenotype' object.
# @param what
# Which data the calling function requires: '"genotype"',
# '"phenotype"', or '"both"'.
# @param fn
# Name of the calling function, used in errors and warnings. Defaults
# to the function that called '.resolveTasselInput()'.
#
# @return
# A 'list' with elements 'jGt' (Java 'GenotypeTable'), 'jPh' (Java
# 'Phenotype'), 'jGp' (Java 'GenotypePhenotype'), and 'original' (the
# input object, for '.wrapLikeInput()'). Components absent from the
# input are returned as Java nulls.
#
# @importFrom rJava is.jnull
.resolveTasselInput <- function(
    x,
    what = c("genotype", "phenotype", "both"),
    fn = NULL
) {
    what <- rlang::arg_match(what)

    if (is.null(fn)) {
        fn <- .callerName()
    }

    if (.isAnyClass(x, TASSEL_INPUT$LEGACY)) {
        lifecycle::deprecate_warn(
            when = "0.14.0",
            what = I(
                sprintf(
                    "Passing a <TasselGenotypePhenotype> object to `%s()`",
                    fn
                )
            ),
            details = c(
                "i" = paste0(
                    "Build inputs with `readGenotype()`, `readPhenotype()`, ",
                    "or `readGenomicDataset()` instead."
                ),
                "i" = paste0(
                    "<TasselGenotypePhenotype> is scheduled for removal in ",
                    "the next major release."
                )
            ),
            # Two frames up is user code, so the warning is reported
            # against the caller of 'fn' rather than rTASSEL itself
            user_env = rlang::caller_env(2)
        )
    } else if (!.isAnyClass(x, TASSEL_INPUT$MODERN)) {
        rlang::abort(c(
            sprintf("Unsupported input object for `%s()`", fn),
            "x" = sprintf(
                "Got an object of class <%s>",
                paste(class(x), collapse = "/")
            ),
            "i" = sprintf(
                "Accepted classes: %s",
                paste0(
                    "<", c(TASSEL_INPUT$MODERN, TASSEL_INPUT$LEGACY), ">",
                    collapse = ", "
                )
            )
        ))
    }

    jComps <- list(
        jGt      = getGenotypeTable(x),
        jPh      = getPhenotypeTable(x),
        jGp      = getGenotypePhenotype(x),
        original = x
    )

    if (what %in% c("genotype", "both") && rJava::is.jnull(jComps$jGt)) {
        rlang::abort(
            sprintf("`%s()` needs genotype data - none found in input", fn)
        )
    }
    if (what %in% c("phenotype", "both") && rJava::is.jnull(jComps$jPh)) {
        rlang::abort(
            sprintf("`%s()` needs phenotype data - none found in input", fn)
        )
    }

    return(jComps)
}


## ----
# Rebuild an R object of the same class as an earlier input
#
# @description
# Implements rTASSEL's "class in, class out" rule: functions that
# return data objects hand their Java result to this helper so that
# modern classes yield modern classes and 'TasselGenotypePhenotype'
# input yields 'TasselGenotypePhenotype' output for back-compatibility.
#
# @param jObj
# A Java object reference returned by TASSEL. The component needed by
# the target class is extracted with the house-keeping getters, so a
# 'GenotypePhenotype' can be handed back for a 'TasselGenotype' input.
# @param original
# The object originally passed to the calling function (element
# 'original' of '.resolveTasselInput()').
#
# @return
# An object of the same class as 'original'.
.wrapLikeInput <- function(jObj, original) {
    if (!methods::is(jObj, "jobjRef") || rJava::is.jnull(jObj)) {
        rlang::abort("`jObj` must be a non-null Java object reference")
    }

    # Ordered most to least specific: 'TasselNumericGenotype' is a
    # 'TasselGenotype' and 'newTasselGenotype()' preserves the subclass.
    if (methods::is(original, "TasselGenotype")) {
        return(newTasselGenotype(getGenotypeTable(jObj), original))
    }
    if (methods::is(original, "TasselPhenotype")) {
        return(createTasselPhenotype(getPhenotypeTable(jObj)))
    }
    if (methods::is(original, "TasselGenomicDataset")) {
        return(createTasselGenomicDataset(getGenotypePhenotype(jObj)))
    }
    if (methods::is(original, "TasselGenotypePhenotype")) {
        return(.tasselObjectConstructor(jObj))
    }

    rlang::abort(sprintf(
        "Cannot rebuild an object of class <%s>",
        paste(class(original), collapse = "/")
    ))
}


## ----
# Name of the function two frames up the call stack
#
# @description
# Gives '.resolveTasselInput()' a sensible default for 'fn' when a
# call site does not name itself.
#
# @return
# A single 'character' value.
.callerName <- function() {
    call <- rlang::caller_call(2)

    if (is.null(call)) {
        return("this function")
    }

    return(rlang::as_label(call[[1]]))
}
