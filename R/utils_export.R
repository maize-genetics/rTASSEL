## ----
#' @title Export Genotype Table to Disk
#'
#' @description Exports genotype tables to various flat file formats.
#'
#' @param tasObj An object of class \code{TasselGenotypePenotype} that
#'   contains a genotype table.
#' @param file Output file name.
#' @param format Export file format. This function current supports the
#'   following:
#'
#'   \itemize{
#'     \item \code{vcf} - A VCF (variant call) file
#'     \item \code{hapmap} - HapMap files
#'     \item \code{plink} - Plink files
#'     \item \code{flapjack} - FlapJack files
#'   }
#'
#' @param keepDepth Whether to keep depth if format supports depth. Defaults
#'   to \code{TRUE}.
#' @param taxaAnnotations Whether to include taxa annotations if format
#'   supports taxa. Defaults to \code{TRUE}.
#' @param branchLengths Whether to include branch lengths for Newick formatted
#'   files. Defaults to \code{TRUE}.
#'
#' @importFrom rJava .jchar
#' @importFrom rJava is.jnull
#' @importFrom rJava J
#'
#' @export
exportGenotypeTable <- function(tasObj,
                                file,
                                format = c("vcf", "hapmap", "plink", "flapjack"),
                                keepDepth = TRUE,
                                taxaAnnotations = TRUE,
                                branchLengths = TRUE) {
    if (class(tasObj) != "TasselGenotypePhenotype") {
        stop("`tasObj` must be of class `TasselGenotypePhenotype`")
    }

    jGenoTable <- getGenotypeTable(tasObj)
    if (rJava::is.jnull(jGenoTable)) {
        stop("TASSEL genotype object not found")
    }

    if (file == "") {
        stop("File name not specified.")
    }

    # Filter type selection
    format <- match.arg(format)

    rJC <- rJava::J("net.maizegenetics.dna.snp.ExportUtils")

    if (format == "vcf") {
        rJC$writeToVCF(
            jGenoTable,
            file,
            keepDepth, # <- keep depth
            NULL  # <- progress listener
        )
    } else if (format == "hapmap") {
        rJC$writeToHapmap(
            jGenoTable,
            FALSE, # <- diploid?
            file,
            rJava::.jchar(9), # <- tab delimited UTF char
            taxaAnnotations,
            NULL
        )
    } else if (format == "plink") {
        rJC$writeToPlink(
            jGenoTable,
            file,
            rJava::.jchar(9) # <- tab delimited UTF char
        )
    } else if (format == "flapjack") {
        rJC <- rJava::J("net.maizegenetics.plugindef.GenerateRCode")
        rJC$exportToFlapjack(
            jGenoTable,
            file
        )
    }
}


