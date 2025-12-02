# Read Genotype Data

This function reads genotype data from a file path or an R matrix. It
supports optional sorting of positions and retaining depth information.

## Usage

``` r
readGenotype(x, sortPositions = FALSE, keepDepth = FALSE)
```

## Arguments

- x:

  A character string representing the file path to the genotype data or
  a matrix containing genotype data.

- sortPositions:

  A logical value indicating whether to sort positions in the genotype
  data. Default is `FALSE`.

- keepDepth:

  A logical value indicating whether to retain depth information in the
  genotype data. Default is `FALSE`.

## Value

A processed genotype object based on the input data.

## Details

- If `x` is a character string:

  - The function checks if the file exists.

  - Reads the genotype data from the file path using
    `readGenotypeFromPath`.

- If `x` is a matrix:

  - The function processes the genotype data using
    `readGenotypeFromRMatrix`.

- If `x` is neither a character string nor a matrix:

  - An error is raised.

## Examples

``` r
if (FALSE) { # \dontrun{
# Read genotype data from a file
readGenotype("path/to/genotype/file.txt", sortPositions = TRUE, keepDepth = TRUE)

# Read genotype data from a matrix
genotypeMatrix <- matrix(data = ..., nrow = ..., ncol = ...)
readGenotype(genotypeMatrix)
} # }
```
