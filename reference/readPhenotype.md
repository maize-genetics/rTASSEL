# Read and convert phenotype data into TASSEL 5 phenotype objects

This function reads phenotype data from either a file path or a data
frame.

## Usage

``` r
readPhenotype(x, attr = NULL)
```

## Arguments

- x:

  A character string representing the file path to the phenotype data or
  a data frame containing the phenotype data.

- attr:

  An optional attribute metadata parameter required when `x` is a data
  frame. Defaults to `NULL`.

## Value

A phenotype object created from the input data.

## Details

- If `x` is a character string, the function assumes it is a file path
  and calls `readPhenotypeFromFile(x)`.

- If `x` is a data frame, the function requires the `attr` parameter to
  provide metadata and calls `readPhenotypeFromDf(x, attr)`.

- If `x` is neither a character string nor a data frame, the function
  throws an error.

## Examples

``` r
if (FALSE) { # \dontrun{
# Reading phenotype data from a file
phenotype <- readPhenotype("path/to/phenotype/file.txt")

# Reading phenotype data from a data frame
attrDf <- tibble::tribble(
    ~"col_id",      ~"tassel_attr",
    "taxa_id",      "taxa",
    "plant_height", "data",
    "PC1",          "covariate",
    "yield",        "data",
)
df <- tibble::tribble(
    ~"taxa_id", ~"plant_height", ~"PC1", ~"yield",
    "line_a",   12.3,            0.5,    2,
    "line_b",   22.8,            -1.5,   3,
)

phenotypeDf <- readPhenotype(df, attr = attrDf)
} # }
```
