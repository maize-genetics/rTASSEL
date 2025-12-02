# Get an R/`DataFrame` phenotype data frame from TASSEL object

This function will extract a `DataFrame`-based R data frame from an
object of `TasselGenotypePhenotype` class that contains phenotypic data.
Column data will be converted to the following types data depending on
TASSEL data type:

- `<taxa>`: `character`

- `<data>`: `numeric`

- `<covariate>`: `numeric`

- `<factor>`: `factor`

## Usage

``` r
getPhenotypeDF(tasObj)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype`.

## Value

Returns an `DataFrame` based data frame
