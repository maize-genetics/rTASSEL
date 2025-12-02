# Wrapper function of TasselGenotypePhenotype class for phenotype data from an R data frame

This function is a wrapper for the `TasselGenotypePhenotype` class. It
is used for storing phenotype information into a class object.

## Usage

``` r
readPhenotypeFromDataFrame(phenotypeDF, taxaID, attributeTypes = NULL)
```

## Arguments

- phenotypeDF:

  A `R` object of class `data.frame`.

- taxaID:

  The column name that represents your taxa data as a string.

- attributeTypes:

  A vector of non-taxa attributes. If `NULL`, all attributes will be
  TASSEL `<data>` types.

## Value

Returns an object of `TasselGenotypePhenotype` class.
