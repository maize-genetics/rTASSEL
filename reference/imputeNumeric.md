# Imputation methods in Numerical Transformations

This method takes an input `GenotypeTable` object with missing values
and imputes the missing values using one of the chosen methods.

## Usage

``` r
imputeNumeric(
  tasObj,
  byMean = TRUE,
  nearestNeighbors = 5,
  distance = c("Euclidean", "Manhattan", "Cosine")
)
```

## Arguments

- tasObj:

  an rTASSEL `TasselGenotypePhenotype` object.

- byMean:

  Will imputation be performed by computing the mean of the respective
  column? Defaults to `TRUE`.

- nearestNeighbors:

  Number of nearest neighbors to be evaluated. Defaults to `5`.

- distance:

  Distance type. Options are `Euclidean`, `Manhattan`, or `Cosine`.
  Defaults to `Euclidean`.
