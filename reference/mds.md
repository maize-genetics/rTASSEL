# Run MDS on `TasselDistanceMatrix` objects

Perform multidimensional scaling (MDS) on `TasselDistanceMatrix`
objects.

## Usage

``` r
mds(distMat, nAxes = 5, removeNaN = TRUE)
```

## Arguments

- distMat:

  A `TasselDistanceMatrix` object.

- nAxes:

  The number of axes or dimensions and associated eigenvalues to be
  returned by the analysis. Defaults to `5`.

- removeNaN:

  Remove `NaNs` from matrix before performing MDS. Defaults to `TRUE`.

## Value

A `DataFrame` object.
