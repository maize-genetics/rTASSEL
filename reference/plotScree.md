# Generate scree plots

This function will generate a line plot of eigenvalues of PCs. These
plots can be used to determine the number of factors to retain for
downstream analyses (e.g GWAS).

## Usage

``` r
plotScree(
  pcaObj,
  nComp = 10,
  interactive = FALSE,
  lineColor = "grey",
  pointColor = "black"
)
```

## Arguments

- pcaObj:

  A `PCAResults` object containing eigenvalue data ( e.g
  `Eigenvalues_Datum`).

- nComp:

  Total number of principal components to plot. Defaults to `10`.

- interactive:

  Do you want to produce an interactive visualization? Defaults to
  `FALSE`.

- lineColor:

  Color of scree line.

- pointColor:

  Color of scree points.
