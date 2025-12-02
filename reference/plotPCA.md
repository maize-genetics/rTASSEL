# Generate PCA plot

This function will generate a general 2D scatterplot for a given set of
principal components.

## Usage

``` r
plotPCA(
  pcaObj,
  x = 1,
  y = 2,
  cluster = FALSE,
  nClust = 2,
  metadata = NULL,
  mCol = NULL,
  interactive = FALSE,
  pointColor = "black"
)
```

## Arguments

- pcaObj:

  A `PCAResults` object containing eigenvalue data ( e.g `PC_Datum`).

- x:

  Principal component number to plot on x-axis.

- y:

  Principal component number to plot on y-axis.

- cluster:

  Plot points by hierachical clustering groups? Defaults to `FALSE`.

- nClust:

  Number of clusters to report.

- metadata:

  A `data.frame` object of additional categorical information for each
  sample/taxa. Defaults to `NULL`.

- mCol:

  What metadata column do you want to plot? Needs metadata object to
  work.

- interactive:

  Do you want to produce an interactive visualization? Defaults to
  `FALSE`.

- pointColor:

  color of non-categorical PCA data points.
