# Run PCA on Genotype Table

This method performs principal components analysis and returns the
requested number of PC axes (components).

## Usage

``` r
pca(
  tasObj,
  useCovariance = TRUE,
  limitBy = c("number_of_components", "min_eigenvalue", "total_variance"),
  nComponents = 5,
  minEigenval = 0,
  totalVar = 0.5,
  reportEigenvalues = TRUE,
  reportEigenvectors = TRUE
)
```

## Arguments

- tasObj:

  an rTASSEL `TasselGenotypePhenotype` object.

- useCovariance:

  If `TRUE`, analysis will do an eigenvalue decomposition of the
  covariance matrix. If `FALSE`, it will use a correlation matrix. NOTE:
  Using the covariance matrix is recommended for genotypes while the
  correlation matrix is often used for phenotypes. Defaults to `TRUE`.

- limitBy:

  This parameter determines the type of value that will be used to limit
  the number of principal components (axes) returned. The possible
  choices are `number_of_components`, `min_eigenvalue`, and
  `total_variance`.

- nComponents:

  The analysis will return this many principal components up to the
  number of taxa.

- minEigenval:

  All principal components with an eigenvalue greater than or equal to
  this value will be returned. NOTE: works only if `min_eigenvalue` is
  set in the `limitBy` parameter.

- totalVar:

  The first principal components that together explain this proportion
  of the total variance will be returned. NOTE: works only if
  `total_variance` is set in the `limitBy` parameter.

- reportEigenvalues:

  Returns a list of eigenvalues sorted high to low.

- reportEigenvectors:

  Returns the eigenvectors calculated from a Singular Value
  Decomposition of the data. The resulting table can be quite large if
  the number of variants and taxa are big.

## Value

A `DataFrame` object.
