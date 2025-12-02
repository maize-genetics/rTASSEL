# R interface for TASSEL's tree creation methods

This function acts as a wrapper for TASSEL's `CreateTreePlugin`.

## Usage

``` r
createTree(tasObj, clustMethod = c("Neighbor_Joining", "UPGMA"))
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype`.

- clustMethod:

  What clustering method should be used? Current options are `UGMA` and
  `Neighbor_Joining`. Defaults to `Neighbor_Joining`.

## Value

Returns a `phylo` tree object. See the
[ape](https://cran.r-project.org/web/packages/ape/ape.pdf) package for
further details.
