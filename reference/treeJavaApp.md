# R interface for Archaeopteryx interactive tree viewer

This function acts as a wrapper for TASSEL's interface to the
Archaeopteryx Java tree Viewer.

## Usage

``` r
treeJavaApp(tasObj, clustMethod = c("Neighbor_Joining", "UPGMA"))
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype`.

- clustMethod:

  What clustering method should be used? Current options are `UGMA` and
  `Neighbor_Joining`. Defaults to `Neighbor_Joining`.

## Value

Returns a Java-based visualization application.
