# PCAResults Class

Class `PCAResults` defines a `rTASSEL` Class for storing TASSEL 5 PCA
results

Prints outs summary information from `PCAResults` objects

## Usage

``` r
# S4 method for class 'PCAResults'
show(object)
```

## Arguments

- object:

  a `PCAResults` object

## Slots

- `results`:

  A list of `data.frame` objects containing summary results

- `jObj`:

  An rJava reference object pointing to PCA results in Java memory
