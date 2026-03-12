# LDResults Class

An S4 class that stores the results of a linkage disequilibrium (LD)
analysis produced by
[`linkageDiseq`](https://rtassel.maizegenetics.net/reference/linkageDiseq.md).
The object contains the full pairwise LD table together with the
parameters used to generate it.

Prints a compact summary of an `LDResults` object.

## Usage

``` r
# S4 method for class 'LDResults'
show(object)
```

## Arguments

- object:

  An `LDResults` object.

## Slots

- `results`:

  A `data.frame` containing the pairwise LD statistics. Expected columns
  include `Locus1`, `Position1`, `Locus2`, `Position2`, `R^2`, `DPrime`,
  `pDiseq`, and `N`, among others.

- `ldType`:

  A single character string indicating the LD calculation type: `"All"`
  or `"SlidingWindow"`.

- `windowSize`:

  A single numeric value for the sliding-window size. `NA_real_` when
  `ldType` is `"All"`.

- `hetCalls`:

  A single character string indicating how heterozygous calls were
  handled: `"missing"`, `"ignore"`, or `"third"`.

## See also

[`linkageDiseq`](https://rtassel.maizegenetics.net/reference/linkageDiseq.md),
[`plotLD`](https://rtassel.maizegenetics.net/reference/plotLD.md)
