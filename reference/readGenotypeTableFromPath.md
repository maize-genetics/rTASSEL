# Wrapper function of TasselGenotypePhenotype class for genotype data

This function is a wrapper for the `TasselGenotypePhenotype` class. It
is used for storing genotype information into a class object.

## Usage

``` r
readGenotypeTableFromPath(path, keepDepth = FALSE, sortPositions = FALSE)
```

## Arguments

- path:

  A genotype data path (e.g. `*.VCF, *.hmp`, etc.).

- keepDepth:

  Should depth be kept? Defaults to `FALSE`.

- sortPositions:

  Should positions be sorted? Defaults to `FALSE`.

## Value

Returns an object of `TasselGenotypePhenotype` class.
