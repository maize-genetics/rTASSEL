# Create Summarized Experiment from a TASSEL Genotype Table

This function will generate an object of `SummarizedExperiment` class
for marker data derived from a `TasselGenotypePhenotype` class object.

## Usage

``` r
getSumExpFromGenotypeTable(tasObj, coerceDosageToInt = TRUE, verbose = FALSE)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype`.

- coerceDosageToInt:

  Should dosage array be returned as `integer` values? If `FALSE`,
  dosage array will be returned as type `raw` byte values. Returning
  `raw` byte values. Will greatly save on memory. Defaults to `TRUE`.

- verbose:

  Should messages be displayed to console? Defaults to `FALSE`.

## Value

Returns a `SummarizedExperiment` of TASSEL genotype data.
