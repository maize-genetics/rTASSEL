# Filter genotype table by site IDs

Filter a genotype table object by specifying literal site names (IDs)
for variant markers.

## Usage

``` r
filterGenotypeTableBySiteName(tasObj, siteNames)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype`.

- siteNames:

  A character vector of site names to filter on.

## Value

Returns an object of `TasselGenotypePhenotype` class.
