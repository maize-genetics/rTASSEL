# Filter genotype table by taxa

This function will filter R objects of `TasselGenotypePhenotype` class
containing genotype tables. The parameters for this function are derived
from TASSEL's `FilterTaxaBuilder` plugin.

## Usage

``` r
filterGenotypeTableTaxa(
  tasObj,
  minNotMissing = 0,
  minHeterozygous = 0,
  maxHeterozygous = 1,
  taxa = NULL
)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype`.

- minNotMissing:

  Minimum proportion of sites not unknown to pass this filter. Value can
  be between 0.0 and 1.0.

- minHeterozygous:

  Minimum proportion of sites that are heterozygous. Value can be
  between 0.0 and 1.0.

- maxHeterozygous:

  Maximum proportion of sites that are heterozygous. Value can be
  between 0.0 and 1.0.

- taxa:

  Vector of taxa IDs (as character) to subset. Defaults to `NULL`.

## Value

Returns an object of `TasselGenotypePhenotype` class.
