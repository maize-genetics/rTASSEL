# Create a TASSEL kinship matrix

This function will calculate a kinship matrix using TASSEL's `Kinship`
plugin and respective parameters.

## Usage

``` r
kinshipMatrix(
  tasObj,
  method = "Centered_IBS",
  maxAlleles = 6,
  algorithmVariation = "Observed_Allele_Freq"
)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype`.

- method:

  A Kinship method. Defaults to `Centered_IBS`. Other options include
  `Normalized_IBS`, `Dominance_Centered_IBS`, and
  `Dominance_Normalized_IBS`.

- maxAlleles:

  Maximum number of alleles. Can be within the range of `2` to `6`.

- algorithmVariation:

  Algorithm variation. If `Dominance_Centered_IBS` is selected, users
  can switch between `Observed_Allele_Freq` and
  `Proportion_Heterozygous`.

## Value

Returns a \`TasselDistanceMatrix\` object.
