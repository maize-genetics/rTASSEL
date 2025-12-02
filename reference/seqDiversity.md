# Calculate sequence diversity

This performs basic diversity analyses on genetic data, specifically,
average pairwise divergence (\\\pi\\), segregating sites, estimated
mutation rate (\\\theta\\, or \\4N\mu\\), and Tajima's D can be
calculated, as well as sliding windows of diversity.

## Usage

``` r
seqDiversity(
  tasObj,
  startSite = 0,
  endSite = NULL,
  slidingWindowAnalysis = FALSE,
  stepSize = 100,
  windowSize = 500
)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype` that contains a genotype
  table.

- startSite:

  Start site.

- endSite:

  End site. Defaults to the maximum index of markers.

- slidingWindowAnalysis:

  Do you want to analyze diversity in a sliding window? Defaults to
  `FALSE`.

- stepSize:

  Step size.

- windowSize:

  Window size.
