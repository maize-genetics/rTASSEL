# Linkage disequilibrium plot

Calculates linkage disequilibrium (LD) and generates a static plot using
`ggplot2` graphics.

## Usage

``` r
ldPlot(
  tasObj,
  ldType = c("All", "SlidingWindow"),
  windowSize = NULL,
  hetCalls = c("missing", "ignore", "third"),
  plotVal = c("r2", "DPrime", "pDiseq"),
  verbose = TRUE
)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype`. That contains a genotype
  table.

- ldType:

  How do you want LD calculated? Currently, the available options are
  `"All"` and `"SlidingWindow"`. If `All` is selected, LD will be
  calculated for every combination of sites in the alignment (NOTE: this
  may produce a massive series of combinations; use only on heavily
  filtered genotype tables). If `SlidingWindow` is selected, LD will be
  calculated for sites within a window of sites surrounding the current
  site. Defaults to `"All"`.

- windowSize:

  What size do you want your LD analysis window? If you have chosen
  `SlidingWindow` for the `ldType` parameter, you will need to specify
  window size.

- hetCalls:

  How should heterozygous calls be handled? Current options are
  `"ignore"` (ignore altogether), `"missing"` (set to missing), and
  `"third"` (treat as third state).

- plotVal:

  What LD value do you want to plot? Options are:

  - `r2`: \\r^{2}\\ (Default parameter)

  - `DPrime`: \\D'\\

  - `pDiseq`: *p*-value

- verbose:

  Display messages? Defaults to `TRUE`.

## Value

Returns a `ggplot2` object.

## Details

Linkage disequilibrium between any set of polymorphisms can be estimated
by initially filtering a genotype dataset and then using this function.
At this time, \\D'\\, \\r^{2}\\ and P-values will be estimated. The
current version calculates LD between haplotypes with known phase only
(unphased diploid genotypes are not supported; see PowerMarker or
Arlequin for genotype support).

- \\D'\\ is the standardized disequilibrium coefficient, a useful
  statistic for determining whether recombination or homoplasy has
  occurred between a pair of alleles.

- \\r^{2}\\ represents the correlation between alleles at two loci,
  which is informative for evaluating the resolution of association
  approaches.

\\D'\\ and \\r^{2}\\ can be calculated when only two alleles are
present. If more than two alleles, only the two most frequent alleles
are used. P-values are determined by a two-sided Fisher's Exact test is
calculated. Since LD is meaningless when scored with very small sample
sizes, a minimum of 20 taxa must be present to calculate LD and there
must be 2 or more minor alleles.

## See also

[`linkageDiseq`](https://rtassel.maizegenetics.net/reference/linkageDiseq.md),
[`ldJavaApp`](https://rtassel.maizegenetics.net/reference/ldJavaApp.md)
