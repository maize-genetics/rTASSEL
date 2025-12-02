# R interface for TASSEL's association methods

This function acts as a front-end for TASSEL's extensive association
analysis methods. Using this function, users can run the following
TASSEL association methods:

- best linear unbiased estimates (BLUEs)

- generalized linear model (GLM)

- mixed linear model

- Fast association (Shabalin 2012)

## Usage

``` r
assocModelFitter(
  tasObj,
  formula,
  fitMarkers = FALSE,
  kinship = NULL,
  fastAssociation = FALSE,
  maxP = 0.001,
  maxThreads = 1,
  minClassSize = 0,
  outputFile = NULL,
  biallelicOnly = FALSE,
  appendAddDom = FALSE
)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype`.

- formula:

  An R-based linear model formula. The general layout of this formula
  uses the following TASSEL data scheme:
  `<data> ~ <factor> and/or <covariate>`. If all traits in a Phenotype
  object should be ran, a simplified formula (`. ~ .`) can be used. This
  scheme can also be used for running all `<data>` or `<factor>` and/or
  `<covariate>` data as well. Single variables are separated witha `+`
  operator. See vignette for further clarification.

- fitMarkers:

  Should marker data be fitted? If `TRUE`, GLM analysis will be
  executed. If `FALSE`, BLUEs will be calculated. Defaults to `FALSE`.

- kinship:

  Should kinship data be accounted for in the model? If so, a TASSEL
  kinship matrix object of class `TasselDistanceMatrix` must be
  submitted. Defaults to `NULL`.

- fastAssociation:

  Should TASSEL's Fast Association plugin be used? Consider setting to
  `TRUE` if you have many phenotypes in your data set.

- maxP:

  Maximum p-value (0 - 1) to be reported. Currently works with fast
  association only. Defaults to a p-value of `0.001` will be used as a
  threshold. **Note:** p-value parameter will not be used for BLUE
  analysis.

- maxThreads:

  Maximum threads to be used when running fast association. If `NULL`,
  all threads on machine will be used.

- minClassSize:

  The minimum acceptable genotype class size. Genotypes in a class with
  a smaller size will be set to missing. Defaults to 0.

- outputFile:

  Output file prefix to be specified in case you want to write data
  directly to disk. Highly recommended for large datasets. If `NULL`, no
  data will be saved to disk. If a character

- biallelicOnly:

  Only test sites that are bi-allelic. The alternative is to test sites
  with two or more alleles. Defaults to `FALSE`

- appendAddDom:

  If true, additive and dominance effect estimates will be added to the
  stats report for bi-allelic sites only. The effect will only be
  estimated when the data source is genotype (not a probability). The
  additive effect will always be non-negative. Defaults to `FALSE`.

## Value

Returns an R list containing `DataFrame`-based data frames
