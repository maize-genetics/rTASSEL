# Stepwise Model Fitter

Fits a model using stepwise selection criteria based on the provided
genotype-phenotype object.

## Usage

``` r
stepwiseModelFitter(
  tasObj,
  formula = . ~ .,
  modelType = c("pvalue", "bic", "mbic", "aic"),
  entryLimit = 0.01,
  exitLimit = 0.01,
  maxNumberOfMarkers = 20,
  nPermutations = 0
)
```

## Arguments

- tasObj:

  A `TasselGenotypePhenotype` object. The input data for model fitting.

- formula:

  A model formula, default is `. ~ .`.

- modelType:

  Character. The model selection criteria used to determine which terms
  enter the model and how many. Must be one of `"pvalue"`, `"bic"`,
  `"mbic"`, or `"aic"`. Default is `"pvalue"`.

- entryLimit:

  Numeric. The enter limit or maximum p-value for which a term can enter
  the model. Must be in range 0.0–1.0. Default is `0.01`.

- exitLimit:

  Numeric. A term exits the model on a backward step if its p-value is
  greater than this value. Must be in range 0.0–1.0. Default is `0.01`.

- maxNumberOfMarkers:

  Integer. The maximum number of markers that will be fit, if the enter
  limit is not reached first. Range 0–10000. Default is `20`.

- nPermutations:

  Integer. Number of permutations for the model to determine an
  empirical alpha. Range 0–100000. Default is `0`.

## Value

A list of association result tables including ANOVA reports and marker
effect estimates, with and without confidence intervals.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- stepwiseModelFitter(tasObj)
} # }
```
