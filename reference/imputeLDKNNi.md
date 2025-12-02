# LD KNNi imputation

This imputation algorithm uses LD to identify good predictors for each
SNP, and then uses the high LD SNPs to identify K-Nearest Neighbors. The
genotype is called with a weighted mode of the KNNs.

## Usage

``` r
imputeLDKNNi(tasObj, highLDSSites = 30, knnTaxa = 10, maxDistance = 1e+07)
```

## Arguments

- tasObj:

  an rTASSEL `TasselGenotypePhenotype` object.

- highLDSSites:

  Number of sites in high LD to use in imputation. Acceptable values are
  between `2` and `2000`. Defaults to `30`.

- knnTaxa:

  Number of neighbors to use in imputation. Acceptable values are
  between `2` and `200`. Defaults to `10`.

- maxDistance:

  Maximum physical distance between sites to search for LD (-1 for no
  distance cutoff - unlinked chromosomes will be tested). Defaults to
  `10e6`.
