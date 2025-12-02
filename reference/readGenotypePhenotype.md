# Wrapper function of TasselGenotypePhenotype class for GenotypePhenotype combined data

Creates a Java GenotypePhenotype object which is used for
`TasselGenotypePhenotype` object construction. The Java
GenotypePhenotype object is created via an `intersect` method from
TASSEL.

## Usage

``` r
readGenotypePhenotype(genoPathOrObj, phenoPathDFOrObj, ...)
```

## Arguments

- genoPathOrObj:

  a path to a genotype file (e.g. VCF, hmp, etc.) or TASSEL Genotype Obj

- phenoPathDFOrObj:

  a path, a data frame of phenotypic data, or TASSEL Phenotype Obj

- ...:

  Additional parameters to be sent to the function. Currently, if an R
  data frame object is passed, additional parameters will also need to
  be entered for this process. These parameters are derived from the
  `readPhenotypeFromDataFrame` function. Mainly, `taxaID` is required.
  If you would like to specify depth retention and position sorting in
  Genotype Tables from a path, indicate them here. See
  [`readGenotypeTableFromPath()`](https://rtassel.maizegenetics.net/reference/readGenotypeTableFromPath.md)
  for more detail.

## Value

Returns an object of `TasselGenotypePhenotype` class.
