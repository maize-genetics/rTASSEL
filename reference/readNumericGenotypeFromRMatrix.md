# Read Genotype Data from R Matrix

This function constructs a \`TasselNumericGenotype\` object from an R
matrix by interfacing with the TASSEL Java API. It creates taxa lists,
position lists, reference probabilities, and a genotype table using the
provided matrix data.

## Usage

``` r
readNumericGenotypeFromRMatrix(m, asTGP = TRUE)
```

## Arguments

- m:

  A numeric matrix where rows represent taxa and columns represent
  positions/sites. The matrix values are used to calculate reference
  probabilities.

- asTGP:

  Should the return object be a "classic" `TasselGenotypePhenotype`
  object (`TRUE`) or should it return a `TasselNumericGenotype` object
  (`FALSE`)? Defaults to `TRUE`.

## Value

An object of class \`TasselNumericGenotype\` containing the genotype
data and metadata.

## Details

The following components are constructed using respective Java builder
classes:

- Taxa list is created using the `TASSEL_JVM$TAXA_LIST_BUILDER` Java
  class.

- Position list is constructed using the
  `TASSEL_JVM$POSITION_LIST_BUILDER` and
  `TASSEL_JVM$GENERAL_POSITION_BUILDER` Java classes.

- Reference probabilities are built using the
  `TASSEL_JVM$REF_PROBABILITY_BUILDER` Java class.

- The genotype table is created using the
  `TASSEL_JVM$GENOTYPE_TABLE_BUILDER` Java class.
