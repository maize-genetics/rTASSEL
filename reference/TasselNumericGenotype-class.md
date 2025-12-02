# TasselNumericGenotype Class Definition

Defines the \`TasselNumericGenotype\` class, which extends the
\`TasselGenotype\` class. This class is used to represent numeric
genotype data in the TASSEL 5 framework.

This method is used to display information about a
`TasselNumericGenotype` object. It prints a summary of the genotype
data, including the number of taxa, number of sites, and memory address
of the Java object.

## Usage

``` r
# S4 method for class 'TasselNumericGenotype'
show(object)
```

## Arguments

- object:

  An object of class `TasselNumericGenotype`.

## Details

The function `printNumGtDisp` is called internally to format and display
the genotype data. The number of taxa and sites are retrieved from the
Java reference object associated with the `TasselNumericGenotype`
instance.
