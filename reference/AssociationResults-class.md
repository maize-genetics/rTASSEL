# AssociationResults Class

Class `AssociationResults` defines a `rTASSEL` Class for storing TASSEL
5 and general-purpose association results from \*WAS studies.

Prints outs summary information from `AssociationResults` objects

## Usage

``` r
# S4 method for class 'AssociationResults'
show(object)
```

## Arguments

- object:

  a `AssociationResults` object

## Slots

- `results`:

  A list of `data.frame` objects containing summary results

- `traits`:

  A vector of type `character` containing the trait IDs modeled.

- `assocType`:

  `character` object describing association type
