# Create a QQ plot from rTASSEL association output

This function allows for quick generation of a QQ plot from rTASSEL
association statistical output data.

## Usage

``` r
plotQQ(assocRes, trait = NULL, overlay = TRUE, interactive = FALSE)
```

## Arguments

- assocRes:

  An object of type `AssociationResults`

- trait:

  Which phenotypic trait do you want to plot? If set to `NULL`, this
  will generate a faceted plot with all mapped traits.

- overlay:

  Do you want trait results faceted or overlayed into one single plot?
  Defaults to `TRUE`.

- interactive:

  Do you want to produce an interactive visualization? Defaults to
  `FALSE`.

## Value

Returns a `ggplot2` object
