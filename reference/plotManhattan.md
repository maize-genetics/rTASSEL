# Create a Manhattan plot from rTASSEL association output

This function allows for quick generation of a Manhattan plot from
rTASSEL association statistical output data.

## Usage

``` r
plotManhattan(
  assocRes,
  trait = NULL,
  threshold = NULL,
  colors = c("#91baff", "#3e619b"),
  interactive = FALSE,
  pltTheme = c("default", "classic")
)
```

## Arguments

- assocRes:

  An object of type `AssociationResults`

- trait:

  Which phenotypic trait do you want to plot? If set to `NULL`, this
  will generate a faceted plot with all mapped traits

- threshold:

  User-defined \\-log\_{10}(p)\\-value threshold for significant marker
  determination. Once specified any marker points higher than this line
  will be highlighted.

- colors:

  A vector of `character` colors used for differentiating multiple
  chromosomes. Defaults to 2 shades of blue.

- interactive:

  Do you want to produce an interactive visualization? Defaults to
  `FALSE`.

- pltTheme:

  What theme would like to display for the plot? Only supports one theme
  currently.

## Value

Returns a `ggplot2` object
