# Create a QC Manhattan plots from rTASSEL association output

This function allows for quick generation of a QC plot from rTASSEL
association statistical output data. The main goal of this function is
to provide "zoomed" in Manhattan plots that typically fall within
genomic ranges of interest plus a flanking window up- and downstream of
the ranges.

## Usage

``` r
plotManhattanQC(
  assocRes,
  trait = NULL,
  gr,
  window = 1e+05,
  threshold = NULL,
  classicNames = FALSE,
  interactive = FALSE,
  verbose = TRUE
)
```

## Arguments

- assocRes:

  An object of type `AssociationResults`

- trait:

  Which phenotypic trait do you want to plot? If set to `NULL`, this
  will generate a faceted plot with all mapped traits.

- gr:

  Genomic ranges of interest. Can be passed as a `GRanges` object or a
  `data.frame` object.

- window:

  Window size (base-pairs) flanking surround reference range. Defaults
  to `100000` (100,000 base-pairs).

- threshold:

  User-defined \\-log\_{10}(p)\\-value threshold for significant marker
  determination. Once specified any marker points higher than this line
  will be highlighted.

- classicNames:

  Do you want to plot classical gene names instead? NOTE: this will need
  a `classical_id` column for the "ranges of interest" data. Defaults to
  `FALSE`.

- interactive:

  Do you want to produce an interactive visualization? Defaults to
  `FALSE`.

- verbose:

  Should messages be printed to console? Defaults to `FALSE`.

## Value

Returns a `ggplot2` object
