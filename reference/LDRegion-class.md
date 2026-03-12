# LDRegion Class

An S4 class to define a genomic region for highlighting LD blocks on a
[`plotLD`](https://rtassel.maizegenetics.net/reference/plotLD.md) plot.
Each region is specified by a start and end position (in base pairs)
with an optional text label, outline color, line thickness, and span
display control.

Prints a compact summary of an `LDRegion` object.

## Usage

``` r
# S4 method for class 'LDRegion'
show(object)
```

## Arguments

- object:

  An `LDRegion` object.

## Slots

- `start`:

  A single numeric value for the region start position (bp).

- `end`:

  A single numeric value for the region end position (bp). Must be
  greater than or equal to `start`.

- `label`:

  A single character string used as the block annotation. Defaults to
  `NA_character_` (no label).

- `color`:

  A single character string specifying the outline color for the block
  highlight. Must be a valid R color. Defaults to `"black"`.

- `linewidth`:

  A single numeric value controlling the outline thickness (in mm) for
  the block highlight. Defaults to `NA_real_`, which uses the plot-level
  auto-calculated width.

- `showSpan`:

  Logical. If `TRUE`, the genomic span (e.g., “342.9 kbp”) is shown in
  the block annotation regardless of whether a `label` is provided. If
  `FALSE`, the span is never shown. Defaults to `TRUE`.

## See also

[`plotLD`](https://rtassel.maizegenetics.net/reference/plotLD.md)
