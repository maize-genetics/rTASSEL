# Create an LDRegion object

Constructor for
[`LDRegion`](https://rtassel.maizegenetics.net/reference/LDRegion-class.md)
objects. Defines a single genomic region for LD block highlighting.

## Usage

``` r
LDRegion(
  start,
  end,
  label = NA_character_,
  color = "black",
  linewidth = NA_real_,
  showSpan = TRUE
)
```

## Arguments

- start:

  Numeric. Start position of the region in base pairs.

- end:

  Numeric. End position of the region in base pairs. Must be `>= start`.

- label:

  Character. Optional text label for the block. Defaults to
  `NA_character_` (no label).

- color:

  Character. Outline color for the block highlight. Must be a valid R
  color. Defaults to `"black"`.

- linewidth:

  Numeric. Outline thickness (mm) for the block highlight. `NA` (the
  default) uses the plot-level auto-calculated width.

- showSpan:

  Logical. If `TRUE` (default), the genomic span (e.g., “342.9 kbp”) is
  appended to the block annotation. When no `label` is provided, the
  span is shown on its own. Set to `FALSE` to suppress the span
  entirely.

## Value

An object of class
[`LDRegion`](https://rtassel.maizegenetics.net/reference/LDRegion-class.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Region with a label (span shown by default)
LDRegion(start = 157104, end = 500000, label = "Block A")

# Region with custom color, no label but span still shown
LDRegion(start = 800000, end = 1200000, color = "blue")

# Suppress the span text
LDRegion(start = 100, end = 500, label = "QTL", showSpan = FALSE)

# Custom line thickness
LDRegion(start = 100, end = 500, linewidth = 1.5)

# Pass a list of regions to plotLD
plotLD(
  tasObj,
  ldBlocks = list(
    LDRegion(start = 157104, end = 500000, label = "A"),
    LDRegion(start = 800000, end = 1200000, label = "B", color = "blue")
  )
)
} # }
```
