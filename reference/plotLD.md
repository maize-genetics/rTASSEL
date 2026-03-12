# Linkage disequilibrium plot

Generates a static LD heatmap from a pre-computed
[`LDResults`](https://rtassel.maizegenetics.net/reference/LDResults-class.md)
object using `ggplot2` graphics.

## Usage

``` r
plotLD(
  ldObj,
  plotVal = c("r2", "DPrime", "pDiseq"),
  colorScheme = c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako",
    "turbo", "haploview"),
  ldBlocks = NULL,
  genomicTrack = FALSE,
  showIndex = TRUE,
  verbose = TRUE
)
```

## Arguments

- ldObj:

  An object of class
  [`LDResults`](https://rtassel.maizegenetics.net/reference/LDResults-class.md)
  as returned by
  [`linkageDiseq`](https://rtassel.maizegenetics.net/reference/linkageDiseq.md).

- plotVal:

  What LD value do you want to plot? Options are:

  - `r2`: \\r^{2}\\ (Default parameter)

  - `DPrime`: \\D'\\

  - `pDiseq`: *p*-value

- colorScheme:

  Color palette for the heatmap cells. `"haploview"` uses the classic
  Haploview scheme where cell color is determined by both LOD score and
  \\D'\\: LOD \< 2 and \\D'\\ \< 1 gives white, LOD \< 2 and \\D'\\ = 1
  gives blue, LOD \\\ge\\ 2 and \\D'\\ \< 1 gives shades of pink/red
  scaled by \\D'\\, and LOD \\\ge\\ 2 and \\D'\\ = 1 gives bright red.
  All other options use continuous viridis-family palettes: `"viridis"`
  (default), `"magma"`, `"inferno"`, `"plasma"`, `"cividis"`,
  `"rocket"`, `"mako"`, and `"turbo"`.

- ldBlocks:

  Optional specification of genomic regions to highlight as LD blocks on
  the plot. Accepted formats:

  - A
    [`LDRegion`](https://rtassel.maizegenetics.net/reference/LDRegion-class.md)
    object (single block) or a `list` of `LDRegion` objects (multiple
    blocks). Each region carries `start`, `end`, an optional `label`,
    and an outline `color`. The chromosome is inferred from the data
    (single-chromosome only).

  - A `GRanges` object where each range defines a chromosome and
    start/end positions. An optional `label` metadata column (in
    `mcols`) adds text annotations. Block outlines default to black.

  Sites falling within each region are outlined with a triangular
  border. Defaults to `NULL`.

- genomicTrack:

  Logical. If `TRUE`, a horizontal genomic track is drawn above the LD
  plot showing the approximate physical positions of each site along the
  chromosome. Site IDs are placed at their physical positions on the
  track, and line segments connect each physical position to the
  corresponding evenly-spaced position above the LD triangle. When data
  span multiple chromosomes, positions are displayed cumulatively.
  Defaults to `FALSE`.

- showIndex:

  Logical. If `TRUE` (default), numeric index labels (1, 2, 3, ...) are
  drawn along the diagonal of the LD triangle. Set to `FALSE` to hide
  them.

- verbose:

  Display messages? Defaults to `TRUE`.

## Value

A `ggplot` object.

## Details

This function visualises a pre-computed LD result set. Use
[`linkageDiseq`](https://rtassel.maizegenetics.net/reference/linkageDiseq.md)
to calculate LD first, then pass the resulting
[`LDResults`](https://rtassel.maizegenetics.net/reference/LDResults-class.md)
object here.

## See also

[`linkageDiseq`](https://rtassel.maizegenetics.net/reference/linkageDiseq.md),
[`LDResults`](https://rtassel.maizegenetics.net/reference/LDResults-class.md)
