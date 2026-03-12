# Plot SNP density across chromosomes

Generates a heatmap-style visualization of SNP density across all
chromosomes in a genotype table. Each chromosome is displayed on the
Y-axis with genomic position on the X-axis. The number of SNPs within
each window determines tile color intensity using the Viridis palette.

## Usage

``` r
plotSnpDensity(
  tasObj,
  windowSize = 1e+06,
  colorOption = c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako",
    "turbo"),
  logNorm = FALSE,
  interactive = FALSE
)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePhenotype` or `TasselGenotype` that
  contains a genotype table.

- windowSize:

  Size of the genomic window (in base pairs) used to bin SNPs for
  density calculation. Defaults to `1e6` (1 Mb).

- colorOption:

  Which viridis color palette to use? Options are: `"viridis"`
  (default), `"magma"`, `"inferno"`, `"plasma"`, `"cividis"`,
  `"rocket"`, `"mako"`, and `"turbo"`.

- logNorm:

  Should SNP counts be log\\\_{10}\\-transformed before mapping to fill
  color? Useful when a few windows have very high counts that compress
  the color scale. Defaults to `FALSE`.

- interactive:

  Do you want to produce an interactive visualization? Defaults to
  `FALSE`.

## Value

Returns a `ggplot2` object or a `plotly` object if `interactive = TRUE`.
