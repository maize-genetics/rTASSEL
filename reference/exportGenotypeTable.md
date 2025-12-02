# Export Genotype Table to Disk

Exports genotype tables to various flat file formats.

## Usage

``` r
exportGenotypeTable(
  tasObj,
  file,
  format = c("vcf", "hapmap", "plink", "flapjack"),
  keepDepth = TRUE,
  taxaAnnotations = TRUE,
  branchLengths = TRUE
)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype` that contains a genotype
  table.

- file:

  Output file name.

- format:

  Export file format. This function current supports the following:

  - `vcf` - A VCF (variant call) file

  - `hapmap` - HapMap files

  - `plink` - Plink files

  - `flapjack` - FlapJack files

- keepDepth:

  Whether to keep depth if format supports depth. Defaults to `TRUE`.

- taxaAnnotations:

  Whether to include taxa annotations if format supports taxa. Defaults
  to `TRUE`.

- branchLengths:

  Whether to include branch lengths for Newick formatted files. Defaults
  to `TRUE`.
