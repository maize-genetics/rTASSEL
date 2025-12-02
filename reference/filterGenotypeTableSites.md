# Filter genotype table by sites

This function will filter R objects of `TasselGenotypePhenotype` class
containing genotype tables. The parameters for this function are derived
from TASSEL's `FilterSiteBuilder` plugin.

## Usage

``` r
filterGenotypeTableSites(
  tasObj,
  siteMinCount = 0,
  siteMinAlleleFreq = 0,
  siteMaxAlleleFreq = 1,
  minHeterozygous = 0,
  maxHeterozygous = 1,
  removeMinorSNPStates = FALSE,
  removeSitesWithIndels = FALSE,
  siteRangeFilterType = c("none", "sites", "position"),
  startSite = NULL,
  endSite = NULL,
  startChr = NULL,
  startPos = NULL,
  endChr = NULL,
  endPos = NULL,
  gRangesObj = NULL,
  chrPosFile = NULL,
  bedFile = NULL
)
```

## Arguments

- tasObj:

  An object of class `TasselGenotypePenotype`.

- siteMinCount:

  Site minimum count of alleles not unknown. Can range from 0 to inf.
  Defaults to 0.

- siteMinAlleleFreq:

  Site minimum minor allele frequency. Can range from 0 to 1.0. Defaults
  to 0.0.

- siteMaxAlleleFreq:

  Site maximum minor allele frequency. Can range from 0 to 1.0. Defaults
  to 1.0.

- minHeterozygous:

  Min heterozygous proportion. Can range from 0 to 1.0. Defaults to 0.0.

- maxHeterozygous:

  Max heterozygous proportion. Can range from 0 to 1.0. Defaults to 1.0.

- removeMinorSNPStates:

  Remove minor SNP states. Defaults to `FALSE`.

- removeSitesWithIndels:

  Remove sites containing an indel (`+` or `-`). Defaults to `FALSE`.

- siteRangeFilterType:

  True if filtering by site numbers. False if filtering by chromosome
  and position. Options are `none`, `sites`, or `position`. Defaults to
  `none`.

- startSite:

  The start site. Defaults to 0.

- endSite:

  The end site. Defaults to 0.

- startChr:

  Start chromosome for site filtration range if `position` is chosen
  from `siteRangeFilterType`. Needs end chromosome (`endChr`) to work.

- startPos:

  Physical start position (bp) for filtration range if `position` is
  chosen from `siteRangeFilterType`. If `NULL`, the first physical
  position in the data set will be chosen.

- endChr:

  End chromosome for site filtration range if `position` is chosen from
  `siteRangeFilterType`. Needs start chromosome (`endChr`) to work.

- endPos:

  Physical end position (bp) for filtration range if `position` is
  chosen from `siteRangeFilterType`. If `NULL`, the last physical
  position in the data set will be chosen.

- gRangesObj:

  Filter genotype table by `GenomicRanges` object. If this parameter is
  selected, you cannot utilize the parameters, `chrPosFile` or
  `bedFile`. Defaults to `NULL`.

- chrPosFile:

  An optional chromosome position file path of `character` class.
  Defaults to `NULL`. **Note:** a chromosome position file must contain
  correct formatting (e.g. a two column file with the header of
  `c("Chromosome", "Position")`).

- bedFile:

  An optional BED coordinate file path of `character` class. Defaults to
  `NULL`.

## Value

Returns an object of `TasselGenotypePhenotype` class.
