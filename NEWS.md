# rTASSEL 0.9.19
* Added two new parameters to `filterGenotypeTableSites()`
  + `removeMinorSNPStates`: Boolean; removes minor SNP states.
  + `removeSitesWithIndels`: Boolean; removes sites with indels.
* Added better descriptions to exceptions when users would enter chromosome IDs
  not present in a genotype table.
* Fixed `siteRangeFilterType` parameter bug in `filterGenotypeTableSites()`. 
  Defaults to `none` when user does not specify filter type.


# rTASSEL 0.9.18
* Added functions to calculate linkage disequilibrium (LD)
* Proposed LD functions:
  + `linkageDiseq()` - Returns TASSEL LD table report as data frame
  + `ldPlot()` - Returns static `ggplot2` plot
  + `ldJavaApp()` - Initiates TASSEL's interactive LD viewer


# rTASSEL 0.9.17
* Added new function:
  + `manhattanPlot()`
* Removed tidyverse dependencies


# rTASSEL 0.9.16
* Added write to file parameters for:
  + assocModelFitter()
* Added p-value threshold parameters for:
  + assocModelFitter()
* Added thread usage paramters for:
  + assocModelFitter()
* Optimized table report to data frame generation
* Added new filtration features for genotype tables via filterGenotypeTableSites()
  + parameters for variant sites
  + parameters for physical positions
  + filtration via chromsomome position files
  + filtration via BED file formats


# rTASSEL 0.9.13
* Added error checks for catching C stack usage errors for the following functions:
  + filterGenotypeTableSites()
  + filterGenotypeTableTaxa()
* Added NEWS file for tracking version updates.


# rTASSEL 0.9.12
* Added new functions:
  + `leaveOneFamilyOut()`
  + `genomicPredction()`
* Fixed a bug where tibbles when passed through `readPhenotypeFromDataFrame()`,
  would cause errors.
