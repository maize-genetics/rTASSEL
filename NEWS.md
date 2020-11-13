# CHANGES IN VERSION 0.9.21
* Added new function:
  + `exportGenotypeTable()`
* Added new vignette:
  + "Filtering Genotype Tables"


# rTASSEL 0.9.20
* Added new parameter to `filterGenotypeTableSites()`
  + `gRangesObj`: Filter genotype tables by using a `GRanges` object.
* Added new parameter to `filterGenotypeTableTaxa()`
  + `taxa`: Pass a vector of taxa IDs to filter genotype table by.
* Fixed `getSumExpFromGenotypeTable()` bug:
  + dosage array now returns `NA`s instead of `128` values.


# rTASSEL 0.9.19
* Added two new parameters to `filterGenotypeTableSites()`
  + `removeMinorSNPStates`: Boolean; removes minor SNP states.
  + `removeSitesWithIndels`: Boolean; removes sites with indels.
* Added better descriptive error handling for `filterGenotypeTableSites()`
* Fixed `siteRangeFilterType` parameter bug in `filterGenotypeTableSites()`. 
  Now defaults to `none` when user does not specify filter type.
* Added two new parameters to `getSumExpFromGenotypeTable()`
  + `coerceDosageToInt`: Returns raw byte dosage array instead of integer from Java.
  + `verbose`: Display console messages for large "memory-intensive" datasets.


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
