# rTASSEL 0.11
* Add fixes to possible issues related to `ggplot2` v4.0
* Bug fixes for `filterGenotypeTableSites()`:
  + Fixed error "missing value where TRUE/FALSE needed" when filtering by 
    position with matching start/end chromosomes and NULL position values
  + Added input validation for negative `startSite` and `endSite` values
  + Added input validation for negative `startPos` and `endPos` values

# rTASSEL 0.10.0
* Updated formula parsing:
  + Drop and keep traits by using `-` and `+`, respectively
  + Added keywords to drop/keep all covariate or factor traits:
    + `I(cov)`
    + `I(fct)`
* Added `stepwiseModelFitter()` function:
  + Runs stepwise regression via TASSEL 5's "Stepwise" plugin
  + Returns `AssociationResults` object
* Added new function `readNumericGenotypeFromRMatrix()`:
  + Converts formatted R matrix to TASSEL 5 numeric genotype
* Added experimental function `readGenotype()`:
  + Reads data into new `TasselGenotype` class (future update)
  + Dynamically read genotype data based data type
* Added experimental function `readPhenotype()`:
  + Reads data into new `TasselPhenotype` class (future update)
  + Dynamically read phenotype data based data type
* Added deprecation warnings to the following methods:
  + `getPhenotypeDF()`
  + `getSumExpFromGenotypeTable()`
  + `readGenotypeTableFromPath()`
  + `readGenotypeTableFromGigwa()`
  + `readGenotypePhenotype()`
  + `readPhenotypeFromDataFrame()`
  + `readPhenotypeFromPath()`
  + `ldJavaApp()`
  + `treeJavaApp()`
  + **NOTE**: these will be removed in future updates
* Updated vignettes:
  + Visualization section updates
  + Formula parsing
  + `AssociationResults` and `PCAResults` object handling
  + Memory allocation guide
  + Numeric genotype handling


# rTASSEL 0.9.33
* Fixed typo in `plotPCA()` error message
* Add new function `filterGenotypeTableBySiteName()`:
  + Filters genotype tables using literal marker names/IDs
* Add new function `mergeGenotypeTables()`:
  + Merges multiple genotype tables by site values


# rTASSEL 0.9.32
* Updated `tableReport()` method dispatch for all `AssociationResults`
  objects:
  + Will now return default statistics output for all association results
    when running `tableReport(assocObj)` where `assocObj` is an object
    of type `AssociationResults`
* Removed HDF5 file export support
* Improved logic support for `plotPCA()`


# rTASSEL 0.9.31
* Added new `PCAResults` class
  + Allows for more controlled access of data and simplified downstream
    functions for end users
* Add new function `plotScree()`:
  + Generates quick scree plots from `PCAResults` objects
* Add new function `plotPCA()`:
  + Generates quick PCA plots from `PCAResults` objects
  + Allows for grouping from generated hierarchical clustering or grouping
    from metadata via the `metadata` parameter and subsequent `mCol`
    parameters.


# rTASSEL 0.9.30
* Added new `AssociationResults` class
  + Allows for more controlled access of data and simplified downstream
    functions for end users
* Added new function `plotManhattan()`:
  + Supercedes older Manhattan plotting methods to work with new
    `AssociationResults` class.
* Added new function `plotQQ()`:
  + Plotting function for QQ results from `AssociationResults` class
* Added new function `plotManhattanQC()`:
  + Plotting function and QC method for zoomed in regions of interest
    across genome
* Prior 3 functions also include interactive component that wraps
  `ggplot2` objects with PlotlyJS components


# rTASSEL 0.9.29
* Added genotype table summary methods:
  + `positionList()`
  + `taxaSummary()`
  + `siteSummary()`
* `TasselGenotypePhenotype` objects containing genotype table data can now
  be coerced into R `matrix` objects using the function `as.matrix()`
  + This will return a taxa x site matrix where taxa is the number of rows and
    sites is the number of columns.
* Added generalized join methods:
  + `intersectJoin()`
  + `unionJoin()`
  + Joins phenotype data data based on taxa ID - similar to the TASSEL API
* Added read method for importing GIGWA data through `QBMS`:
  + `readGenotypeTableFromGigwa()`


# rTASSEL 0.9.28
* Fixed `log4j` warning issue
  + This removes `log4j` warning messages when the `startLogger()` function
    is called.
* Removed `useRef` parameter from `getSumExpFromGenotypeTable()` function.
  + This is now automatically detected from the file input.
  + This fixes ref/alt allele vs major/minor allele encoding issues.
* Added Journal of Open Source Software citation for the `rTASSEL` package.
  + For citation information, use `utils::citation("rTASSEL")`
* Added data object, `rtPaths`
  + Includes paths to external toy data for `rTASSEL`


# rTASSEL 0.9.27
* No significant updates in this version. This version is virtually identical
  to `0.9.26` and is for linking to Zenodo for archival purposes.


# rTASSEL 0.9.26
* Bug fixes:
  + Fixed `r2` parameter bug in `ldPlot()`
  + Fixed space bugs in certain column names of data frame objects. 
    `_` values now replace spaces.
  + Fixed `show()` method for `TasselDistanceMatrix` objects.
* Add new function:
  + `seqDiversity()`
  + Calculates diversity basic diversity metrics on genetic data.


# rTASSEL 0.9.25
* Bug fixes:
  + Fixed character conversion bug in `DataFrame` object returns.
* `pca()` can optionally report eigenvalues and eigenvectors as a list object.
* Added new function:
  + `imputeNumeric()`
  + Allows for numeric imputation of `GenotypeTable` objects.
* Added new function:
  + `imputLDKNNi()`
  + Allows for LD KNNi imputation of `GenotypeTable` objects.


# rTASSEL 0.9.24
* Added new function:
  + `pca()`
  + Allows for user to run PCA on rTASSEL objects containing a `GenotypeTable`
    object.
* Added new function:
  + `mds()`
  + Allows for user to run MDS on `TasselDistanceMatrix` objects.
* Enhancements:
  + New summary print output for `TasselDistanceMatrix` objects.


# rTASSEL 0.9.23
* Added new `TasselDistanceMatrix` class
  + Specified function (`kinshipMatrix()` and `distanceMatrix()`) now return
    an object of type `TasselDistanceMatrix`.
  + Prevents console overload and freezing as seen with large distance matrix
    objects.
  + Now shows summary overview of matrix instead of Java object reference.
  + Generic functions `colnames()`, `rownames()`, `ncol()`, and `nrow()` will
    return relative information similar to how these operate with `matrix`
    type objects.
  + Primitive function `as.matrix()` now supersedes deprecated functions
    `kinshipToRMatrix()` and `distanceToRMatrix()`.
  + Prior functions that take in a kinship object will now take in this new
    class.
* Added new function:
  + `readTasselDistanceMatrix()`
  + Allows for user to read in delimited distance matrix stored in a flat
    file.
* Added new function:
  + `asTasselDistanceMatrix()`
  + Coerces a pairwise matrix (e.g. m x m dimensions) with the same column
    and row names of type `matrix` to an object of type `TasselDistanceMatrix`.
* Added new function:
  + `createTree()`
  + interface to TASSEL's tree creation methods
  + Allows for `Neighbor_Joining` and `UPGMA` methods
* Added new function:
  + `treeJavaApp()`
  + wrapper for TASSEL's interface to the Archaeopteryx Java tree Viewer
  + Implements same methods for tree creation as `createTree()`


# rTASSEL 0.9.22
* Fix `manhattanPlot()` aesthetics:
  + Remove redundant marker labels from x-axis
  + Change x-axis label to `SNP Positions`


# rTASSEL 0.9.21
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
