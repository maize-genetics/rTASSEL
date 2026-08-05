# rTASSEL 0.13.0
* Added installation of TASSEL from the standalone archives published on
  GitHub, which is the only source of nightly builds:
  + `setupTASSEL(source = "github")` installs the newest nightly build cut
    from TASSEL's `develop` branch
  + A specific build can be pinned by version (`"5.2.98-dev.20260801"`) or
    by release tag (`"dev-20260801"`), and tagged releases can be installed
    the same way with `"latest"` or `"v5.2.97"`
  + Archives bundle every dependency, so no dependency resolution is
    needed, and each download is verified against its published SHA-256
    checksum
  + Nightly builds are cached beside released versions under their full
    version, so both can be installed at once and switched between with
    `options(rTASSEL.tassel.version = ...)`
  + `GITHUB_PAT` or `GITHUB_TOKEN` is used for release lookups when set,
    which avoids the anonymous GitHub rate limit
* The startup message now names which channel the loaded JARs came from,
  distinguishing a nightly build from a released version
* Update checks follow the installed channel: a nightly install is checked
  against the newest nightly build rather than Maven Central, so it is
  never prompted to "upgrade" to an older release. The channel can also be
  queried directly with `checkForTASSELUpdate(channel = "nightly")`
* Added automatic TASSEL update checks:
  + `library(rTASSEL)` now checks Maven Central for a newer TASSEL release
    and reports it in the startup message
  + The check runs at most once per day and caches its result on disk
  + Only versions that rTASSEL can actually install are reported, so an
    incomplete upload upstream will not prompt an unusable upgrade
  + Skipped in non-interactive sessions, during `R CMD check`, and on CI;
    network failures are ignored
  + Disable with `options(rTASSEL.check.updates = FALSE)` or the
    `RTASSEL_NO_VERSION_CHECK` environment variable
* Added new function `checkForTASSELUpdate()` to run the check on demand
* Added Maven dependency resolution so `setupTASSEL()` can install TASSEL
  releases published without a bundled dependency JAR:
  + Detects whether a release ships a fat JAR or a thin JAR plus a POM
  + Resolves the transitive dependency tree from the POM, honoring parent
    inheritance, dependency management, BOM imports, scopes, optional
    dependencies, exclusions, and nearest-wins version conflicts
  + Fetches artifacts from Maven Central, SciJava, and JBoss, matching the
    repositories used by TASSEL's own build
  + Drops artifacts superseded by a newer coordinate under a different
    name, such as `google-collections` alongside `guava`, which would
    otherwise shadow each other on the classpath
* The startup message now reports the TASSEL version the JVM actually
  loaded rather than the version rTASSEL expected, which differed when
  JARs came from a user-supplied `rTASSEL.java.path`
* Added tracking of the active TASSEL version:
  + `setupTASSEL()` records the version it installs, so a newly installed
    version is loaded on the next `library(rTASSEL)`
  + Override with `options(rTASSEL.tassel.version = ...)`
  + Fixes a bug where installing any version other than the pinned default
    left rTASSEL reporting "TASSEL JARs not found"
* Every downloaded artifact is now verified against its published SHA-1
  checksum, rather than only the pinned default version
* `setupTASSEL(force = TRUE)` now clears the cached version first, so stale
  dependency JARs cannot linger on the classpath
* Moved `digest` from Suggests to Imports and added `jsonlite` to Imports

# rTASSEL 0.12.0
* Added new `LDResults` class:
  + Stores pairwise LD statistics and analysis parameters from
    `linkageDiseq()`
* Added new `LDRegion` class:
  + Defines genomic regions for highlighting LD blocks on `plotLD()` plots
* Added new function `plotLD()`:
  + Replaces deprecated `ldPlot()` with a redesigned LD heatmap
  + Supports Haploview and viridis-family color schemes
  + LD block highlighting via `LDRegion` objects
  + Optional genomic position track above the heatmap
  + Toggleable site index labels along the diagonal
* Added new function `plotSnpDensity()`:
  + Generates heatmap-style SNP density plots across chromosomes
  + Configurable window size, viridis color palettes, and log-scaled counts
* Removed deprecated `ldPlot()` and `ldJavaApp()` functions
* Removed deprecated interactive Java visualizations section from vignette
* Updated vignette with new LD and SNP density visualization sections

# rTASSEL 0.11.0
* Add fixes to possible issues related to `ggplot2` v4.0
* Bug fixes for `filterGenotypeTableSites()`:
  + Fixed error "missing value where TRUE/FALSE needed" when filtering by 
    position with matching start/end chromosomes and NULL position values
  + Added input validation for negative `startSite` and `endSite` values
  + Added input validation for negative `startPos` and `endPos` values
* Removed bundled JAR files from `inst/java/` to reduce package size
* Added `setupTASSEL()` function to download TASSEL JARs from Maven Central
  and cache them locally
* Added JAR resolution system with priority: user option, Maven cache, then
  bundled (legacy fallback)
* Replaced all `class() ==` checks with `inherits()` for proper S4 class
  handling
* Fixed implicit S4 generic imports from `BiocGenerics` (`colnames`,
  `rownames`, `ncol`, `nrow`)
* Improved package startup messages using `cli` formatting
* Consolidated GitHub Actions workflows into a single `ci.yaml`
* Added unit tests for setup, constants, and package loading
* Added `BiocGenerics` to Imports; `digest` and `withr` to Suggests

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
