# rTASSEL TODO list

## Overview
* We currently have in total, 46 functions
* Current package size is ~40 MB

## TASSEL Side
- [ ] Console log output - silence output on TASSEL side - Terry
- [ ] Formula generation for GWAS - Ed 
- [ ] Dealing with random or fixed effects - Peter
- [ ] Reduce size of jar files? - current size of package is around 40 MB
      can we get this down to around 5 MB?

## R Side
- [ ] Implement `createPhenoGenoBasedOnFormula()` function
- [x] `assocModelDesign()` error - "could not find function data_frame"
- [x] `assocModelDesign()` error - "could not find function add_case"
- [ ] `assocModelDesign()` error - `:=` must be a string or a symbol - dplyr
- [ ] `fixedEffectLMPlugin()` - implement formula input
- [ ] GWAS wrapper
- [ ] TasselGWAS class? - implement methods for visualization, extraction, etc.
- [ ] Finalize function generation
- [ ] Code cleanup / purging / function naming
- [ ] Finish Roxygen2 documentation for all user functions
- [ ] Finish R checks - no warnings or notes
- [ ] Finish Bioconductor checks - no warnings or notes
- [ ] Finalize vignettes

## Other things to consider
- [ ] Rocker/Docker integration
