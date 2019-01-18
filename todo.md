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
- [ ] Filtration methods
- [ ] implement methods for setting co-variate - Ed
- [ ] if only two columns and one is Taxa - assume other is data - Ed
- [ ] `createPhenoGenoBasedOnFormula()` - implement formula input
- [ ] `fixedEffectLMPlugin()` - implement formula input
- [x] `assocModelDesign()` error - "could not find function data_frame"
- [x] `assocModelDesign()` error - "could not find function add_case"
- [ ] `assocModelDesign()` error - The LHS of `:=` must be a string or a symbol (dplyr)
- [ ] GWAS wrapper
- [ ] LD / Distance matrix class - Ed and Brandon
- [ ] Finalize function generation
- [ ] Code cleanup / purging / function naming - Brandon
- [ ] Finish Roxygen2 documentation for all user functions - Brandon
- [ ] Finish R checks - no warnings or notes - Brandon
- [ ] Finish Bioconductor checks - no warnings or notes - Brandon
- [ ] Finalize vignettes - Group
- [ ] Submit to Bioconductor - Group

## Other things to consider
- [ ] Rocker/Docker integration
