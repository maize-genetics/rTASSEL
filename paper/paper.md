---
title: 'rTASSEL: An R interface to TASSEL for analyzing genomic diversity'
tags:
  - R
  - GWAS
  - genetics
  - genomic diversity
  - genomic prediction
authors:
  - name: Brandon Monier
    orcid: 0000-0001-6797-1221
    affiliation: 1
  - name: Terry M. Casstevens
    orcid: 0000-0001-7602-0487
    affiliation: 1
  - name: Peter J. Bradbury
    orcid: 0000-0003-3825-8480
    affiliation: "1, 2"
  - name: Edward S. Buckler
    orcid: 0000-0002-3100-371X
    affiliation: "1, 2"
affiliations:
 - name: Institute for Genomic Diversity, Cornell University, Ithaca, NY 14853
   index: 1
 - name: United States Department of Agriculture-Agricultural Research Service, Robert W. Holley Center for Agriculture and Health, Ithaca, NY 14853
   index: 2
date: 2 June 2022
bibliography: paper.bib
---

# Summary

The need for efficient tools and applications for analyzing genomic diversity is 
essential for any genetics research or breeding program. One commonly used tool, 
TASSEL (**T**rait **A**nalysis by a**SS**ociation, **E**volution, and 
**L**inkage), provides many core methods for genomic analyses. Despite its 
efficiency, TASSEL has limited means to use scripting languages for reproducible 
research and to interact with other analytical tools. Here we present an R 
package `rTASSEL`, a front-end to connect to a variety of highly used TASSEL 
methods and analytical tools. The goal of this package is to create a unified 
scripting workflow that leverages the analytical prowess of TASSEL, in 
conjunction with R's data handling and visualization capabilities without ever 
having the user switch between these two environments.



# Statement of need

As breakthroughs in genotyping technologies allow for evermore available variant 
resources, methods and implementations to analyze complex traits in a diverse 
array of organisms are needed. One such resource is TASSEL (**T**rait 
**A**nalysis by a**SS**ociation, **E**volution, and **L**inkage). This software 
suite contains functionality for analyses in association studies, linkage 
disequilibrium (LD), kinship, and dimensionality reduction (e.g., PCA and 
MDS) [@bradbury_tassel_2007]. While initially released in 2001, the fifth version, 
TASSEL 5, has been optimized for handling large data sets and has added newer 
approaches to association analyses for many thousands of traits [@shabalin_matrix_2012]. 
Despite these improvements, interacting with TASSEL has been limited to either a 
graphical user interface with limited workflow reproducibility or a command-line 
interface with a higher learning curve that can dissuade novice researchers and
provide unnecessary intermediate files in an analytics workflow 
[@zhang_software_2009]. To remediate this issue, we have created an R package, 
`rTASSEL`. This package interfaces the analytical power of TASSEL with R’s data 
formats and intuitive function handling.



# Approach

## Implementation

![_Overview of the `rTASSEL` workflow_. Genotypic and phenotypic data (A) are used to create an R S4 object (B). From this object, TASSEL functionalities can be called to run various association, linkage disequilibrium, and relatedness functions (C). Outputs from these TASSEL analyses are returned to the R environment as data frame objects (D), Manhattan plot visualizations (E), or interactive visualizations for linkage disequilibrium analysis (F).\label{fig:rtassel_concepts}](img/rtassel_concepts.png)

`rTASSEL` combines TASSEL’s abilities to store genotype data as half bytes, 
bitwise arithmetic for kinship analyses, genotype filtration, extensive forms 
of linear modeling, multithreading, and access to a range of native libraries 
while providing access to R’s prominent scripting capabilities and commonly 
used Bioconductor classes [@gentleman_bioconductor_2004; @lawrence_software_2013; @morgan_summarizedexperiment_2020]. 
Since TASSEL is written in Java, a Java to R interface is implemented via the 
`rJava` package [@urbanek_rjava_2013].

`rTASSEL` allows for the rapid import, analysis, visualization, and export of 
various genomic data structures. Diverse formats of genotypic 
information can be used as inputs for `rTASSEL`. These include variant call 
format (`.vcf`), HapMap (`.hmp.txt`), and Flapjack (`.flpjk.*`). Phenotype data can 
also be supplied in multiple formats. These include TASSEL formatted data sets 
or R data frame objects (\autoref{fig:rtassel_concepts}A).

Once data is imported, an S4 object is constructed that is used for all 
downstream analyses (\autoref{fig:rtassel_concepts}B, \autoref{fig:rtassel_concepts}C). 
To construct this object, the function, `readGenotypePhenotype`, is used. This 
S4 object contains slots that exclusively hold references to objects held in 
the Java virtual machine, which can be called with downstream functions. Prior 
to analysis, genotype objects can be quickly imported and filtered in several 
ways to help in the reduction of confounding errors. `rTASSEL` can filter 
genotype objects by either variant site properties (`filterGenotypeTableSites`) 
or by individuals (`filterGenotypeTableTaxa`).


## Association functions

One of TASSEL’s most dynamic functionalities is its capability to perform 
various association modeling techniques. `rTASSEL` allows several types of 
association studies to be conducted using one primary function, 
`assocModelFitter`, with different parameter inputs. This allows for implementing 
both least-squares fixed-effects general linear models (GLM) and mixed linear 
models (MLM) via the $Q + K$ method [@yu_unified_2006]. If no genotypic data is 
provided to the GLM model, `assocModelFitter` can calculate best linear unbiased 
estimates (BLUEs). Additionally, fast GLM approaches are implemented in 
`rTASSEL`, which allow for the rapid analysis of many phenotypic traits 
[@shabalin_matrix_2012].

The data model for analysis can be specified by a formula like R’s `lm` function, 
which is shown as follows:

$$y \sim A_{1} + A_{2} + \cdots + A_{n}$$

Where $y$ is phenotype data, and $A_{n}$ is any covariate or factor data. This 
formula parameter and several other parameters allow the user to run BLUE, GLM, 
or MLM modeling. Once association analysis is completed, TASSEL table reports 
of association statistics are generated as an R list which can then be exported 
as flat files or converted to data frames (\autoref{fig:rtassel_concepts}D). 
`rTASSEL` can also visualize association statistics with the function, 
`manhattanPlot`, which utilizes the graphical capabilities of the package, 
`ggplot2` (Wickham 2016) [@wickham_ggplot2_2016] (\autoref{fig:rtassel_concepts}E).


## Linkage disequilibrium

rTASSEL can also generate linkage disequilibrium (LD) from genotype data via 
the function `linkageDiseq`. LD is estimated by the standardized disequilibrium 
coefficient, $D'$, correlation between alleles at two loci ($r^2$), and 
subsequent $p$-values via a two-sided Fisher’s Exact test. TASSEL Table reports 
for all pairwise comparisons are generated as `data.frame` objects and heatmap 
visualizations for each given metric are generated via TASSEL’s legacy LD 
Java viewer or `ggplot2` (\autoref{fig:rtassel_concepts}F).


## Relatedness functions

For users to run MLM methods, relatedness estimates need to be calculated. `rTASSEL` can 
efficiently compute this on large data sets by processing blocks of sites at a time using bitwise 
operations. This can be accomplished using the function `kinshipMatrix`, which will generate 
a kinship matrix from genotype data. Several methods for calculating kinship in TASSEL are 
implemented. By default, a “centered” identity by state (IBS) approach is used [@endelman_shrinkage_2012]. 
Additionally, normalized IBS [@yang_gcta_2011], dominance-centered IBS [@munoz_unraveling_2014], 
and dominance normalized IBS [@zhu_dominance_2015] can be used. `rTASSEL` can either generate a 
reference object for association analysis or an R `matrix` object via R's `as.matrix` function 
for additional analyses. In addition to kinship generation, principal components
analysis and multidimensional scaling can be used on genotype data using 
`rTASSEL` methods, `pca` and `mds`, respectively. Finally, phylogenetic analysis
can be performed on genotype data using the `createTree` method which will
generate `phylo` objects commonly used by the `ape` package [@paradis_and_schliep_2019].


## Genomic prediction

The function `genomicPrediction` can be used for predicting phenotypes from 
genotypes. To do this, `genomicPrediction` uses genomic best linear unbiased 
predictors (gBLUPs). It proceeds by fitting a mixed model that uses kinship to 
capture covariance between taxa. The mixed model can calculate BLUPs for taxa 
that do not have phenotypes based on the phenotypes of lines with relationship 
information. 



# Additional resources

More information about various functionalities and workflows can be found on our 
[project webpage](https://maize-genetics.github.io/rTASSEL/). Source code can be 
found on our [GitHub repository](https://github.com/maize-genetics/rTASSEL). An 
interactive Jupyter notebook session detailing additional `rTASSEL` workflows 
can be found on [Binder](https://mybinder.org/v2/gh/btmonier/rTASSEL_sandbox/HEAD?labpath=getting_started.ipynb).



# Acknowledgements

This project is supported by the USDA-ARS, the Bill and Melinda Gates 
Foundation, and NSF IOS #1822330. We thank Sara J. Miller, Guillaume Ramstein, 
and Joseph Gage for their insightful suggestions on this manuscript and pipeline 
testing. 



# References
