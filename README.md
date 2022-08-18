# rTASSEL <img src="man/figures/logo.png" align="right" width="120" />

[![Build Status](https://travis-ci.com/maize-genetics/rTASSEL.svg?branch=master)](https://travis-ci.com/maize-genetics/rTASSEL) [![Life Cycle Status](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/btmonier/rTASSEL_sandbox/HEAD?labpath=getting_started.ipynb) [![DOI](https://joss.theoj.org/papers/10.21105/joss.04530/status.svg)](https://doi.org/10.21105/joss.04530)


## Objective
The main goal of developing this package is to construct an R-based front-end
to connect to a variety of highly used TASSEL methods and analytical tools.
By using R as a front-end, we aim to utilize a *unified* scripting workflow that
exploits the analytical prowess of TASSEL in conjunction with R's popular
data handling and parsing capabilities without ever having the user to switch
between these two environments. `rTASSEL`  also provide feature and speed
advantages compared to other commonly used R packages. Take a look [here](https://maize-genetics.github.io/rTASSEL/articles/rtassel_benchmarks.html)
for more information.


## Citation
To cite `rTASSEL`, please use the following citation:

> Monier et al., (2022). rTASSEL: An R interface to TASSEL for analyzing genomic diversity. Journal of Open Source Software, 7(76), 4530, https://doi.org/10.21105/joss.04530


## Demo
If you want to test out what this package does but do not want to install it 
locally, we have set up an interactive Jupyter notebook detailing the
walkthrough of `rTASSEL` on Binder. The Binder link can be accessed through
the Binder icon on this page or by clicking 
[here](https://mybinder.org/v2/gh/btmonier/rTASSEL_sandbox/HEAD?labpath=getting_started.ipynb).


## Installation and usage
If you do not have experience working with and setting up `rJava` with your
R installation, *it is recommended that you read the long-form documentation*.
This walkthrough can be found [here](https://maize-genetics.github.io/rTASSEL/articles/rtassel_walkthrough.html).
If you are already fairly comfortable working with Java JDK and `rJava`, you
can follow the following commands.

Package source code can be installed directly from this BitBucket repository
using the `devtools` package:

```
if (!require("devtools")) install.packages("devtools")
devtools::install_github("maize-genetics/rTASSEL")
```

Vignettes (`build_vignettes`) are optional since there are constantly updated 
article links on our [website](https://maize-genetics.github.io/rTASSEL/). If 
you do want to build vignettes locally, please use the following instructions:

```
if (!require("devtools")) install.packages("devtools")
devtools::install_github(
    repo = "maize-genetics/rTASSEL",
    ref = "master",
    build_vignettes = TRUE,
    dependencies = TRUE
)
```


## Getting help
For an overview of available functions, use the following command:

```
help(package = "rTASSEL")
```

If you need a walkthrough for potential pipelines, long-form documentation can 
be found on our [website](https://maize-genetics.github.io/rTASSEL/) including
a [getting started article](https://maize-genetics.github.io/rTASSEL/articles/rtassel_walkthrough.html).

If you would like to study a function in full, refer to the R documentation
by using `?<function>` in the console, where `<function>` is an
`rTASSEL`-based function.


