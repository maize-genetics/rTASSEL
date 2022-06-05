## rTASSEL

[![Build Status](https://travis-ci.com/maize-genetics/rTASSEL.svg?branch=master)](https://travis-ci.com/maize-genetics/rTASSEL) [![Life Cycle Status](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/btmonier/rTASSEL_sandbox/HEAD?labpath=getting_started.ipynb)

**Website**: https://maize-genetics.github.io/rTASSEL


### Objective
The main goal of developing this package is to construct an R-based front-end
to connect to a variety of highly used TASSEL methods and analytical tools.
By using R as a front-end, we aim to utilize a *unified* scripting workflow that
exploits the analytical prowess of TASSEL in conjunction with R's popular
data handling and parsing capabilities without ever having the user to switch
between these two environments. `rTASSEL`  also provide feature and speed
advantages compared to other commonly used R packages. Take a look [here](https://maize-genetics.github.io/rTASSEL/articles/rtassel_benchmarks.html)
for more information.


### Demo
If you want to test out what this package does but do not want to install it 
locally, we have set up an interactive Jupyter notebook detailing the
walkthrough of `rTASSEL` on Binder. The Binder link can be accessed through
the Binder icon on this page or by clicking 
[here](https://mybinder.org/v2/gh/btmonier/rTASSEL_sandbox/HEAD?labpath=getting_started.ipynb).


### Installation
If you do not have experience working with and setting up `rJava` with your
R installation, *it is recommended that you read the long-form documentation*.
This walkthrough can be found [here](https://maize-genetics.github.io/rTASSEL/articles/rtassel_walkthrough.html).
If you are already fairly comfortable working with Java JDK and `rJava`, you
can follow the following commands.

Package source code can be installed directly from this BitBucket repository
using the `devtools` package:

```
if (!require("devtools")) install.packages("devtools")
devtools::install_bitbucket(
    repo = "bucklerlab/rTASSEL",
    ref = "master",
    build_vignettes = FALSE
)
```

Vignettes (`build_vignettes`) are optional since there are constantly updated 
article links on our [website](https://maize-genetics.github.io/rTASSEL/). If 
you do want to build vignettes locally, make sure you have the required 
packages and programs available.


### Usage
For an overview of available functions, use the following command:

```
help(package = "rTASSEL")
```

If you need a walkthrough for potential pipelines, long-form documentation can 
be found on our [website](https://maize-genetics.github.io/rTASSEL/).
If you prefer to compile a vignette locally, you can set the `build_vignettes`
parameter to `TRUE` when you download from Bitbucket (*Note:* compiling the
vignette may take some time to process and analyze test data.)


If you would like to study a function in full, refer to the R documentation
by using `?<function>` in the console, where `<function>` is an
`rTASSEL`-based function.


### Contacts
* Brandon Monier (bm646@cornell.edu) - _Author, Maintainer_
