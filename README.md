## rTASSEL

[![Build Status](https://travis-ci.com/maize-genetics/rTASSEL.svg?branch=master)](https://travis-ci.com/maize-genetics/rTASSEL) [![Life Cycle Status](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/maize-genetics/rTASSEL_sandbox/master?filepath=getting_started.ipynb)

**Website**: https://maize-genetics.github.io/rTASSEL


### Objective
The main goal of developing this package is to construct an R-based front-end
to connect to a variety of highly used TASSEL methods and analytical tools.
By using R as a front-end, we aim to utilize a *unified* scripting workflow that
exploits the analytical prowess of TASSEL in conjunction with R's popular
data handling and parsing capabilities without ever having the user to switch
between these two environments.


### Demo
If you want to test out what this package does but do not want to install it 
locally, we have set up an interactive Jupyter notebook detailing the
walkthrough of `rTASSEL` on Binder. The Binder link can be accessed through
the Binder icon on this page or by clicking 
[here](https://mybinder.org/v2/gh/maize-genetics/rTASSEL_sandbox/master?filepath=getting_started.ipynb).


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


### Citation
If you plan on using this software for your publications, please cite:

**rTASSEL: an R interface to TASSEL for association mapping of complex traits**\
Brandon Monier, Terry M. Casstevens, Peter J. Bradbury, Edward S. Buckler\
bioRxiv 2020.07.21.209114; doi: https://doi.org/10.1101/2020.07.21.209114

or for BibTex users:

```
@article {
	author = {Monier, Brandon and Casstevens, Terry M. and Bradbury, Peter J. and Buckler, Edward S.},
	title = {rTASSEL: an R interface to TASSEL for association mapping of complex traits},
	elocation-id = {2020.07.21.209114},
	year = {2021},
	doi = {10.1101/2020.07.21.209114},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/03/31/2020.07.21.209114},
	eprint = {https://www.biorxiv.org/content/early/2021/03/31/2020.07.21.209114.full.pdf},
	journal = {bioRxiv}
}
```



### Contacts
* Brandon Monier (bm646@cornell.edu) - _Author, Maintainer_
