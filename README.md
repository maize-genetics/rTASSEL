# rTASSEL <img src="man/figures/logo.png" align="right" width="120" />

[![Life Cycle Status](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![R-CMD-check](https://github.com/maize-genetics/rTASSEL/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/maize-genetics/rTASSEL/actions/workflows/check-standard.yaml) [![codecov](https://codecov.io/gh/maize-genetics/rTASSEL/branch/master/graph/badge.svg?token=4D0JSKT0UC)](https://codecov.io/gh/maize-genetics/rTASSEL) [![DOI](https://joss.theoj.org/papers/10.21105/joss.04530/status.svg)](https://doi.org/10.21105/joss.04530) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/btmonier/rTASSEL_sandbox/HEAD?labpath=getting_started.ipynb) 

`rTASSEL` is an R-based front-end for accessing key TASSEL 5 methods and tools.
This allows users to run powerful TASSEL 5 analyses within a unified R workflow,
without switching between environments. `rTASSEL` also offers performance and
feature advantages over other R packages. See [this benchmark](https://maize-genetics.github.io/rTASSEL/articles/rtassel_benchmarks.html)
for details.


## Citation
To cite `rTASSEL`, please use the following citation:

> Monier et al., (2022). rTASSEL: An R interface to TASSEL for analyzing 
> genomic diversity. Journal of Open Source Software, 7(76), 4530, 
> https://doi.org/10.21105/joss.04530


## Installation

### Stable release
```{r}
# install.packages("pak")
pak::pak("maize-genetics/rTASSEL@v0.10.0")
```

### Developmental version
```{r}
# install.packages("pak")
pak::pak("maize-genetics/rTASSEL")
```

> [!NOTE]  
> Since `rTASSEL` uses `rJava`, you will need a working version of Java 
> ($\geq$ 8). 
>
> **Mac/Linux users**: you may need to run `R CMD javareconf` if
> you run into issues with installing `rJava` via `pak`. More installation tips
> can be found [here](https://rtassel.maizegenetics.net/articles/rtassel_installation.html).


## Demo
If you want to test out what this package does but do not want to install it 
locally, we have set up an interactive Jupyter notebook detailing the
walkthrough of `rTASSEL` on Binder. The Binder link can be accessed through
the Binder icon on this page or by clicking 
[here](https://mybinder.org/v2/gh/btmonier/rTASSEL_sandbox/HEAD?labpath=getting_started.ipynb).


## Getting help
For an overview of available functions, use the following command:

```
help(package = "rTASSEL")
```

If you need a walkthrough for potential pipelines, long-form documentation can 
be found on our [website](https://rtassel.maizegenetics.net) including
a [getting started article](https://rtassel.maizegenetics.net/articles/rTASSEL.html).

If you would like to study a function in full, refer to the R documentation
by using `?<function>` in the console, where `<function>` is an
`rTASSEL`-based function.


