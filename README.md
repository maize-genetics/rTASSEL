## rtassel

### Objective
The main goal of developing this package is to construct an R-based front-end
to connect to a variety of highly used TASSEL methods and analytical tools.
By using R as a front-end, we aim to utilize a *unified* scripting workflow that
exploits the analytical prowess of TASSEL in conjunction with R's popular
data handling and parsing capabilities without ever having the user to switch
between these two environments.

### Installation
If you have the `devtools` package on your machine, you can use the following
commands:

```
if (!require("devtools")) install.packages("devtools")
devtools::install_bitbucket(
    repo = "bucklerlab/rtassel",
    ref = "master",
    build_vignettes = TRUE
)
```

### Usage
To load example data for this package (after it has been loaded into the
environment), please use the following locations and
commands:

```
# Phenotype path example
phenoPath <- system.file("extdata", "mdp_traits.txt", package = "rTASSEL")
```

If you would like to study a function in full,
refer to the R documentation by using `?<function>` in the console, where 
`<function>` is an `rTASSEL`-based function.

### Contacts
* Brandon Monier (bm646@cornell.edu)
* Ed Buckler

- - -
*Last Modified:* 2019-04-05 at 21:40:14
