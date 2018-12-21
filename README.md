## rtassel

### Objective
To creat an R front-end for the Java program, `TASSEL` using `rJava`

**NOTE**: Main branch functionalities are currently broken as we implement new
methods.

### Installation
If you have the `devtools` package on your machine, you can use the following
commands:

```
if (!require("devtools")) install.packages("devtools")
devtools::install_bitbucket("bucklerlab/rtassel")
```

### Usage
To load example data for this package (after it has been loaded into the
environment), please use the following locations and
commands:

```
# Phenotype path example
phenoPath <- system.file("extdata", "mdp_traits.txt", package = "rTASSEL")
```

- - -
*Last Modified:* 2018-12-21 at 15:26:30
