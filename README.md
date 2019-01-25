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

### Contacts
* Brandon Monier (bm646@cornell.edu)
* Ed Buckler

- - -
*Last Modified:* 2019-01-25 at 13:48:14
