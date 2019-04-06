## rTASSEL

### Objective
The main goal of developing this package is to construct an R-based front-end
to connect to a variety of highly used TASSEL methods and analytical tools.
By using R as a front-end, we aim to utilize a *unified* scripting workflow that
exploits the analytical prowess of TASSEL in conjunction with R's popular
data handling and parsing capabilities without ever having the user to switch
between these two environments.

### Installation
Package source code can be installed directly from this BitBucket repository
using the `remotes` package:

```
if (!require("remotes")) install.packages("remotes")
remotes::install_bitbucket(
    repo = "bucklerlab/rtassel",
    build = TRUE
    build_opts = c("--no-resave-data", "--no-manual")
)
```

### Usage
For an overview of available functions, use the following command:

```
help(package = "rTASSEL")
```

Once the package is installed, long-form documentation and walkthroughs can
be loaded with:

```
utils::vignette("rtassel_walkthrough")
```


If you would like to study a function in full, refer to the R documentation 
by using `?<function>` in the console, where `<function>` is an 
`rTASSEL`-based function.

### Contacts
* Brandon Monier (bm646@cornell.edu)
* Terry Casstevens
* Ed Buckler

- - -
*Last Modified:* 2019-04-05 at 22:45:50
