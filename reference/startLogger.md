# Start TASSEL logging information

This function will create a file for storing logging output from TASSEL.

## Usage

``` r
startLogger(filePath = NULL, verbose = TRUE)
```

## Arguments

- filePath:

  File path and name of log file location. If `NULL`, a logging file
  (`rtasel_log.txt` will be added to current working directory.

- verbose:

  Display logging file creation message? Defaults to `TRUE`.
