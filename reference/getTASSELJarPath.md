# Get path to cached TASSEL JAR files

Returns the file path to the cached TASSEL JAR directory, or `NULL` if
no cached JARs are found. JARs are cached per-version in the standard R
user cache directory.

Every packaging layout is recognised: a single fat JAR, a thin JAR
accompanied by separately resolved dependency JARs, and the
`sTASSEL.jar` of an unpacked GitHub standalone archive.

## Usage

``` r
getTASSELJarPath(version = getActiveTASSELVersion())
```

## Arguments

- version:

  TASSEL version string. Defaults to the currently active version.

## Value

A character string path to the JAR cache directory, or `NULL` if no
cached JARs exist.
