# Get path to cached TASSEL JAR files

Returns the file path to the cached TASSEL JAR directory, or `NULL` if
no cached JARs are found. JARs are cached per-version in the standard R
user cache directory.

## Usage

``` r
getTASSELJarPath(version = TASSEL_MAVEN$VERSION)
```

## Arguments

- version:

  TASSEL version string. Defaults to the version bundled with this
  release of rTASSEL.

## Value

A character string path to the JAR cache directory, or `NULL` if no
cached JARs exist.
