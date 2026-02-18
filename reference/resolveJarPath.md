# Resolve the JAR path from all available sources

Resolves the JAR path in priority order:

1.  User-defined path via `options(rTASSEL.java.path = ...)`

2.  Maven cache (from
    [`setupTASSEL`](https://rtassel.maizegenetics.net/reference/setupTASSEL.md))

3.  Bundled `inst/java/` (legacy fallback)

## Usage

``` r
resolveJarPath(pkgname = "rTASSEL", libname = NULL)
```

## Arguments

- pkgname:

  Package name (used for bundled path lookup).

- libname:

  Library location (used for bundled path lookup).

## Value

A list with elements `path` (character or `NULL`) and `source`
(`"option"`, `"maven cache"`, `"bundled"`, or `NULL`).
