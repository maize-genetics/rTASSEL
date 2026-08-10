# Get the TASSEL version rTASSEL should load

Resolves the active version in priority order:

1.  `options(rTASSEL.tassel.version = ...)`

2.  The version recorded by the most recent
    [`setupTASSEL`](https://rtassel.maizegenetics.net/reference/setupTASSEL.md)

3.  The version pinned by this release of rTASSEL

## Usage

``` r
getActiveTASSELVersion()
```

## Value

A character string version.
