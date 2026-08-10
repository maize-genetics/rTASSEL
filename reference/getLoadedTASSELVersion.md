# Get the TASSEL version actually loaded in the JVM

Asks TASSEL which version it is, which is the only reliable answer when
JARs come from a user-supplied path that need not match any version
rTASSEL has recorded. Falls back to the configured version if TASSEL
cannot be queried.

## Usage

``` r
getLoadedTASSELVersion(fallback = getActiveTASSELVersion())
```

## Arguments

- fallback:

  Version to return if the JVM cannot be queried.

## Value

A character string version.
