# Resolve the transitive runtime dependencies of a Maven artifact

Performs a breadth-first walk of the dependency graph rooted at the
given coordinates, applying Maven's nearest-wins conflict rule: the
version encountered at the shallowest depth is the one that is kept.

## Usage

``` r
resolveDependencies(groupId, artifactId, version, verbose = TRUE)
```

## Arguments

- groupId:

  Maven group identifier.

- artifactId:

  Maven artifact identifier.

- version:

  Artifact version.

- verbose:

  Whether to report progress. Defaults to `TRUE`.

## Value

A data frame with columns `groupId`, `artifactId`, and `version`,
excluding the root artifact itself.
