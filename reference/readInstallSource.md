# Read where a cached TASSEL version was installed from

Installations predating source tracking carry no marker, and only Maven
installs existed then, so an absent marker means Maven.

## Usage

``` r
readInstallSource(version)
```

## Arguments

- version:

  TASSEL version string.

## Value

Either `"maven"` or `"github"`.
