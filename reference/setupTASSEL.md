# Download and configure TASSEL JAR files from Maven Central

Downloads the TASSEL fat JAR (with all dependencies) from Maven Central
and caches it locally. This only needs to be run once per TASSEL
version. Subsequent package loads will use the cached JAR automatically.

## Usage

``` r
setupTASSEL(version = TASSEL_MAVEN$VERSION, force = FALSE)
```

## Arguments

- version:

  TASSEL version to download. Defaults to `"5.2.96"`.

- force:

  If `TRUE`, re-download even if a cached version exists. Defaults to
  `FALSE`.

## Value

The path to the JAR cache directory (invisibly).

## Details

The JAR is downloaded from Maven Central at:
<https://mvnrepository.com/artifact/net.maizegenetics/tassel>

Files are cached under the standard R user cache directory (see
[`R_user_dir`](https://rdrr.io/r/tools/userdir.html)) at
`~/.cache/R/rTASSEL/java/<version>/` (Linux),
`~/Library/Caches/org.R-project.R/R/rTASSEL/java/<version>/` (macOS), or
the equivalent on Windows.

A SHA-1 checksum is verified after download to ensure file integrity.

## Examples

``` r
if (FALSE) { # \dontrun{
## Download default TASSEL version
setupTASSEL()

## Force re-download
setupTASSEL(force = TRUE)
} # }
```
