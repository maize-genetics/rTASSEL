# Check for a newer version of TASSEL

Queries Maven Central for the newest release of TASSEL that rTASSEL is
able to install, and reports it alongside the version currently in use.

With `channel = "nightly"`, the TASSEL GitHub releases page is queried
for the newest nightly build instead.

This check also runs automatically when the package is attached, at most
once per day, against whichever channel the installed version came from.
The result is cached on disk so that repeated calls to
[`library(rTASSEL)`](https://github.com/maize-genetics/rTASSEL) do not
repeatedly hit the network.

## Usage

``` r
checkForTASSELUpdate(force = FALSE, channel = c("maven", "nightly"))
```

## Arguments

- force:

  If `TRUE`, ignore both the cached result and the settings that
  normally suppress automatic checks, and query the release channel
  directly. Defaults to `FALSE`.

- channel:

  Which release channel to query: `"maven"` for Maven Central releases
  (the default) or `"nightly"` for GitHub nightly builds.

## Value

A list with element `latest` (the newest installable version) and
`checkedAt` (the time of the check), or `NULL` if no result is
available. Maven results also carry `layout` (`"fat"` or `"thin"`);
nightly results carry `tag`, `publishedAt`, and `channel`.

## Details

The automatic check is skipped in non-interactive sessions, during
`R CMD check`, and on continuous integration. It can be disabled
outright with either of:


    options(rTASSEL.check.updates = FALSE)
    Sys.setenv(RTASSEL_NO_VERSION_CHECK = "true")

Versions are probed newest-first and only reported when the artifacts
rTASSEL needs are actually present, since a version can appear in Maven
metadata without being installable.

Network failures are never propagated: if the release channel cannot be
reached, the last cached result is returned, or `NULL` if there is none.

## Examples

``` r
if (FALSE) { # \dontrun{
## Check for updates, bypassing the daily cache
checkForTASSELUpdate(force = TRUE)

## Check for a newer nightly build
checkForTASSELUpdate(channel = "nightly")
} # }
```
