# Download and configure TASSEL JAR files

Downloads TASSEL and caches it locally. This only needs to be run once
per TASSEL version; subsequent package loads use the cached JARs
automatically.

Releases come from Maven Central by default. The standalone archives on
the TASSEL GitHub releases page can be used instead, which is the only
way to obtain a nightly build.

## Usage

``` r
setupTASSEL(version = NULL, force = FALSE, source = c("maven", "github"))
```

## Arguments

- version:

  Which TASSEL release to install. For `source = "maven"`, a version
  string, defaulting to the version pinned by this release of rTASSEL.
  For `source = "github"`, either `"nightly"` (the newest nightly build,
  the default), `"latest"` (the newest tagged release), a release tag
  such as `"dev-20260801"`, or a nightly version such as
  `"5.2.98-dev.20260801"`.

- force:

  If `TRUE`, discard any cached copy and re-download. Defaults to
  `FALSE`.

- source:

  Where to download from: `"maven"` for Maven Central (the default) or
  `"github"` for a standalone archive from the TASSEL releases page.

## Value

The path to the JAR cache directory (invisibly).

## Details

Maven Central publishes TASSEL in one of two layouts, and the
appropriate one is detected automatically:

- Fat JAR:

  Releases up to 5.2.96 ship a single `jar-with-dependencies` archive
  containing every dependency.

- Thin JAR:

  Later releases ship only TASSEL's own classes, so its dependencies are
  resolved from the POM and downloaded alongside it. This pulls
  artifacts from Maven Central, SciJava, and JBoss, since no single
  repository serves the whole graph.

GitHub standalone archives bundle `sTASSEL.jar` together with every
dependency, so no resolution is needed. Nightly builds are cut from the
`develop` branch and are only published there; they are unstable by
construction and are best pinned to an exact version once a session
depends on one.

Files are cached under the standard R user cache directory (see
[`R_user_dir`](https://rdrr.io/r/tools/userdir.html)) at
`~/.cache/R/rTASSEL/java/<version>/` (Linux),
`~/Library/Caches/org.R-project.R/R/rTASSEL/java/<version>/` (macOS), or
the equivalent on Windows. Nightly builds live beside released versions
under their full version, as in `java/5.2.98-dev.20260801/`, so both can
be installed at once and switched between with
`options(rTASSEL.tassel.version = ...)`.

Every download is verified against its published checksum: SHA-1 for
Maven artifacts and SHA-256 for GitHub assets.

On success the version is recorded as active, so the next
[`library(rTASSEL)`](https://github.com/maize-genetics/rTASSEL) loads
it.

## Examples

``` r
if (FALSE) { # \dontrun{
## Download the default TASSEL version
setupTASSEL()

## Install a specific version
setupTASSEL(version = "5.2.96")

## Force re-download
setupTASSEL(force = TRUE)

## Install the newest nightly build from GitHub
setupTASSEL(source = "github")

## Pin a specific nightly build
setupTASSEL(version = "5.2.98-dev.20260801", source = "github")
} # }
```
