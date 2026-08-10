# Installing rTASSEL

## Prerequisite - installing rJava

Since TASSEL is written primarily in Java, a Java JDK will need to be
installed on your machine. Additionally, for R to communicate with Java,
the R package [`rJava`](https://www.rforge.net/rJava/) will need to be
installed. In order to use `rTASSEL`, ensure that you have:

- A `JDK` (Java Development Kit $`\geq`$`8`) installed on your system.
- Your system environment variable `JAVA_HOME` is configured
  appropriately and points to your `JDK` of choice. This will usually be
  included in your PATH environment variable as well. Options and system
  environmental variables that are available from R can be seen with
  [`Sys.getenv()`](https://rdrr.io/r/base/Sys.getenv.html) and more
  specifically `Sys.getenv("JAVA_HOME")`.

**NOTE**: If you are using a UNIX system (e.g. Ubuntu) and are
experiencing issues, you may need to reconfigure R with Java. To perform
this, open a terminal and enter the command:

    R CMD javareconf

You may need to have root privileges when performing this so you may
need to add `sudo` to the prior command.

If you need additional steps on how to perform these actions, detailed
information can be found using the following links, depending on your
OS:

- [Linux](https://datawookie.netlify.com/blog/2018/02/installing-rjava-on-ubuntu/)
- [macOS](https://zhiyzuo.github.io/installation-rJava/)
- [Windows](https://cimentadaj.github.io/blog/2018-05-25-installing-rjava-on-windows-10/installing-rjava-on-windows-10/)

## Install from GitHub

### Building with vignettes

After you have `rJava` up and running on your machine, install `rTASSEL`
by installing the source code from our GitHub repository using the
`devtools` package. Here, we show how you can install the package and
build vignettes locally:

``` r

if (!require("devtools")) install.packages("devtools")
devtools::install_github(
    repo = "maize-genetics/rTASSEL",
    ref = "master",
    build_vignettes = TRUE,
    dependencies = TRUE
)
```

The `dependencies = TRUE` parameter will have to be set if you do not
have the suggested packages described in the
[`DESCRIPTION`](https://github.com/maize-genetics/rTASSEL/blob/master/DESCRIPTION)
file of this package.

### Building without vignettes

If you wish to **not** build vignettes, the prior method can be
simplified as shown below:

``` r

if (!require("devtools")) install.packages("devtools")
devtools::install_github("maize-genetics/rTASSEL")
```

## Setting up TASSEL JARs

### Overview

Starting with this version of `rTASSEL`, the TASSEL Java libraries are
no longer bundled with the R package. Instead, they are downloaded from
[Maven
Central](https://mvnrepository.com/artifact/net.maizegenetics/tassel)
and cached locally on your machine. This approach greatly reduces the
package size and makes it easier to update TASSEL independently of the R
package.

Released versions come from Maven Central by default. If you want
changes that have not been released yet, the standalone archives on the
[TASSEL releases
page](https://github.com/maize-genetics/tassel/releases) can be
installed instead - see [Installing nightly
builds](#installing-nightly-builds).

### One-time setup

After installing `rTASSEL`, you must run
[`setupTASSEL()`](https://rtassel.maizegenetics.net/reference/setupTASSEL.md)
**once** to download and cache the required JAR files:

``` r

setupTASSEL()
```

This will:

1.  Download TASSEL (~70 MB) from Maven Central
2.  Verify the file integrity via a SHA-1 checksum
3.  Cache it under the standard R user cache directory:
    - **Linux**: `~/.cache/R/rTASSEL/java/<version>/`
    - **macOS**:
      `~/Library/Caches/org.R-project.R/R/rTASSEL/java/<version>/`
    - **Windows**: `%LOCALAPPDATA%/R/cache/R/rTASSEL/java/<version>/`

Subsequent calls to
[`library(rTASSEL)`](https://github.com/maize-genetics/rTASSEL) will
automatically detect and use the cached JARs - no re-download is needed.

TASSEL is published in two different layouts, and
[`setupTASSEL()`](https://rtassel.maizegenetics.net/reference/setupTASSEL.md)
detects which one applies. Releases up to 5.2.96 ship a single “fat” JAR
that already contains every dependency. Later releases ship only
TASSEL’s own classes, so `rTASSEL` reads the published POM, resolves the
full dependency tree, and downloads those JARs alongside it. Because no
single repository serves the whole tree, dependencies are fetched from
Maven Central, [SciJava](https://maven.scijava.org), and
[JBoss](https://repository.jboss.org) - the same three repositories used
by TASSEL’s own build.

### Installing nightly builds

TASSEL also publishes standalone archives on its [GitHub releases
page](https://github.com/maize-genetics/tassel/releases), including a
**nightly build** cut from the `develop` branch. These are the only way
to get changes that have not yet been released to Maven Central. To
install the newest one:

``` r

setupTASSEL(source = "github")
```

Nightly builds are unstable by construction, so once your work depends
on one it is worth pinning the exact build rather than tracking whatever
is newest:

``` r

setupTASSEL(version = "5.2.98-dev.20260801", source = "github")
```

The same `source` also installs the standalone archive of a tagged
release, either by name or by tag:

``` r

setupTASSEL(version = "latest", source = "github")
setupTASSEL(version = "v5.2.97", source = "github")
```

These archives bundle `sTASSEL.jar` together with every dependency, so
no dependency resolution is needed, and each download is verified
against the SHA-256 checksum GitHub publishes for it. Nightly builds are
cached beside released versions under their full version, as in
`java/5.2.98-dev.20260801/`, so a nightly and a release can be installed
at the same time and switched between with
`options(rTASSEL.tassel.version = ...)`.

If you hit a GitHub rate limit while looking up releases, set
`GITHUB_PAT` (or `GITHUB_TOKEN`) in your environment and the lookup will
use it.

### Checking for new TASSEL versions

When `rTASSEL` is attached in an interactive session, it checks Maven
Central at most once per day for a newer TASSEL release and notes it in
the startup message. A version is only reported once `rTASSEL` has
confirmed it can actually be installed, so a partial or broken upload
upstream will not prompt you to upgrade to something unusable.

If you are running a nightly build, the daily check looks for a newer
nightly build on GitHub instead, so a nightly install is never nagged to
“upgrade” to an older release.

You can run the check yourself at any time:

``` r

checkForTASSELUpdate(force = TRUE)
```

Or query the nightly channel explicitly:

``` r

checkForTASSELUpdate(channel = "nightly")
```

To turn the automatic check off, use either of:

``` r

options(rTASSEL.check.updates = FALSE)
Sys.setenv(RTASSEL_NO_VERSION_CHECK = "true")
```

The check never runs in non-interactive sessions, during `R CMD check`,
or on continuous integration, and a network failure is silently ignored.

### Updating or re-downloading

To move to a newer TASSEL release, pass its version:

``` r

setupTASSEL(version = "5.2.96")
```

The version you install last is recorded as the active one, so the next
[`library(rTASSEL)`](https://github.com/maize-genetics/rTASSEL) picks it
up. To pin a specific cached version instead, set an option before
loading the package:

``` r

options(rTASSEL.tassel.version = "5.2.96")
library(rTASSEL)
```

If you need to re-download the JARs (e.g. a corrupted cache), use the
`force` parameter, which clears the cached copy first:

``` r

setupTASSEL(force = TRUE)
```

### Custom JAR path

Advanced users who maintain their own TASSEL builds can bypass the Maven
cache entirely by setting an R option **before** loading `rTASSEL`:

``` r

options(rTASSEL.java.path = "/path/to/my/tassel/jars")
library(rTASSEL)
```

When this option is set, `rTASSEL` will load JARs from the specified
directory instead of the Maven cache or bundled location.

### JAR resolution order

When `rTASSEL` is loaded, the TASSEL JARs are resolved in the following
priority order:

1.  **User-defined path** via `options(rTASSEL.java.path = ...)`
2.  **Local cache** (from
    [`setupTASSEL()`](https://rtassel.maizegenetics.net/reference/setupTASSEL.md),
    whether from Maven Central or GitHub)
3.  **Bundled `inst/java/`** (legacy fallback for older installations)

If no JARs are found from any source, `rTASSEL` will load without
initializing the JVM and display a message prompting you to run
[`setupTASSEL()`](https://rtassel.maizegenetics.net/reference/setupTASSEL.md).

Which version the cache resolves to is determined, in order, by
`options(rTASSEL.tassel.version = ...)`, the version recorded by your
most recent
[`setupTASSEL()`](https://rtassel.maizegenetics.net/reference/setupTASSEL.md)
call, and finally the version pinned by this release of `rTASSEL`. The
startup message names where the JARs came from, so a nightly install is
easy to tell apart from a released one.

## Loading `rTASSEL`

After installation and the one-time
[`setupTASSEL()`](https://rtassel.maizegenetics.net/reference/setupTASSEL.md)
step, the package can be loaded using:

``` r

library(rTASSEL)
```

    ## ── Welcome to rTASSEL (version 0.13.0) ──
    ## ℹ Running TASSEL version "5.2.96" (maven cache)
    ## ℹ Consider starting a TASSEL log file (see startLogger()
    ##   (`?rTASSEL::startLogger()`))

Or, if you want to use a function without violating your environment you
can use `rTASSEL::<function>`, where `<function>` is an `rTASSEL`
function.

## Running from Docker

If you wish to run a containerized version of `rTASSEL`, build an image
from the
[`docker/Dockerfile`](https://github.com/maize-genetics/rTASSEL/blob/master/docker/Dockerfile)
in the rTASSEL repository. There is no pre-built image published from
this project; you build the image on your machine (or in your own
registry workflow).

**1. Choose an rTASSEL Git tag** (for example `v0.12.0`) for both the
build context and `RTASSEL_TAG`. The URL fragment `#v0.12.0` tells
Docker to fetch that ref from GitHub; keep `--build-arg RTASSEL_TAG`
aligned with the same tag so the container installs the matching
release.

**2. Build the image** from the GitHub URL (no local clone required).
Tag the image so you can use it in `docker run`:

    docker build \
      --build-arg RTASSEL_TAG=v0.12.0 \
      -t rtassel:v0.12.0 \
      -f docker/Dockerfile \
      https://github.com/maize-genetics/rTASSEL.git#v0.12.0

You may use any image name and tag you prefer
(e.g. `-t my-rtassel:local`); the examples below assume
`rtassel:v0.12.0`. If you already have a clone of the repository, you
can run the same command with `.` as the final argument instead of the
GitHub URL.

### With the terminal

After the build finishes, start R inside the container:

    docker run --rm -ti rtassel:v0.12.0 R

### With RStudio Server

This image also contains an RStudio Server instance. To run this, you
will need to publish the container’s port(s) to the host (`-p`). For
example:

    docker run --rm -ti -p 8787:8787 rtassel:v0.12.0

From here, you can go to `localhost:8787` on a web browser and enter a:

- Username (by default, this will be `rstudio`)
- Password (this will be a randomly generated password displayed in the
  terminal output)

## Setting `rTASSEL`/Java memory

### Local overview

Since `rTASSEL` leverages the TASSEL 5 Java API via the `rJava` package,
it is important to allocate sufficient memory to the Java Virtual
Machine (JVM) before it starts. This is done using the
`options(java.parameters = "-Xmx...")` command in R, which sets JVM
parameters such as the maximum heap size (e.g., `-Xmx4g` for **4 GB**).
The reason this must be set *before* loading `rTASSEL` is because the
JVM can only be configured at startup. Once initialized, its memory
settings cannot be changed without restarting the R session. This
becomes especially important when working with large datasets or
computationally intensive method calls, which can quickly exceed the
default memory allocation and lead to `OutOfMemoryError`s. By increasing
the available heap space proactively, we ensure that Java operations can
be performed efficiently and without interruption due to memory
constraints.

In short, if you are loading large genotype datasets and/or phenotype
data, it is adamant that you specify the memory allocated ***before***
loading the `rTASSEL` package via the
[`options()`](https://rdrr.io/r/base/options.html) function:

``` r

# Allocate 4 GB of memory to the JVM
options(java.parameters = "-Xmx4g")

# Load rTASSEL
library(rTASSEL)
```

### Running `rTASSEL` on RStudio Server

Certain instances of RStudio Server on computing clusters can override
what you specify in the prior example (i.e., running the
[`options()`](https://rdrr.io/r/base/options.html) function before
loading `rTASSEL`) due to when the JVM is initialized and the
[`options()`](https://rdrr.io/r/base/options.html) function is called in
the recently initialized R session. If the JVM is initialized, any value
provided to the `java.parameters` key in the
[`options()`](https://rdrr.io/r/base/options.html) call **will be
silently ignored**. To prevent this from happening, make sure to set up
a `.Rprofile` configuration file with the aforementioned
[`options()`](https://rdrr.io/r/base/options.html) call:

    ## Example .Rprofile entry

    # Allocate 4 GB of memory to the JVM
    options(java.parameters = "-Xmx4g")

Since setting up a `.Rprofile` configuration are out of scope for this
package, please refer to [Posit’s excellent write up on the
subject](https://docs.posit.co/ide/user/ide/guide/environments/r/managing-r.html#rprofile).

### Helpful tips

#### Verify memory has been set

If you are running into `OutOfMemory` exceptions, verify if you have
specified enough memory via the prior
[`options()`](https://rdrr.io/r/base/options.html) call. This can help
verify if you have properly set enough memory at startup. By default,
`rJava` will allocate **500 MB (0.5 GB)** of memory to your session. At
any time during your R session you can report the total memory allocated
using a couple of `rJava` calls:

``` r

# Call Java Runtime class
runtime  <- rJava::.jcall("java/lang/Runtime", "Ljava/lang/Runtime;", "getRuntime")

# Get total memory allocation (reported in gigabytes [GB])
gbConv <- 1024^3 # e.g. ~1e9 (billion) bytes
totMem <- round(round(rJava::.jcall(runtime, "J", "totalMemory") / gbConv, 3))

# Show in console
totMem
```

If the `java.parameters` value in
[`options()`](https://rdrr.io/r/base/options.html) was set up properly,
the value specified in `totMem` should be the same value you specify in
[`options()`](https://rdrr.io/r/base/options.html).

#### Ensure enough memory is allocated

In most instances, genotype data will be the main determinant of how
much memory you should allocate to the JVM. In most cases, the amount of
memory you should allocate for genotype data is at least:

    (# taxa) * (# sites) * (1 byte)

For example, if you have genotype data consisting of **250** taxa and
**3000** sites, this would be **750000** bytes or **0.75** megabytes
(MB).

## Prior issues and possible resolutions

### Problems installing rJava on macOS with M1 CPU architecture

If you are running into issues with installing `rJava` using the newer
Mac chip architecture, Oracle JDK currently (as of writing this) does
not work. Consider an alternative JDK source such as
[OpenJDK](https://openjdk.org/) or [Azul
JDK](https://www.azul.com/downloads/?version=java-8-lts&package=jdk).

More detailed information about a possible workaround can be found in
this [Stack Overflow
post](https://stackoverflow.com/questions/67849830/how-to-install-rjava-package-in-mac-with-m1-architecture).

### Problems installing if you have both 32- and 64-bit architecture installed for R

If you are using a machine that has **both architectures installed for
R**, you might run into problems pulling code using `devtools`. If this
is the case, one solution would be to add the parameter `--no-multiarch`
option in `INSTALL_opts`. This will force building the package for your
currently running R version:

``` r

devtools::install_github(
    repo = "maize-genetics/rTASSEL",
    ref = "master",
    build_vignettes = FALSE,
    INSTALL_opts = "--no-multiarch"
)
```

### Problems with `rJava` if you have upgraded Java

On macOS: if you previously had `rJava` working through RStudio, then
you upgraded your Java and it now longer works, try the following:

At the command line type:

    R CMD javareconf

Then check for a left over symbolic link via:

    ls -ltr /usr/local/lib/libjvm.dylib

If the link exists, remove it, then create it fresh via these commands:

    rm /usr/local/lib/libjvm.dylib
    sudo ln -s $(/usr/libexec/java_home)/lib/server/libjvm.dylib /usr/local/lib

You should now be able to enter RStudio and setup `rJava`.
