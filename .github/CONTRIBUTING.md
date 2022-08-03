# Contributing to rTASSEL

Thank you for taking time to contribute!

Since `rTASSEL` is an open source project and is maintained by a small team, we 
encourage interested members of the research community to seek out issues and 
propose changes where needed. In this document, we detail outlines to how you
may contribute to the `rTASSEL` project.


## Minor changes

You can fix typos, spelling mistakes, or grammatical errors in the 
documentation directly using the GitHub web interface, as long as the changes 
are made in the _source_ file. Since we use 
[`roxygen2`](https://roxygen2.r-lib.org/articles/roxygen2.html) for document
generation, please alter only typos found in `.R` files and _not_ `.Rd` files.
`roxygen2` comments can be identified within `.R` files with the `#'` header.


## Bigger changes and bug reporting

If you want to make a bigger change, it's a good idea to first file an
[issue](https://github.com/maize-genetics/rTASSEL/issues)
and make sure someone from the team agrees that it’s needed. If you’ve found a 
bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write 
a unit test, if needed).

### Submitting a pull request

* Fork the `rTASSEL` package from GitHub and clone onto your computer. If you
  are not familiar with `git` actions, the `usethis` package provides a series
  of `git`-based methods to help you out. For forking, you can use:
  
  ```
  usethis::create_from_github("maize-genetics/rTASSEL", fork = TRUE)
  ```

* From this fork, make sure you have all dependencies installed for `rTASSEL`. 
  This can be accomplished using `devtools::install_dev_deps()`.

* Once all dependencies are installed, run an `R CMD` check using either 
  RStudio's integrated check features under the "Build" tab, or by running 
  `devtools::check()` in the console.

* **If build check does not fail**, create a Git branch for your pull request 
  (PR). The `usethis` package has a method for this that you may wish to use:
  
  ```
  usethis::pr_init("brief-description-of-change")
  ```
* Make your changes, commit to git, and then create a PR by running 
  `usethis::pr_push()` or other `git`-related commands/software, and following 
  the prompts in your browser. **Please provide a descriptive PR**. Your
  PR should have the following:
  + Descriptive title (briefly document the change)
  + What number does this PR fix? Please provide this in the body of the PR
    with the language `Fixes #<issue number>`.

* If accepted, this will be documented in the NEWS material (`NEWS.md`) with
  the next version update and acknowledgement of the user.

### Code style

* We try to follow the guidelines proposed by
  [Bioconductor](https://contributions.bioconductor.org/r-code.html#r-code) in
  terms of code style. When styling code, please only style your own code and
  not restyle code that falls out of your PR scope.

* We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
  [Markdown syntax](https://roxygen2.r-lib.org/articles/rd-formatting.html),
  for documentation.  

### Unit testing
* **Please test your code!** Like all good scientific projects, rigorous testing
  is a must. We use [testthat](https://cran.r-project.org/package=testthat) 
  for unit tests. Contributions with test cases included are easier to accept 
  and will have higher priority for integration.


## Code of Conduct

Please note that the rTASSEL project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
