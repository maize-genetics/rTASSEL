# Contributing to rTASSEL

Thank you for taking time to contribute!

Since `rTASSEL` is an open source project and is maintained by a small team, we 
encourage interested members of the research community to seek out issues and 
propose changes where needed. In this document, we detail outlines to how you
may contribute to the `rTASSEL` project. Because our build and check process
is virtually similar to `usethis` and other `tidyverse` packages, 


## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the 
documentation directly using the GitHub web interface, as long as the changes 
are made in the _source_ file. Since we use 
[`roxygen2`](https://roxygen2.r-lib.org/articles/roxygen2.html) for document
generation, please alter only typos found in `.R` files and _not_ `.Rd` files.
`roxygen2` comments can be identified within `.R` files with the `#'` header.


## Bigger changes

If you want to make a bigger change, it's a good idea to first file an
[issue](https://github.com/maize-genetics/rTASSEL/issues)
and make sure someone from the team agrees that it’s needed. If you’ve found a 
bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write 
a unit test, if needed).

### Pull request process

* Fork the package and clone onto your computer. If you haven't done this 
  before, we recommend using 
  `usethis::create_from_github("maize-genetics/rTASSEL", fork = TRUE)`.

* Install all development dependencies with `devtools::install_dev_deps()`, and 
  then make sure the package passes R CMD check by running `devtools::check()`. 
  If R CMD check doesn't pass cleanly, it's a good idea to ask for help before 
  continuing. 

* Create a Git branch for your pull request (PR). We recommend using
  `usethis::pr_init("brief-description-of-change")`.

* Make your changes, commit to git, and then create a PR by running 
  `usethis::pr_push()`, and following the prompts in your browser. The title of
  your PR should briefly describe the change. The body of your PR should 
  contain `Fixes #issue-number`.

* For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just 
  below the first header). Follow the style described in 
  <https://style.tidyverse.org/news.html>.

### Code style

* We try to follow the guidelines proposed by
  [Bioconductor](https://contributions.bioconductor.org/r-code.html#r-code) in
  terms of code style. When styling code, please only style your own code and
  not restyle code that falls out of your PR scope.

* We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
  [Markdown syntax](https://roxygen2.r-lib.org/articles/rd-formatting.html),
  for documentation.  

* **Please test your code!** Like all good scientific projects, rigorous testing
  is a must. We use [testthat](https://cran.r-project.org/package=testthat) 
  for unit tests. Contributions with test cases included are easier to accept 
  and will have higher priority for integration.


## Code of Conduct

Please note that the rTASSEL project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
