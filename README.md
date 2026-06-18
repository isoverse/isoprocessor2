
<!-- README.md is generated from README.Rmd. Please edit that file -->

# isoprocessor2

<!-- badges: start -->

[![Documentation](https://img.shields.io/badge/docs-online-green.svg)](https://isoprocessor2.isoverse.org/)
[![R-CMD-check](https://github.com/isoverse/isoprocessor2/workflows/R-CMD-check/badge.svg)](https://github.com/isoverse/isoprocessor2/actions)
[![Codecov test
coverage](https://codecov.io/gh/isoverse/isoprocessor2/graph/badge.svg)](https://app.codecov.io/gh/isoverse/isoprocessor2)
<!-- badges: end -->

## Overview

[isoprocessor2](https://isoprocessor2.isoverse.org/) provides data
processing and reduction pipelines for stable isotope data read with
[isoreader2](https://isoreader2.isoverse.org/). It succeeds the
[isoprocessor](https://isoprocessor.isoverse.org/) package, rebuilt to
work directly with the aggregated data from
[isoreader2](https://isoreader2.isoverse.org/).

This package is in early development.

## Installation

[isoprocessor2](https://isoprocessor2.isoverse.org/) is not yet on the
Comprehensive R Archive Network (CRAN) but you can install the latest
version from [GitHub](https://github.com/isoverse/isoprocessor2) as
shown below. If you are on Windows, make sure to install the equivalent
version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for
your version of R.

``` r
# checks that you are set up to build R packages from source
if (!requireNamespace("pkgbuild", quietly = TRUE)) {
  install.packages("pkgbuild")
}
pkgbuild::check_build_tools()

# installs the latest isoprocessor2 package (and isoreader2) from GitHub
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("isoverse/isoprocessor2")
```
