# GCIMS <a href="https://sipss.github.io/GCIMS"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Build Status](https://github.com/sipss/GCIMS/workflows/R-CMD-check-bioc/badge.svg?branch=master)](https://github.com/sipss/GCIMS/actions/)
[![Documentation](https://img.shields.io/badge/documentation-pkgdown-informational)](https://sipss.github.io/GCIMS/)

<!-- badges: end -->

GCIMS is an R package implementing a preprocessing pipeline for Gas Chromatography – Ion
Mobility Spectrometry samples.

## Installation

`GCIMS` can be installed with the `remotes` package. Ensure you have it installed
with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
```

You can then install the package with:

```r
remotes::install_github("sipss/GCIMS")
```

Checkout our [Introduction to GCIMS](https://sipss.github.io/GCIMS/articles/introduction-to-gcims.html) to start.
