# GCIMS


[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Build Status](https://github.com/sipss/GCIMS/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/sipss/GCIMS/actions/)

The goal of `GCIMS` is to offer a data analysis preprocessing pipeline for
GC-IMS samples.

## Installation

`GCIMS` can be installed with the `remotes` package.
It uses depebdencies from CRAN and BioConductor repositories:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
    

remotes::install_github("sipss/GCIMS")
```


========================

