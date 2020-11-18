# GCIMS


[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Build Status](https://github.com/Luis-Fernandez-Romero/GCIMS/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/Luis-Fernandez-Romero/GCIMS/actions/)

The goal of `GCIMS` is to offer a data analysis preprocessing pipeline for
GC-IMS samples.

## Installation

`GCIMS` can be installed with the `devtools` package. For this is needed Rtools
and note that it uses packages from CRAN, from BioConductor and from git
repositories:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    

BiocManager::install(c("xx", "xxx"))

devtools::install_github("Luis-Fernandez-Romero/GCIMS")
```


========================

