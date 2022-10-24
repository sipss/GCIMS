# GCIMS


[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Build Status](https://github.com/sipss/GCIMS/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/sipss/GCIMS/actions/)

The goal of `GCIMS` is to offer a data analysis preprocessing pipeline for
GC-IMS samples.

## Installation

`GCIMS` can be installed with the `remotes` package. Ensure you have it installed
with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
```

If you have cloned this repository in a local folder, you can run from that folder:

``` r
remotes::install_local(".")
```

Alternatively, you can clone and install it from github. Since this repository is 
for now private, you will need a github
[personal access token](https://github.com/settings/tokens) to access it:

```r
# Create your gh_PAT at: https://github.com/settings/tokens
gh_PAT <- "your-personal-access-token"
remotes::install_github("sipss/GCIMS", auth_token = gh_PAT)
```

Once the repository becomes public, you can install the package with:

```r
remotes::install_github("sipss/GCIMS") # only if the repository is public
```

