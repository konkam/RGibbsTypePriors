---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# RGibbsTypePriors

<!-- badges: start -->
<!-- badges: end -->

The goal of RGibbsTypePriors is to Compute clusters prior distribution for Gibbs-type processes


## Installation

You can install the released version of RGibbsTypePriors from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RGibbsTypePriors")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("konkam/RGibbsTypePriors")
```
## Example

This is a basic example which shows you how to solve a common problem:


```r
library(RGibbsTypePriors)
res=pkn_ngg(1:100, 100, 1.2, 0.8, prec=1000)
plot(res$k, as.numeric(res$pkn))
```

<img src="man/figures/README-example-1.png" title="plot of chunk example" alt="plot of chunk example" width="100%" />


