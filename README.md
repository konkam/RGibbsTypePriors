
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RGibbsTypePriors

<!-- badges: start -->

<!-- badges: end -->

Computing clusters prior distribution for Gibbs-type processes.

## Introduction


The following reference gives an overview of Gibbs-type priors and their importance for Bayesian Nonparametrics:

De Blasi, Pierpaolo, Stefano Favaro, Antonio Lijoi, Ramsés H. Mena, Igor Prünster, and Matteo Ruggiero. “Are Gibbs-Type Priors the Most Natural Generalization of the Dirichlet Process?” IEEE Transactions on Pattern Analysis and Machine Intelligence 37, no. 2 (2015): 212–29. https://doi.org/10.1109/TPAMI.2013.217.


An application of the functions implemented in this package was presented here:

Bystrova, D., Arbel, J., Kon Kam King, G., Deslandes, F. (2021) Approximating the clusters' prior distribution in Bayesian nonparametric models,  Third Symposium on Advances in Approximate Bayesian Inference, https://openreview.net/forum?id=J0SSW5XeWUY

Please cite it as a reference if you use the package.


## Installation


You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("konkam/RGibbsTypePriors")
```


The package depends on the [Arb](https://arblib.org/) library by F. Johansson:

F. Johansson. “Arb: efficient arbitrary-precision midpoint-radius interval arithmetic”, IEEE Transactions on Computers, 66(8):1281-1292, 2017. DOI: 10.1109/TC.2017.2690633.
