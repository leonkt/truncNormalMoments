# truncnormbayes

<!-- badges: start -->
<!-- badges: end -->

`truncnormbayes` provides functionality to estimate mean and standard
deviation for a truncated normal distribution.

Specifically, this package finds the posterior modes for the mean and standard
deviation for a truncated normal distribution with one or two known truncation
points. The method used extends Bayesian methods for parameter estimation for a
singly truncated normal distribution under the Jeffreys prior [see Zhou, X.,
Giacometti, R., Fabozzi, F. J., & Tucker, A. H. (2014). Bayesian estimation of
truncated data with applications to operational risk measurement. Quantitative
Finance, 14(5), 863-888. [DOI 10.1080/14697688.2012.752103](https://doi.org/10.1080/14697688.2012.752103)]. This package
additionally allows for a doubly truncated normal distribution.

## Installation

You can install the development version of truncnormbayes from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("leonkt/truncnormbayes")
```

## Example

``` r
library(truncnormbayes)

# generate data from a truncated normal
x <- truncnorm::rtruncnorm(100, a = -1, b = 2, mean = 0.5, sd = 0.5)

# estimate its paramaeters
trunc_est(x, a = -1, b = 2)
```
