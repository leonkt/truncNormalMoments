# truncNormalMoments

<!-- badges: start -->
<!-- badges: end -->

`truncNormalMoments` provides functionality to estimates mean and standard
deviation for a truncated normal distribution.

Specifically, this package finds the MAP parameter estimates for the mean and
standard deviation for a normal distribution, given data from the truncated
distribution. The method used extends Bayesian methods for parameter estimation
for a singly truncated normal distribution [see Xiaoping Zhou, Rosella
Giacometti, Frank J. Fabozzi & Ann H. Tucker (2014). Bayesian estimation of
truncated data with applications to operational risk measurement, Quantitative
Finance, 14:5, 863-888, DOI: 10.1080/14697688.2012.752103]. We extend it by
considering the doubly truncated normal distribution. We use the Jeffreys prior
and find MAP estimates of the mean and standard deviation of a doubly truncated
normal, with left and right truncations specified beforehand.

## Installation

You can install the development version of truncNormalMoments from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("leonkt/truncNormalMoments")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(truncNormalMoments)

# generate data from a truncated normal
x <- truncnorm::rtruncnorm(100, a = -1, b = 2, mean = 0.5, sd = 0.5)

# estimate its paramaeters
trunc_est(x, a = -1, b = 2)
```
