
# R Package "spfda" - Sparse function-on-scalar regression

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/dipterix/spfda.svg?branch=master)](https://travis-ci.com/dipterix/spfda)
[![CRAN-version](https://www.r-pkg.org/badges/version/spfda)](https://CRAN.R-project.org/package=spfda)
[![license-MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](https://github.com/dipterix/spfda/blob/master/LICENSE)
<!-- [![Cran-version](http://cranlogs.r-pkg.org/badges/grand-total/spfda)](https://CRAN.R-project.org/package=spfda) -->
<!-- badges: end -->

> This package implements paper "Functional Group Bridge for Simultaneous Regression and Support Estimation" 

[[PDF](https://arxiv.org/abs/2006.10163)] [[Github](https://github.com/dipterix/spfda)] [[CRAN](https://cran.r-project.org/package=spfda)]

![Demo example](https://raw.githubusercontent.com/dipterix/spfda/master/docs/cover.png)

## Installation

You can install the released version of `spfda` from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("spfda")
```

The development version can be installed via

``` r
# install.packages("remotes")
remotes::install_github("dipterix/spfda")
```

## Example

``` r
library(spfda)
dat <- spfda_simulate()
x <- dat$X
y <- dat$Y

## basic example code

fit <- spfda(y, x, lambda = 5, CI = TRUE)

## Generics

BIC(fit)
plot(fit, col = c("orange", "dodgerblue3", "darkgreen"),
     main = "Fitted with 95% CI", aty = c(0, 0.5, 1), atx = c(0,0.2,0.8,1))
matpoints(fit$time, t(dat$env$beta), type = 'l', col = 'black', lty = 2)
legend('topleft', c("Fitted", "Underlying"), lty = c(1,2), bty = 'n')
print(fit)
coefficients(fit)
```

## Citation

Use `citation('spfda')` to generate citation information.

```
Wang, Z, Magnotti, JF, Beauchamp, MS. Li, M, Functional Group Bridge for
  Simultaneous Regression and Support Estimation.
```
