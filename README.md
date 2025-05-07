
# R Package "spfda" - Sparse function-on-scalar regression

<!-- badges: start -->
[![CRAN-version](https://www.r-pkg.org/badges/version/spfda)](https://CRAN.R-project.org/package=spfda)
[![license-MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](https://github.com/dipterix/spfda/blob/master/LICENSE)
<!-- [![Cran-version](http://cranlogs.r-pkg.org/badges/grand-total/spfda)](https://CRAN.R-project.org/package=spfda) -->
[![R-check](https://github.com/dipterix/spfda/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dipterix/spfda/actions/workflows/R-CMD-check.yaml)
[![DOI-YAEL](https://img.shields.io/badge/DOI-10.1111%2Fbiom.13684-blue?link=https%3A%2F%2Fdoi.org%2F10.1111%2Fbiom.13684)](https://doi.org/10.1111/biom.13684)
<!-- badges: end -->

> This package implements paper "Functional Group Bridge for Simultaneous Regression and Support Estimation" 

[[Preprint](https://arxiv.org/abs/2006.10163)] [[Github](https://github.com/dipterix/spfda)] [[CRAN](https://cran.r-project.org/package=spfda)] [[Examples](https://doi.org/10.5281/zenodo.6363319)]

![Demo example](https://raw.githubusercontent.com/dipterix/spfda/master/inst/cover.png)

## Installation

You can install the **released** version of `spfda` ![CRAN-version](https://www.r-pkg.org/badges/version/spfda) from [CRAN](https://CRAN.R-project.org/package=spfda) with:

``` r
install.packages("spfda")
```

The **experimental** version can be installed via

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

Use `citation('spfda')` to generate citation information, or check [this link](https://doi.org/10.1111/biom.13684).

```
Zhengjia Wang, John Magnotti, Michael S. Beauchamp, Meng Li, 
  Functional Group Bridge for Simultaneous Regression and Support Estimation, 
  Biometrics, Volume 79, Issue 2, June 2023, Pages 1226–1238, 
  https://doi.org/10.1111/biom.13684
```

