
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `pbbo` – Prior by Bayesian Optimisation

<!-- badges: start -->

[![R-CMD-check](https://github.com/hhau/pbbo/workflows/R-CMD-check/badge.svg)](https://github.com/hhau/pbbo/actions)
<!-- badges: end -->

`pbbo` uses information you supply about the prior predictive
distribution to help you find a prior for the parameters in your
Bayesian model.

## Installation

You can install the development version of pbbo from
[GitHub](https://github.com/hhau/pbbbo) with:

``` r
# install.packages("devtools")
devtools::install_github("hhau/pbbo")
```

## An example

Suppose the target distribution is `N(0, 0.5^2)`. We define the target
LCDF and function to draws samples from the target:

``` r
library(pbbo)
suppressPackageStartupMessages(library(ParamHelpers))

target_lcdf <- function(x) {
  Rmpfr::pnorm(x, mean = 2, sd = 0.5, log.p = TRUE)
}

target_sampler <- function(n) {
  rnorm(n = n, mean = 2, sd = 0.5)
}
```

and define a the prior predictive distribution of a model (and its
hyperparameters):

``` r
prior_predictive_sampler <- function(n, lambda) {
  rnorm(n = n, mean = lambda["mu"], sd = lambda["sigma"])
}

param_set <- makeParamSet(
  makeNumericParam(id = "mu", default = 0.2, lower = -50, upper = 50),
  makeNumericParam(id = "sigma", lower = 0, upper = 20, default = 0.2)
)
```

We can then call `pbbo` to get optimal values of λ = (μ, σ)

``` r
library(futile.logger)

flog.threshold(futile.logger::INFO, name = "pbbo")
#> NULL
flog.appender(
  appender.file("../../logs/pbbo-test.log"),
  "pbbo"
)
#> NULL

pbbo_res <- suppressWarnings(pbbo(
  target_lcdf = target_lcdf,
  target_sampler = target_sampler,
  prior_predictive_sampler = prior_predictive_sampler,
  param_set = param_set,
  n_crs2_iters = 300,
  n_internal_prior_draws = 5e3,
  n_internal_importance_draws = 1e3,
  importance_method = "uniform",
  bayes_opt_batches = 1,
  bayes_opt_iters_per_batch = 50,
  bayes_opt_design_points_per_batch = 40
))
#> INFO [2022-06-08 13:57:01] Starting stage one CRS2 optimiser
#> INFO [2022-06-08 13:57:06] Starting Bayes opt batch 1

opt_lambda <- pbbo_res[[1]]$x
print(opt_lambda)
#> $mu
#> [1] 2.001372
#> 
#> $sigma
#> [1] 0.5057144
```

We can compare the prior predictive distribution at the optima against
the target:

``` r
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

n_samples <- 2e3

plot_tbl <- tibble(
  x = c(
    prior_predictive_sampler(n_samples, unlist(opt_lambda)),
    target_sampler(n_samples)
  ),
  grp = rep(c("optima", "target"), each = n_samples)
)

ggplot(plot_tbl) +
  geom_line(
    mapping = aes(x = x, colour = grp),
    stat = "density"
  ) +
  labs(colour = "Type")
```

<img src="man/figures/README-compare-1.png" width="100%" />
