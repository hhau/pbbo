# `pbbo` -- Prior by Bayesian Optimisation

<!-- badges: start -->
[![R-CMD-check](https://github.com/hhau/pbbo/workflows/R-CMD-check/badge.svg)](https://github.com/hhau/pbbo/actions)
<!-- badges: end -->

This package accompanies a forthcoming paper on translation plausible prior predictive distributions in a prior for parameters.
The intention is to acquire an optimal set of hyperparameters λ by matching the prior predictive distribution from some model to a "target" prior predictive distribution.

# An example

Suppose the target distribution is `N(0, 0.5^2)`.
We define the target LCDF and function to draws samples from the target:

```{r}
target_lcdf <- function(x) {
  pnorm(x, mean = 2, sd = 0.5, log.p = TRUE)
}

target_sampler <- function(n) {
  rnorm(n = n, mean = 2, sd = 0.5)
}
```
and define a the prior predictive distribution of a model (and its hyperparameters):

```{r}
prior_predictive_sampler <- function(n, lambda) {
  rnorm(n = n, mean = lambda["mu"], sd = lambda["sigma"])
}

param_set <- makeParamSet(
  makeNumericParam(id = "mu", default = 0.2, lower = -50, upper = 50),
  makeNumericParam(id = "sigma", lower = 0, upper = 20, default = 0.2)
)
```

We can then call `pbbo` to get optimal values of λ = (μ, σ)

```{r}
pbbo_res <- pbbo(
  target_lcdf = target_lcdf,
  target_sampler = target_sampler,
  prior_predictive_sampler = prior_predictive_sampler,
  param_set = param_set,
  n_internal_prior_draws = 5e3,
  n_internal_importance_draws = 1e3
)
```
