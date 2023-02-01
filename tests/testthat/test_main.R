suppressPackageStartupMessages(library(mlrMBO))

# necessary things
target_lcdf <- function(x) {
  pnorm(x, mean = 2, sd = 0.5, log.p = TRUE)
}

target_lpdf <- function(x) {
  dnorm(x, mean = 2, sd = 0.5, log = TRUE)
}

target_sampler <- function(n) {
  rnorm(n = n, mean = 2, sd = 0.5)
}

prior_predictive_sampler <- function(n, lambda) {
  rnorm(n = n, mean = lambda["mu"], sd = lambda["sigma"])
}

param_set <- makeParamSet(
  makeNumericParam(id = "mu", default = 0.2, lower = -50, upper = 50),
  makeNumericParam(id = "sigma", lower = 0, upper = 20, default = 0.2)
)

## actual tests
test_that("main function can run error free", {
  mlr_res <- suppressWarnings(pbbo(
    model_name = "test_normal",
    target_lcdf = target_lcdf,
    target_lpdf = target_lpdf,
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    discrepancy = "log_cvm",
    n_crs2_iters = 20,
    param_set = param_set,
    n_internal_prior_draws = 20,
    n_internal_importance_draws = 10,
    bayes_opt_iters_per_batch = 5,
    bayes_opt_print = FALSE
  ))

  expect_s3_class(mlr_res[[1]], class = "MBOSingleObjResult")
})

test_that("bad args give errors", {
  expect_error(suppressWarnings(pbbo(
    model_name = "test_normal",
    target_lcdf = target_lcdf,
    target_lpdf = target_lpdf,
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    discrepancy = "definitely_not_a_discrep",
    param_set = param_set
  )))
})

test_that("batching can run error free", {
  pbbo_res <- suppressWarnings(pbbo(
    model_name = "test_normal",
    target_lcdf = target_lcdf,
    target_lpdf = target_lpdf,
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    discrepancy = "log_cvm",
    param_set = param_set,
    n_crs2_iters = 10,
    n_internal_prior_draws = 40,
    n_internal_importance_draws = 10,
    bayes_opt_batches = 2,
    bayes_opt_iters_per_batch = 5,
    bayes_opt_print = FALSE
  ))

  expect_s3_class(pbbo_res[[2]], class = "MBOSingleObjResult")
})
