suppressPackageStartupMessages(library(mlrMBO))

# necessary things
target_lcdf <- function(x) {
  pnorm(x, mean = 2, sd = 0.5, log.p = TRUE)
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

extra_term <- function(lambda) {
  mean_log_sd <- log(sqrt(lambda["sigma"] ^ 2))
  return(mean_log_sd)
}

## actual tests
test_that("main with multiple objectives can run error free", {
  mlr_res <- suppressWarnings(pbbo(
    model_name = "test_normal",
    target_lcdf = target_lcdf,
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    discrepancy = "log_cvm",
    n_crs2_iters = 20,
    param_set = param_set,
    n_internal_prior_draws = 20,
    n_internal_importance_draws = 10,
    bayes_opt_iters_per_batch = 5,
    bayes_opt_print = FALSE,
    extra_objective_term = extra_term
  ))

  expect_s3_class(mlr_res[[1]], class = "MBOSingleObjResult")
})


test_that("main with multiple objectives _and batches_ can run error free", {
  mlr_res <- suppressWarnings(pbbo(
    model_name = "test_normal",
    target_lcdf = target_lcdf,
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    discrepancy = "log_cvm",
    n_crs2_iters = 20,
    param_set = param_set,
    n_internal_prior_draws = 20,
    n_internal_importance_draws = 10,
    bayes_opt_batches = 2,
    bayes_opt_iters_per_batch = 5,
    bayes_opt_print = FALSE,
    extra_objective_term = extra_term
  ))

  expect_s3_class(mlr_res[[1]], class = "MBOSingleObjResult")
})
