library(mlrMBO)

# covariate version
cov_values <- c(-2, 2)

target_lcdf <- function(x, cov) {
  pnorm(x, mean = cov, sd = 0.5, log.p = TRUE)
}

target_sampler <- function(n, cov) {
  rnorm(n = n, mean = cov, sd = 0.5)
}

prior_predictive_sampler <- function(n, lambda, cov) {
  rnorm(n = n, mean = lambda["mu"] + cov, sd = lambda["sigma"])
}

param_set <- makeParamSet(
  makeNumericParam(id = "mu", default = 0.2, lower = -50, upper = 50),
  makeNumericParam(id = "sigma", lower = 0, upper = 20, default = 0.2)
)

test_that("main function with covariates can run error free", {
  mlr_res <- suppressWarnings(
    pbbo(
      model_name = "test_normal_covariate",
      target_lcdf = target_lcdf,
      target_sampler = target_sampler,
      prior_predictive_sampler = prior_predictive_sampler,
      covariate_values = cov_values,
      discrepancy = "log_cvm",
      param_set = param_set,
      n_internal_prior_draws = 20,
      n_internal_importance_draws = 50,
      bayes_opt_iters_per_batch = 3,
      bayes_opt_print = FALSE
    )
  )

  expect_s3_class(mlr_res[[1]], class = "MBOSingleObjResult")
})

test_that("bad covariate args cause an error", {
  bad_covariate_form <- list(x1 = 1, x2 = c(2, 3))
  expect_error(
    pbbo(
      model_name = "test_normal_covariate",
      target_lcdf = target_lcdf,
      target_sampler = target_sampler,
      prior_predictive_sampler = prior_predictive_sampler,
      covariate_values = bad_covariate_form,
      discrepancy = "log_cvm",
      param_set = param_set,
      n_internal_prior_draws = 20,
      n_internal_importance_draws = 50,
      bayes_opt_iters_per_batch = 3,
      bayes_opt_print = FALSE
    )
  )
})
