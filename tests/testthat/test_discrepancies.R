target_lcdf <- function(x) {
  pnorm(as.numeric(x), mean = 2, sd = 0.5, log.p = TRUE)
}

target_sampler <- function(n) {
  rnorm(n = n, mean = 2, sd = 0.5)
}

prior_predictive_sampler <- function(n, lambda) {
  rnorm(n = n, mean = lambda["mu"], sd = lambda["sigma"])
}

test_lambda <- c(mu = 1, sigma = 2.5)
internal_discrepancy_f <- pbbo:::log_cvm_discrepancy

## baseline
test_discrepancy <- build_discrep(
  target_lcdf = target_lcdf,
  target_sampler = target_sampler,
  prior_predictive_sampler = prior_predictive_sampler,
  internal_discrepancy_f = internal_discrepancy_f,
  n_internal_prior_draws = 500,
  importance_method = "uniform",
  importance_args = list(
    uniform_lower = NULL,
    uniform_upper = NULL
  ),
  n_internal_importance_draws = 200
)

test_val <- test_discrepancy(lambda_mlrform = test_lambda)

test_that("log_cvm discrep + uniform importance yields sensible result", {
  expect_true(is.numeric(test_val))
  expect_true((test_val > -5e3) & (test_val < 5e3))
})

# test the log_ad discrep
ad_test_target_lcdf <- function(x) {
  Rmpfr::pnorm(x, mean = 2, sd = 0.5, log.p = TRUE)
}

ad_test_target_sampler <- function(n) {
  rnorm(n, mean = 2, sd = 0.5)
}

ad_test_discrepancy <- build_discrep(
  target_lcdf = ad_test_target_lcdf,
  target_sampler = ad_test_target_sampler,
  prior_predictive_sampler = prior_predictive_sampler,
  internal_discrepancy_f = pbbo:::log_ad_discrepancy,
  n_internal_prior_draws = 5e3,
  importance_method = "uniform",
  importance_args = list(
    uniform_lower = NULL,
    uniform_upper = NULL
  ),
  n_internal_importance_draws = 500
)

ad_test_val <- ad_test_discrepancy(c(mu = 0, sigma = 6))
test_that("log_ad discrep + uniform importance yields sensible result", {
  expect_true(!is.nan(ad_test_val))
  expect_true(is.numeric(ad_test_val))
})
