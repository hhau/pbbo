target_cdf <- function(x) {
  pnorm(x, mean = 2, sd = 0.5)
}

target_sampler <- function(n) {
  rnorm(n = n, mean = 2, sd = 0.5)
}

prior_predictive_sampler <- function(n, lambda) {
  rnorm(n = n, mean = lambda['mu'], sd = lambda['sigma'])
}

test_lambda <- c(mu = 1, sigma = 2.5)

local_importance <- pbbo:::uniform_importance
internal_discrepancy_f <- pbbo:::cvm_discrepancy

## baseline
test_discrepancy <- build_discrep(
  target_cdf = target_cdf,
  target_sampler = target_sampler,
  prior_predictive_sampler = prior_predictive_sampler,
  internal_discrepancy_f = internal_discrepancy_f,
  n_internal_prior_draws = 500,
  importance_method = 'uniform',
  n_internal_importance_draws = 200
)

test_val <- test_discrepancy(lambda_mlrform = test_lambda)

test_that('cvm discrep + uniform importance yeilds sensible result', {
  expect_true(is.numeric(test_val))
  expect_true((test_val > 0) & (test_val < 5e3))
})
