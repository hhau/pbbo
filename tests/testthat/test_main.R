library(mlrMBO)

# necessary things
target_cdf <- function(x) {
  pnorm(x, mean = 2, sd = 0.5)
}

target_sampler <- function(n) {
  rnorm(n = n, mean = 2, sd = 0.5)
}

prior_predictive_sampler <- function(n, lambda) {
  rnorm(n = n, mean = lambda['mu'], sd = lambda['sigma'])
}

param_set <- makeParamSet(
  makeNumericParam(id = 'mu', default = 0.2, lower = -50, upper = 50),
  makeNumericParam(id = 'sigma', lower = 0, upper = 20, default = 0.2)
)

## actual tests
test_that('main function can run error free', {
  mlr_res <- suppressWarnings(pbbo(
    model_name = 'test_normal',
    target_cdf = target_cdf,
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    discrepancy = 'cvm',
    param_set = param_set,
    n_internal_prior_draws = 20,
    n_internal_importance_draws = 10,
    bayes_opt_iters = 5,
    bayes_opt_print = FALSE
  ))

  expect_s3_class(mlr_res, class = 'MBOSingleObjResult')
})

test_that('bad args give errors', {
  expect_error(suppressWarnings(pbbo(
    model_name = 'test_normal',
    target_cdf = target_cdf,
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    discrepancy = 'definitely_not_a_discrep',
    param_set = param_set
  )))
})
