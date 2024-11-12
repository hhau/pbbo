suppressPackageStartupMessages(library(mlrMBO))

# covariate version
cov_values <- c(-2, 2)

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

test_that("build_approx_kl_discrep passes known good case", {
  approx_kl <- build_approx_kl_discrep(
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    covariate_list = cov_values,
    n_samples_for_approx = 2e6
  )

  val <- approx_kl(c("mu" = 1, "sigma" = 0.25))
  # this should pass approx ((1 - 2e-5) * 100)% of the time
  expect_equal(object = val, expected = 17.61254, tolerance = 0.089)
})
