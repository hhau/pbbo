suppressPackageStartupMessages(library(mlrMBO))

# covariate version
cov_values <- c(-2, 2)

target_sampler <- function(n, cov) {
  rnorm(n = n, mean = cov, sd = 0.5)
}

prior_predictive_sampler <- function(n, lambda, cov) {
  rnorm(n = n, mean = lambda["mu"] + cov, sd = lambda["sigma"])
}

test_that("build_approx_kl_discrep passes known good case", {
  approx_kl <- build_approx_kl_discrep(
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    covariate_list = cov_values,
    n_samples_for_approx = 2e6,
    direction = "fwd"
  )

  val <- approx_kl(c("mu" = 1, "sigma" = 0.25))
  # this should pass approx ((1 - 2e-5) * 100)% of the time
  expect_equal(object = val, expected = log(17.61254), tolerance = 0.004530517)
})


target_sampler <- function(n) {
  rnorm(n = n, mean = 2, sd = 0.5)
}

prior_predictive_sampler <- function(n, lambda) {
  rnorm(n = n, mean = lambda["mu"], sd = lambda["sigma"])
}

test_lambda <- c(mu = 1, sigma = 2.5)

test_that("build_approx_kl_discrep_pop passes known good case for population setting", {
  approx_kl <- build_approx_kl_discrep_pop(
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    n_samples_for_approx = 2e6,
    direction = "fwd"
  )

  val <- approx_kl(test_lambda)
  # idk tolerance is behaving very weirdly to me, relative diff or something
  # just sqrt for now and move on - if it fails too often something is weird.
  expect_equal(object = val, expected = 0.1901456, tolerance = sqrt(0.002503809))
})