suppressPackageStartupMessages(library(mlrMBO))

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
  makeNumericParam(id = "sigma", lower = 0.2, upper = 20, default = 0.5)
)

test_discrepancy <- pbbo:::build_discrep_covariate(
  target_lcdf = target_lcdf,
  target_sampler = target_sampler,
  prior_predictive_sampler = prior_predictive_sampler,
  covariate_list = as.list(cov_values),
  internal_discrepancy_f = pbbo:::log_cvm_discrepancy,
  n_internal_prior_draws = 500,
  importance_method = "student_t_mixture",
  importance_args = list(
    uniform_lower = NULL,
    uniform_upper = NULL,
    gamma_sd_multiplier = 1.05,
    student_t_sd_multiplier = 1.05,
    surv_mixture_sd_multiplier = 1.05,
    surv_mixture_cont_frac = 0.95,
    surv_mixture_cens_times = rep(NULL, length(cov_values))
  ),
  n_internal_importance_draws = 200
)

t_res <- pbbo:::design_from_crs2(
  discrep = test_discrepancy,
  param_set = param_set,
  n_design = 20,
  n_crs2_iters = 33
)

test_that("crs2_desgin returns a tibble and does not error", {
  expect_s3_class(
    object = t_res,
    "tbl"
  )

  expect_equal(
    object = nrow(t_res),
    expected = 20
  )
})
