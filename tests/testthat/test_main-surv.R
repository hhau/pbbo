suppressPackageStartupMessages(library(mlrMBO))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(futile.logger))

flog.threshold(WARN, name = "pbbo")

n_indiv <- 3
n_cov <- 2
censoring_times <- c(1, 2, 3.5)
cov_vals <- rnorm(n = n_cov * n_indiv) %>%
  matrix(nrow = n_indiv)

tar_cdf_cts_portion <- 0.8

# the annoying thing here is that I have to cbind in the censoring times to get
# them to the target_lcdf and etc functions. Annoying. Reenforces that they
# really are covariate information though.

cov_mat <- cbind(censoring_times, cov_vals)

target_lcdf <- function(x, cov_val) {
  res <- array(dim = length(x))
  cens_time <- cov_val["censoring_times"]
  cens_indices <- which(x == cens_time)
  unces_indices <- base::setdiff(seq_len(length(x)), cens_indices)

  res[cens_indices] <- log(1)
  res[unces_indices] <- log((x[unces_indices] / cens_time) * tar_cdf_cts_portion)

  return(res)
}

target_lpdf <- function(x, cov_val) {
  res <- array(dim = length(x))
  cens_time <- cov_val["censoring_times"]
  cens_indices <- which(x == cens_time)
  unces_indices <- base::setdiff(seq_len(length(x)), cens_indices)

  res[cens_indices] <- log(1 - tar_cdf_cts_portion)
  res[unces_indices] <- log(tar_cdf_cts_portion * (1 / cens_time))

  return(res)
}

target_sampler <- function(n, cov_val) {
  cens_time <- as.numeric(cov_val["censoring_times"])
  cts_or_cens_indicators <- sample(
    x = 1 : 2,
    size = n,
    prob = c(tar_cdf_cts_portion, 1 - tar_cdf_cts_portion),
    replace = TRUE
  ) %>%
     tabulate(nbins = 2)

  cts_times <- runif(n = cts_or_cens_indicators[1], min = 0, max = cens_time)
  cens_times <- rep(cens_time, cts_or_cens_indicators[2])

  res <- c(cts_times, cens_times)
  return(res)
}

prior_predictive_sampler <- function(n, lambda, cov_val) {
  cens_time <- as.numeric(cov_val[1])
  local_cov_vec <- cov_val[-1]
  haz_gamma <- rgamma(n = n, shape = lambda["gamma_shape"], rate = lambda["gamma_rate"])
  theta_zero <- rnorm(n = n, mean = lambda["theta_zero_loc"], sd = lambda["theta_zero_scale"])
  theta_loc_names <- names(lambda) %>%
    grep("^theta_loc[0-9+]$", x = ., value = TRUE)

  theta_scale_names <- names(lambda) %>%
    grep("^theta_scale[0-9+]$", x = ., value = TRUE)

  theta_mat <- rnorm(
    n = n * length(local_cov_vec),
    mean = rep(lambda[theta_loc_names], each = n),
    sd = rep(lambda[theta_scale_names], each = n)
  ) %>%
    matrix(nrow = n, ncol = length(local_cov_vec))

  neg_log_u <- -log(runif(n = n))
  cov_term <- as.numeric(theta_zero + (local_cov_vec %*% t(theta_mat)))
  event_times <- exp((1 / haz_gamma) * (log(neg_log_u) - cov_term))
  cens_indicators <- which(event_times >= cens_time)
  event_times[cens_indicators] <- cens_time
  return(event_times)
}

param_set <- makeParamSet(
  makeNumericParam(id = "gamma_shape", default = 1, lower = 0, upper = 10),
  makeNumericParam(id = "gamma_rate", default = 1, lower = 0, upper = 10),
  makeNumericParam(id = "theta_zero_loc", default = 0, lower = -10, upper = 10),
  makeNumericParam(id = "theta_zero_scale", default = 1, lower = 0, upper = 10),
  makeNumericVectorParam(id = "theta_loc", len = n_cov, lower = -10, upper = 10, default = rep(0, n_cov)),
  makeNumericVectorParam(id = "theta_scale", len = n_cov, lower = 0, upper = 10, default = rep(1, n_cov))
)

test_that("surv example can run error free", {
  pbbo_res <- suppressWarnings(pbbo(
    target_lcdf = target_lcdf,
    target_lpdf = target_lpdf,
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    param_set = param_set,
    covariate_values = cov_mat,
    discrepancy = "log_cvm",
    n_crs2_iters = 50,
    importance_method = "surv_mixture",
    importance_args = list(
      surv_mixture_sd_multiplier = 1.05,
      surv_mixture_cont_frac = 0.95,
      surv_mixture_cens_times = censoring_times
    ),
    bayes_opt_iters_per_batch = 10,
    n_internal_prior_draws = 20,
    n_internal_importance_draws = 10
  ))

  expect_s3_class(pbbo_res[[1]], "MBOSingleObjResult")
})
