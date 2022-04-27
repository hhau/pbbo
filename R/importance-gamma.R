gamma_mixture_importance <- function(
  sample_one,
  sample_two,
  n_internal_importance_draws,
  importance_args,
  ...
) {
  stopifnot(!is.null(importance_args$gamma_sd_multiplier))
  params <- make_params_mix_gamma(
    sample_one = sample_one,
    sample_two = sample_two,
    sd_multiplier = importance_args$gamma_sd_multiplier
  )
  points <- sample_mix_gamma(n = n_internal_importance_draws, params)
  weights <- density_mix_gamma(points, params)
  res <- list(points = points, weights = weights)
  return(res)
}

make_params_mix_gamma <- function(
  sample_one,
  sample_two,
  sd_multiplier
) {
  mean_s1 <- mean(sample_one)
  mean_s2 <- mean(sample_two)
  var_s1 <- min((sd_multiplier ^ 2) * var(sample_one), 1e5)
  var_s2 <- min((sd_multiplier ^ 2) * var(sample_two), 1e5)
  mix_weight <- length(sample_one) / (length(sample_one) + length(sample_two))

  params <- list(
    shape_1 = (mean_s1 ^ 2) / var_s1,
    shape_2 = (mean_s2 ^ 2) / var_s2,
    rate_1 = mean_s1 / var_s1,
    rate_2 = mean_s2 / var_s2,
    mix_weight = mix_weight
  )

  futile.logger::flog.trace(
    paste(
      "mix_gamma params contains",
      paste(names(params), params, sep = " = ", collapse = "; ")
    )
  )

  return(params)
}

density_mix_gamma <- function(x, params, log_scale = FALSE) {
  log_dens_vals <- matrix(data = NA, nrow = length(x), ncol = 2)
  log_dens_vals[, 1] <- log(params$mix_weight) +
    dgamma(
      x = x,
      shape = params$shape_1,
      rate = params$rate_1,
      log = TRUE
    )

  log_dens_vals[, 2] <- log(1 - params$mix_weight) +
    dgamma(
      x = x,
      shape = params$shape_2,
      rate = params$rate_2,
      log = TRUE
    )

  log_dens <- apply(log_dens_vals, 1, matrixStats::logSumExp)

  if (log_scale) {
    return(log_dens)
  } else {
    return(exp(log_dens))
  }
}

sample_mix_gamma <- function(n, params, numerical_lb = 1e-12) {
  samples_per_component <- as.numeric(table(sample(
    x = c(1, 2),
    size = n,
    replace = TRUE,
    prob = c(params$mix_weight, 1 - params$mix_weight)
  )))

  futile.logger::flog.trace(
    "[sample_mix_gamma] samples_per_component = (%d, %d)",
    samples_per_component[1],
    samples_per_component[2]
  )

  points_one <- rgamma(
    n = samples_per_component[1],
    shape = params$shape_1,
    rate = params$rate_1
  )

  if (any(points_one < numerical_lb)) {
    points_one <- resample_gamma(
      sample = points_one,
      shape = params$shape_1,
      rate = params$rate_1,
      numerical_lb = numerical_lb
    )
  }

  points_two <- rgamma(
    n = samples_per_component[2],
    shape = params$shape_2,
    rate = params$rate_2
  )

  if (any(points_two < numerical_lb)) {
    points_two <- resample_gamma(
      sample = points_two,
      shape = params$shape_2,
      rate = params$rate_2,
      numerical_lb = numerical_lb
    )
  }

  res <- c(points_one, points_two)
  return(res)
}

resample_gamma <- function(sample, shape, rate, numerical_lb) {
  futile.logger::flog.info(
    "gamma_mixture resampling triggered. hopefully this just an issue with the
    early stage numerics. If this warning occurs in the late stages of
    optimisation, then results are probably invalid."
  )

  n_target <- length(sample)
  n_invalid <- sum(sample < numerical_lb)
  valid_samples <- sample[sample > numerical_lb]

  extra_samples <- rgamma(
    n = n_invalid * 5,
    shape = shape,
    rate = rate
  )

  valid_extras <- extra_samples[extra_samples > numerical_lb]

  if (length(valid_extras) < n_invalid) {
    stop("Something is very wrong with the numerics for this problem.")
  }

  final_extras <- sample(x = valid_extras, size = n_invalid, replace = FALSE)
  res <- c(valid_samples, final_extras)
  return(res)
}
