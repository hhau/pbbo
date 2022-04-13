# This might be a poor idea? As it effectively truncates the integral to
# the values spanned by [lower, upper].
uniform_importance <- function(
  sample_one,
  sample_two,
  n_internal_importance_draws,
  lower,
  upper,
  ...
) {
  if (is.null(lower)) {
    lower <- min(c(sample_one, sample_two))
  }

  if (is.null(upper)) {
    upper <- max(c(sample_one, sample_two))
  }

  stopifnot(lower < upper)

  diff <- upper - lower
  points <- runif(n = n_internal_importance_draws, min = lower, max = upper)
  weights <- rep(1 / diff, times = n_internal_importance_draws)
  res <- list(points = points, weights = weights)
  return(res)
}

# fit a two component mixture model (of t_5 densities)
# evaluate fitted mixture using khat statistics
student_t_mixture_importance <- function(
  sample_one,
  sample_two,
  n_internal_importance_draws,
  ...
) {
  params <- make_params_mix_student_t(sample_one, sample_two)
  points <- sample_mix_student_t(n = n_internal_importance_draws, params)
  weights <- density_mix_student_t(points, params)
  res <- list(points = points, weights = weights)
  return(res)
}

make_params_mix_student_t <- function(
  sample_one,
  sample_two,
  sd_multiplier = 1.05
) {
  mean_s1 <- mean(sample_one)
  mean_s2 <- mean(sample_two)
  sd_s1 <- sd_multiplier * sd(sample_one)
  sd_s2 <- sd_multiplier * sd(sample_two)
  mix_weight <- length(sample_one) / (length(sample_one) + length(sample_two))

  params <- list(
    loc_1 = mean_s1,
    loc_2 = mean_s2,
    scale_1 = sd_s1,
    scale_2 = sd_s2,
    mix_weight = mix_weight
  )

  return(params)
}

density_mix_student_t <- function(x, params, log_scale = FALSE) {
  log_dens_vals <- matrix(data = NA, nrow = length(x), ncol = 2)

  # log jacobian term here necessary!!
  log_dens_vals[, 1] <- log(params$mix_weight) +
    dt(
      x = (x - params$loc_1) / params$scale_1,
      df = 5,
      log = TRUE
    ) -
      log(params$scale_1)

  log_dens_vals[, 2] <- log(1 - params$mix_weight) + dt(
    x = (x - params$loc_2) / params$scale_2,
    df = 5,
    log = TRUE
  ) -
    log(params$scale_2)

  log_dens <- apply(log_dens_vals, 1, matrixStats::logSumExp)

  if (log_scale) {
    return(log_dens)
  } else {
    return(exp(log_dens))
  }
}

sample_mix_student_t <- function(n, params) {
  samples_per_component <- table(sample(
    x = c(1, 2),
    size = n,
    replace = TRUE,
    prob = c(params$mix_weight, 1 - params$mix_weight)
  ))

  points_one <- params$loc_1 +
    rt(n = samples_per_component[1], df = 5) * params$scale_1

  points_two <- params$loc_2 +
    rt(n = samples_per_component[2], df = 5) * params$scale_2

  res <- c(points_one, points_two)
  return(res)
}

# same as student_t but with gammas? enforce ordering somehow.
# easiest to do with stan.
gamma_mixture_importance <- function(
  sample_one,
  sample_two,
  n_internal_importance_draws,
  ...
) {
  params <- make_params_mix_gamma(sample_one, sample_two)
  points <- sample_mix_gamma(n = n_internal_importance_draws, params)
  weights <- density_mix_gamma(points, params)
  res <- list(points = points, weights = weights)
  return(res)
}

make_params_mix_gamma <- function(
  sample_one,
  sample_two,
  sd_multiplier = 1.05
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

density_mix_gamma <- function(
  x,
  params,
  log_scale = FALSE
) {
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

sample_mix_gamma <- function(
  n,
  params
) {
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

  points_two <- rgamma(
    n = samples_per_component[2],
    shape = params$shape_2,
    rate = params$rate_2
  )

  res <- c(points_one, points_two)
  return(res)
}
