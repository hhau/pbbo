student_t_mixture_importance <- function(
  sample_one,
  sample_two,
  n_internal_importance_draws,
  importance_args,
  ...
) {
  stopifnot(!is.null(importance_args$student_t_sd_multiplier))
  params <- make_params_mix_student_t(
    sample_one = sample_one,
    sample_two = sample_two,
    sd_multiplier = importance_args$student_t_sd_multiplier
  )
  points <- sample_mix_student_t(n = n_internal_importance_draws, params)
  weights <- density_mix_student_t(points, params)
  res <- list(points = points, weights = weights)
  return(res)
}

make_params_mix_student_t <- function(
  sample_one,
  sample_two,
  sd_multiplier
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
