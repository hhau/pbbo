#' @importFrom stats rbeta dbeta var sd
surv_mixture_importance <- function(
  sample_one,
  sample_two,
  n_internal_importance_draws,
  importance_args,
  ...
) {
  stopifnot(
    !is.null(importance_args$surv_mixture_sd_multiplier),
    !is.null(importance_args$surv_mixture_cont_frac),
    !is.null(importance_args$censoring_time)
  )

  sd_multiplier <- importance_args$surv_mixture_sd_multiplier
  continuous_fraction <- importance_args$surv_mixture_cont_frac
  censoring_time <- importance_args$censoring_time

  params <- make_params_surv_mix(
    sample_one = sample_one,
    sample_two = sample_two,
    sd_multiplier = sd_multiplier,
    continuous_fraction = continuous_fraction,
    censoring_time = censoring_time
  )

  points <- sample_surv_mix(
    n = n_internal_importance_draws,
    params = params,
    censoring_time = censoring_time
  )

  weights <- density_surv_mix(points, params, censoring_time)
  res <- list(points = points, weights = weights)
  return(res)
}

make_params_surv_mix <- function(
  sample_one,
  sample_two,
  sd_multiplier,
  continuous_fraction,
  censoring_time
) {
  scaled_sample_one <- sample_one / censoring_time
  scaled_sample_two <- sample_two / censoring_time

  uncens_indices_one <- which(!(sample_one == censoring_time))
  uncens_indices_two <- which(!(sample_two == censoring_time))

  pars_one <- make_beta_params(
    scaled_sample_one[uncens_indices_one],
    sd_multiplier
  )

  pars_two <- make_beta_params(
    scaled_sample_two[uncens_indices_two],
    sd_multiplier
  )

  params <- list(
    mix_weight_1 = continuous_fraction / 2,
    shape1_1 = pars_one$shape1,
    shape2_1 = pars_one$shape2,
    mix_weight_2 = continuous_fraction / 2,
    shape1_2 = pars_two$shape1,
    shape2_2 = pars_two$shape2,
    continuous_fraction = continuous_fraction
  )

  return(params)
}

make_beta_params <- function(sample, sd_multiplier) {
  samp_mean <- mean(sample)
  samp_var <- var(sample) * (sd_multiplier ^ 2)
  samp_var <- max(samp_var, 1e-6)

  t1 <- (samp_mean) * (1 - samp_mean) / (samp_var)
  shape1 <- samp_mean * (t1 - 1)
  shape2 <- (1 - samp_mean) * (t1 - 1)

  if (is.na(shape1) | shape1 <= 0) {
    futile.logger::flog.info(
      "surv_mixture_importance failed to find appropriate beta shape1.
      defaulting to shape1 = 1. This is hopefully just a function of early
      optimisation numerics"
    )

    shape1 <- 1
  }

  if (is.na(shape2) | shape2 <= 0) {
     futile.logger::flog.info(
      "surv_mixture_importance failed to find appropriate beta shape2.
      defaulting to shape2 = 1. This is hopefully just a function of early
      optimisation numerics"
    )

     shape2 <- 1
  }

  res <- list(shape1 = shape1, shape2 = shape2)
  return(res)
}

sample_surv_mix <- function(n, params, censoring_time) {
  samples_per_component <- sample(
    x = c(1, 2, 3),
    size = n,
    replace = TRUE,
    prob = c(
      params$mix_weight_1,
      params$mix_weight_2,
      1 - params$continuous_fraction
    )
  ) %>%
    tabulate(nbins = 3)

  futile.logger::flog.trace(
    "samples_per_component contents: %s",
    paste(samples_per_component, collapse = ", ")
  )

  points_one <- rbeta(
    n = samples_per_component[1],
    shape1 = params$shape1_1,
    shape2 = params$shape2_1
  ) * censoring_time

  points_two <- rbeta(
    n = samples_per_component[2],
    shape1 = params$shape1_2,
    shape2 = params$shape2_2
  ) * censoring_time

  points_three <- rep(censoring_time, samples_per_component[3])

  res <- c(points_one, points_two, points_three)
  return(res)
}

density_surv_mix <- function(x, params, censoring_time, log_scale = FALSE) {
  n_x <- length(x)
  cens_indices <- which(x == censoring_time)
  uncens_indices <- !(seq_len(n_x) %in% cens_indices)
  uncens_x <- x[uncens_indices]

  log_dens_vals <- matrix(data = NA, nrow = n_x, ncol = 3)
  log_dens_vals[uncens_indices, 1] <- log(params$mix_weight_1) +
    dbeta(
      x = uncens_x / censoring_time,
      shape1 = params$shape1_1,
      shape2 = params$shape2_1,
      log = TRUE
    ) -
    log(censoring_time)

  log_dens_vals[uncens_indices, 2] <- log(params$mix_weight_2) +
    dbeta(
      x = uncens_x / censoring_time,
      shape1 = params$shape1_2,
      shape2 = params$shape2_2,
      log = TRUE
    ) -
    log(censoring_time)

  log_dens_vals[cens_indices, 3] <- log(1 - params$continuous_fraction)

  log_dens <- apply(log_dens_vals, 1, matrixStats::logSumExp, na.rm = TRUE)
  stopifnot(all(!is.na(log_dens)))

  if (log_scale) {
    return(log_dens)
  } else {
    return(exp(log_dens))
  }
}
