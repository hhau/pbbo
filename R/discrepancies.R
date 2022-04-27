build_discrep <- function(
  target_lcdf,
  target_sampler,
  prior_predictive_sampler,
  internal_discrepancy_f,
  n_internal_prior_draws,
  importance_method,
  importance_args,
  n_internal_importance_draws
) {
  local_importance <- get(
    paste0(importance_method, "_importance"),
    envir = environment(pbbo)
  )

  res <- function(lambda_mlrform) {
    prior_sample <- prior_predictive_sampler(
      n_internal_prior_draws,
      lambda_mlrform
    )

    futile.logger::flog.trace(
      "lamba_mlrform contains: %s",
      data.frame(
        names = names(lambda_mlrform),
        vals = as.numeric(lambda_mlrform)
      ) %>%
        tidyr::unite(col = res, sep = ': ') %>%
        magrittr::extract2('res') %>%
        paste0(collapse = ', ')
    )

    futile.logger::flog.trace(
      "prior_sample contains: %s",
      paste(format(prior_sample, digits = 2), collapse = ', ')
    )

    prior_ecdf <- ecdf(prior_sample)
    target_sample <- target_sampler(n_internal_prior_draws)
    importance_draws <- local_importance(
      sample_one = prior_sample,
      sample_two = target_sample,
      n_internal_importance_draws = n_internal_importance_draws,
      importance_args = importance_args
    )

    internal_result <- internal_discrepancy_f(
      prior_ecdf,
      target_lcdf,
      points = importance_draws$points,
      weights = importance_draws$weights
    )

    return(internal_result)
  }

  return(res)
}

build_discrep_covariate <- function(
  target_lcdf,
  target_sampler,
  prior_predictive_sampler,
  covariate_list,
  internal_discrepancy_f,
  n_internal_prior_draws,
  importance_method,
  importance_args,
  n_internal_importance_draws
) {
  n_covariate_obs <- length(covariate_list)

  discrep_list <- lapply(covariate_list, function(a_covariate_draw) {
    local_target_lcdf <- function(x) target_lcdf(x, a_covariate_draw)
    local_target_sampler <- function(n) target_sampler(n, a_covariate_draw)
    local_pp_sampler <- function(n, l) prior_predictive_sampler(n, l, a_covariate_draw)

    inner_discrep <- build_discrep(
      target_lcdf = local_target_lcdf,
      target_sampler = local_target_sampler,
      prior_predictive_sampler = local_pp_sampler,
      internal_discrepancy_f = internal_discrepancy_f,
      n_internal_prior_draws = round(n_internal_prior_draws / n_covariate_obs),
      importance_method = importance_method,
      importance_args = importance_args,
      n_internal_importance_draws = round(n_internal_importance_draws / n_covariate_obs)
    )
  })

  full_discrep <- function(lambda) {
    inner_res <- lapply(discrep_list, function(sub_discrep) {
      sub_discrep(lambda)
    }) %>%
      unlist()

    # switch to taking the mean of the total discrep, so the logSumExp of the
    # log total discrep.
    total_discrep <- -log(n_covariate_obs) + matrixStats::logSumExp(inner_res)
    return(total_discrep)
  }

  return(full_discrep)
}


# cramer-von mises
log_cvm_discrepancy <- function(ecdf_1, log_cdf_2, points, weights) {
  log_n <- log(length(points))
  vals_1 <- ecdf_1(points)
  log_vals_2 <- log_cdf_2(points)

  futile.logger::flog.trace(
    "points contains: %s",
    paste(format(points, digits = 4), collapse = ', ')
  )

  futile.logger::flog.trace(
    "log_vals_2 contains: %s",
    paste(format(log_vals_2, digits = 4), collapse = ', ')
  )

  stopifnot(
    length(points) > 0,
    all(weights > 0),
    all(log_vals_2 <= 0)
  )

  log_top <- 2 * log(abs(vals_1 - exp(log_vals_2)))
  log_weights <- log(weights)
  log_integrand <- log_top - log_weights
  max_c <- max(log_integrand)
  log_res <- log(sum(exp(log_integrand - max_c))) + max_c
  final_res <- -log_n + log_res
  return(as.numeric(final_res))
}

# andersen darling
# note that the second argument is the **log_cdf**
log_ad_discrepancy <- function(ecdf_1, log_cdf_2, points, weights) {
  log_n <- log(length(points))
  high_prec_points <- Rmpfr::mpfr(points, 120)
  high_prec_weights <- Rmpfr::mpfr(weights, 120)

  vals_1 <- ecdf_1(as.numeric(points))
  log_vals_2 <- log_cdf_2(high_prec_points)

  # method assumes that log_cdf_2 returns lcdf values
  # leq to deal with rounding up to zero? for the right hand tail????
  # high prec is supposed to get around that.
  futile.logger::flog.trace(
    "log_vals_2 contains: %s",
    paste(format(log_vals_2, digits = 4), collapse = ', ')
  )

  stopifnot(all(log_vals_2 <= 0))

  # numerator term -- will underflow, which will make me return 0
  # use high prec log_vals_2
  log_top <- 2 * log(abs(vals_1 - exp(log_vals_2)))

  # denominator term
  log_bot <- (log_vals_2) + Rmpfr::log1mexp(a = -log_vals_2)
  log_weights_term <- log(high_prec_weights)
  log_integrand <- log_top - log_bot - log_weights_term
  max_c <- max(log_integrand)
  log_res <- log(sum(exp(log_integrand - max_c))) + max_c
  final_res <- -log_n + log_res
  return(as.numeric(final_res))
}
