build_discrep <- function(
  target_cdf,
  target_sampler,
  prior_predictive_sampler,
  internal_discrepancy_f,
  n_internal_prior_draws,
  importance_method,
  n_internal_importance_draws
) {
  # will return a function f(lambda) that evaluates the discrepancy
  # given a certain value of lambda.

  # this may need to be inside the function for dispatch to work.
  local_importance <- match.fun(paste0(importance_method, "_importance"))

  res <- function(lambda_mlrform) {
    prior_sample <- prior_predictive_sampler(
      n_internal_prior_draws,
      lambda_mlrform
    )

    prior_ecdf <- ecdf(prior_sample)

    # ideally one could use the samples from the other prior to figure out the
    # range of the target prior, but that seems impossible in general,
    # so for the momemt we just require a target_sampler.
    target_sample <- target_sampler(n_internal_prior_draws)

    importance_draws <- local_importance(
      prior_sample,
      target_sample,
      n_internal_importance_draws
    )

    internal_result <- internal_discrepancy_f(
      cdf_1 = prior_ecdf,
      cdf_2 = target_cdf,
      points = importance_draws$points,
      weights = importance_draws$weights
    )

    return(internal_result)
  }

  return(res)
}

# cramer-von mises
cvm_discrepancy <- function(cdf_1, cdf_2, points, weights) {
  base <- (cdf_1(points) - cdf_2(points))^2
  res <- sum(base / weights)
  return(res)
}

# andersen darling
# the issue is in the denominator, and does not go away with a trivial move
# to the log scale. Computing this accurately will involve:
# 1. computing the lcdf of both samples accurately (in the tails)
# 2. computing the log numerator accurately (log-sum-exp or equiv)
# 3. computing log1m(cdf_2(x)) accurately (look at stan code?)
ad_discrepancy <- function(cdf_1, cdf_2, points, weights) {
  stop(
    "Computing the Anderson–Darling distance is numerically unstable for
    CDFs too far apart."
  )

  vals_1 <- cdf_1(points)
  vals_2 <- cdf_2(points)

  base <- (vals_1 - vals_2)^2 / (vals_2 * (1 - vals_2))
  res <- sum(base / weights)
  return(res)
}

# other functions from @drovandi_comparions_2021? MMD+Laplace, Wasserstein
# these do not require the importance sampling approach, as they are entirely
# sample specific.
# in fact, we could alternate depending on which of target_cdf or target_sampler
# was provided
