# This might be a poor idea? As it effectively truncates the integral to
# the values spanned by [lower, upper]. Empirically this works fine, just
# inefficient (high variance) when samples disagree a lot.
uniform_importance <- function(
  sample_one,
  sample_two,
  n_internal_importance_draws,
  importance_args,
  ...
) {
  lower <- importance_args$uniform_lower
  upper <- importance_args$uniform_upper

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
