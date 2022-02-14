uniform_importance <- function(
  sample_one,
  sample_two,
  n_internal_importance_draws
) {
  lower <- min(c(sample_one, sample_two))
  upper <- max(c(sample_one, sample_two))

  stopifnot(lower < upper)

  diff <- upper - lower
  points <- runif(n = n_internal_importance_draws, min = lower, max = upper)
  weights <- rep(1 / diff, times = n_internal_importance_draws)
  res <- list(points = points, weights = weights)
  return(res)
}

# fit a two component mixture model (of t_5 densities)
# evaluate fitted mixture using khat statistics
student_t_mixture_importance <- function() {
  NULL
}

# same as student_t but with gammas? enforce ordering somehow.
# easiest to do with stan.
gamma_mixture_importance <- function() {
  NULL
}
