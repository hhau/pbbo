n_internal_importance_draws <- 100

sample_one <- rgamma(n = 100, shape = 3, rate = 0.2)
sample_two <- rgamma(n = 100, shape = 5, rate = 1.2)

gamma_importance_res <- gamma_mixture_importance(
  sample_one = sample_one,
  sample_two = sample_two,
  n_internal_importance_draws = n_internal_importance_draws
)

test_that("gamma importance result is correct shape with valid entries", {
  expect_length(
    gamma_importance_res$points,
    n_internal_importance_draws
  )

  expect_length(
    gamma_importance_res$weights,
    n_internal_importance_draws
  )

  checkmate::expect_numeric(
    x = gamma_importance_res$points,
    finite = TRUE,
    any.missing = FALSE
  )

  checkmate::expect_numeric(
    x = gamma_importance_res$weights,
    finite = TRUE,
    any.missing = FALSE
  )
})

test_that("gamma importance 'weights' are all less than 1 for wide target", {
  checkmate::expect_numeric(
    x = gamma_importance_res$weights,
    lower = 0, # at the moment these aren't weights, but rather density values
    upper = 1
  )
})