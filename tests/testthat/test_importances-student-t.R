n_internal_importance_draws <- 100

sample_one <- rnorm(n = 300, mean = 10, sd = 5)
sample_two <- 20 + rt(n = 300, df = 8) * 8

student_t_importance_res <- student_t_mixture_importance(
  sample_one = sample_one,
  sample_two = sample_two,
  n_internal_importance_draws = n_internal_importance_draws
)

test_that("student_t importance result is correct shape with valid entries", {
  expect_length(
    student_t_importance_res$points,
    n_internal_importance_draws
  )

  expect_length(
    student_t_importance_res$weights,
    n_internal_importance_draws
  )

  checkmate::expect_numeric(
    x = student_t_importance_res$points,
    finite = TRUE,
    any.missing = FALSE
  )

  checkmate::expect_numeric(
    x = student_t_importance_res$weights,
    finite = TRUE,
    any.missing = FALSE
  )
})

test_that("student_t importance 'weights' are all less than 1 for wide target", {
  checkmate::expect_numeric(
    x = student_t_importance_res$weights,
    lower = 0, # at the moment these aren't weights, but rather density values
    upper = 1
  )
})