n_internal_importance_draws <- 100

outer_cens_time <- 3

sample_one <- c(
  rep(outer_cens_time, 10),
  rbeta(n = 300, 3, 2) * outer_cens_time
)

sample_two <- c(
  rep(outer_cens_time, 58),
  rbeta(310 - 58, 5, 10) * outer_cens_time
)

importance_args = list(
  surv_mixture_sd_multiplier = 1.05,
  surv_mixture_cont_frac = 0.95,
  censoring_time = outer_cens_time
)

surv_mixture_importance_res <- surv_mixture_importance(
  sample_one,
  sample_two,
  n_internal_importance_draws,
  importance_args
)

test_that("surv mix importance correct shape/valid", {
  expect_length(
    surv_mixture_importance_res$points,
    n_internal_importance_draws
  )

  expect_length(
    surv_mixture_importance_res$weights,
    n_internal_importance_draws
  )

  if (requireNamespace("checkmate", quietly = TRUE)) {
    checkmate::expect_numeric(
      x = surv_mixture_importance_res$points,
      finite = TRUE,
      any.missing = FALSE
    )

    checkmate::expect_numeric(
      x = surv_mixture_importance_res$weights,
      finite = TRUE,
      any.missing = FALSE,
      lower = 0
    )
  }
})
