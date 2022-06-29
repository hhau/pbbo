#' Get the best value of the hyperparameters from a \code{pbbo} result
#'
#' @param pbbo_res An object produced by \code{\link{pbbo}}.
#' @param pbbo_kappa Numeric: A positive number to linearise the dual objective
#'    functions (should \code{pbbo_res} have multiple objectives). Defaults to
#'    \code{NULL}, and will error if not provided and required.
#'
#' @return Named numeric vector, of the form accepted by
#'    \code{prior_predictive_sampler} in \code{\link{pbbo}}.
#' @export
get_best_lambda <- function(pbbo_res, pbbo_kappa = NULL) {
  n_batches <- length(pbbo_res)
  has_multiple_objectives <- !is.null(pbbo_res[[n_batches]][["pareto.front"]])

  if (has_multiple_objectives) {
    if (is.null(pbbo_kappa)) {
      stop(
        paste(
          "Result has multiple objectives but kappa is NULL.",
          "Specify a kappa value to choose an optimum value for lambda."
        )
      )
    } else {
      all_losses <- lapply(seq_len(n_batches), function(batch) {
        pbbo_res[[batch]][["pareto.front"]] %>%
          dplyr::as_tibble() %>%
          dplyr::mutate(
            batch = batch,
            batch_value_id = 1 : dplyr::n(),
            loss = y_1 + pbbo_kappa * y_2
          )
      }) %>%
        dplyr::bind_rows()

      best_index <- which.min(all_losses$loss)
      best_batch <- dplyr::pull(all_losses[best_index, "batch"])
      best_pf_index <- dplyr::pull(all_losses[best_index, "batch_value_id"])
      best_lambda <- pbbo_res[[best_batch]][["pareto.set"]][[best_pf_index]] %>%
        unlist()
    }
  } else {
    best_batch <- lapply(seq_len(n_batches), function(batch) {
      pbbo_res[[batch]][["y"]]
    }) %>%
      unlist() %>%
      which.min()

    best_lambda <- unlist(pbbo_res[[best_batch]][["x"]])
  }

  return(best_lambda)
}