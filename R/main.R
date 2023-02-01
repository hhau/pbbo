#' \code{pbbo}: prior by Bayesian optimisation
#' @encoding UTF-8
#'
#' @param model_name String: name of the model/optimisation function
#' @param target_lcdf Function: Vectorised target log-CDF. Takes one argument, a
#'   vector of points at which we wish to evaluate the log-CDF. This function
#'   must accept \code{\link[Rmpfr]{mpfr}} points. This requirement can be
#'   trivially handled by calling \code{as.numeric} on the vector of points, but
#'   performance is improved by using high precision versions of functions
#'   compatible with \code{\link[Rmpfr]{mpfr}} points. See
#'   \code{\link[Rmpfr]{mpfr-class}} for a range of such functions.
#' @param target_lpdf Function: Vectorised target log-PDF (density/mass function
#'   ). Corresponds to \code{target_lcdf}, and is currently required by the
#'   importance sampling algorithm.
#' @param target_sampler Function: Generates samples from the target
#'   distribution. Takes one argument \code{n}, the desired number of samples.
#' @param prior_predictive_sampler Function: Generates samples from the prior
#'   predictive distribution. Takes two arguments: \code{n}, the desired number
#'   of samples; and \code{lambda}, a named vector whose names correspond to the
#'   hyperparameters of interest (and defined by \code{param_set}).
#' @param param_set A \code{\link[ParamHelpers]{makeParamSet}}: Parameter set
#'   corresponding to the hyperparameters of interest. Best understood by
#'   inspecting the example. For an example with vectors of parameters, we could
#'   take
#'
#'   \preformatted{param_set <- makeParamSet(
#'     makeNumericVectorParam(id = 'mu', len = 2, lower = -50, upper = 50),
#'     makeNumericParam(id = 'sigma', lower = 0, upper = 20)
#' )}
#'
#'   the elements of which would be accessible inside
#'   \code{prior_predictive_sampler} as \code{c(lambda['mu1'], lambda['mu2'],
#'   lambda['sigma'])}.
#' @param covariate_values An optional vector or matrix/data.frame of covariate
#'   values: \bold{including covariate values changes the signatures required
#'   for} \code{target_lcdf}, \code{target_sampler}, \bold{and}
#'   \code{prior_predictive_sampler}. Each of these function must take an
#'   additional third argument, perhaps called \code{cov}, that dictates the
#'   behaviour of the function at a specific covariate value, i.e. value in the
#'   covariate vector or row in the covariate matrix.
#' @param discrepancy String (or Function): One of either \code{'log_cvm'} or
#'   \code{'log_ad'} for Cramér–von Mises (default) or Anderson–Darling distance
#'   respectively. The latter is currently unsupported as it is too numerically
#'   difficult to compute accurately. Alternatively, one can supply a function
#'   with the signature \code{function(ecdf_1, lcdf_2, points, weights)}, where
#'   the first argument is the ECDF of the PPD, the second is the LCDF of the
#'   target, \code{points} is a vector of evaluation points, and \code{weights}
#'   is the corresponding vector of importance weights.
#' @param n_crs2_iters Numeric: Number of iterations to run the CRS2 algorithm
#'   from \code{\link[nloptr]{nloptr}}. These are used to initialise the Bayes
#'   Opt stages, as it can be better at finding global minima. It only handles
#'   single objectives, so it does not include the \code{extra_objective_term}.
#' @param n_internal_prior_draws Numeric: Number of draws to generate from the
#'   prior predictive distribution for the given value of \code{lambda} in each
#'   of the optimisation iterations. More draws result in better estimates of
#'   the ECDF. Particularly importance if \code{discrepancy = "log_ad"}.
#' @param importance_method String: Defaults to \code{'uniform'}. Specifies
#'   which importance sampling method should be used. Options include \code{
#'   gamma_mixture} for data on the positive real numbers and \code{
#'   student_t_mixture} for data on the whole real line.
#' @param importance_args List: Additional optional arguments to the importance
#'   functions. These should be in the form of \code{importance_param}, e.g.
#'   \code{uniform_lower} to set the lower bound for the uniform importance
#'   method. Possible arguments include:
#'   \itemize{
#'   \item{\code{uniform_lower}}
#'   {Left as NULL (the default) the lower bound of the uniform importance
#'   distribution is set to be \code{min(sample_from_target ,
#'   sample_from_current_prior_predictive)}. If you wish to set this lower bound
#'   manually, perhaps if the draws from the prior predictive distribution have
#'   a compact support, then set \code{uniform_lower} to a specific value}
#'   \item{\code{uniform_upper}}
#'   {Corresponding upper bound to \code{
#'   uniform_lower}} \item{\code{gamma_sd_multiplier}}{Constant with which to
#'   multiply the estimated standard deviations for each of the mixture
#'   components. Can help in ensuring the tails are sufficiently heavy.}
#'   \item{\code{student_t_multiplier}}
#'   {Same purpose as \code{
#'   gamma_sd_multiplier} but for the mixture of student-t densities
#'   importance.}
#'   }
#' @param n_internal_importance_draws Numeric: Number of draws to generate from
#'   the importance distribution, and subsequently used to evaluate the
#'   discrepancy integral.
#' @param n_design_pad Numeric: Number of additional samples added to designs
#'   for numerical stability/robustness. These are drawn from a Latin hypercube.
#' @param bayes_opt_batches Numeric: Number of batches to run the Bayesian
#'   optimisation algorithm for. Minimum 1.
#' @param bayes_opt_iters_per_batch Numeric: Number of iterations of Bayesian
#'   optimisation to use per batch. Passed to
#'   \code{\link[mlrMBO]{setMBOControlTermination}} as \code{iters}.
#' @param bayes_opt_design_points_per_batch Numeric: Number of points from
#'   previous batch to carry forward into the next batch. Points are selected by
#'   sampling, with weights proportional to the value of the objective function
#'   at each point.
#' @param bayes_opt_print Boolean: if \code{TRUE}, print the progress/status
#'   from \code{\link[mlrMBO]{mbo}}. Defaults to \code{FALSE}.
#' @param extra_objective_term Function: univariate real-value function that
#'   takes \code{lambda} as an argument and the result is added on to the
#'   objective function.
#' @param ... Currently unused.
#'
#' @return A list of length \code{bayes_opt_batches}, with each element being a
#'   \code{\link[mlrMBO]{MBOSingleObjResult}}.
#' @export
#'
#' @examples
#' \dontrun{
#' library(pbbo)
#' library(mlrMBO)
#'
#' target_lcdf <- function(x) {
#'   pnorm(x, mean = 2, sd = 0.5, log.p = TRUE)
#' }
#'
#' target_sampler <- function(n) {
#'   rnorm(n = n, mean = 2, sd = 0.5)
#' }
#'
#' prior_predictive_sampler <- function(n, lambda) {
#'   rnorm(n = n, mean = lambda['mu'], sd = lambda['sigma'])
#' }
#'
#' param_set <- makeParamSet(
#'   makeNumericParam(id = 'mu', default = 0.2, lower = -50, upper = 50),
#'   makeNumericParam(id = 'sigma', lower = 0, upper = 20, default = 0.2)
#' )
#'
#' pbbo(
#'   model_name = 'test_normal',
#'   target_lcdf = target_lcdf,
#'   target_sampler = target_sampler,
#'   prior_predictive_sampler = prior_predictive_sampler,
#'   discrepancy = 'cvm',
#'   param_set = param_set,
#'   n_internal_prior_draws = 750,
#'   n_internal_importance_draws = 200,
#'   bayes_opt_iters_per_batch = 50,
#'   bayes_opt_print = FALSE
#' )}
pbbo <- function(
  model_name = "default",
  target_lcdf,
  target_lpdf,
  target_sampler,
  prior_predictive_sampler,
  param_set,
  covariate_values = NULL,
  discrepancy = "log_cvm",
  n_crs2_iters = 2000,
  n_internal_prior_draws = 5e4,
  importance_method = "uniform",
  importance_args = list(
    uniform_lower = NULL,
    uniform_upper = NULL,
    gamma_sd_multiplier = 1.05,
    student_t_sd_multiplier = 1.05,
    surv_mixture_sd_multiplier = 1.05,
    surv_mixture_cont_frac = 0.95,
    surv_mixture_cens_times = rep(NULL, length(covariate_values))
  ),
  n_internal_importance_draws = 5e3,
  n_design_pad = 10,
  bayes_opt_batches = 1,
  bayes_opt_iters_per_batch = 300,
  bayes_opt_design_points_per_batch = min(
    4 * length(param_set$pars),
    bayes_opt_iters_per_batch
  ),
  bayes_opt_print = FALSE,
  extra_objective_term = NULL,
  ...
) {
  stopifnot(
    is.function(target_lcdf),
    is.function(target_lpdf),
    is.function(target_sampler),
    is.function(prior_predictive_sampler),
    is.list(param_set),
    discrepancy %in% c("log_ad", "log_cvm") | is.function(discrepancy),
    is.numeric(n_internal_prior_draws),
    is.numeric(n_internal_importance_draws),
    1 < n_internal_prior_draws,
    importance_method %in% c(
      "uniform", "gamma_mixture", "student_t_mixture", "surv_mixture"
    ),
    is.list(importance_args),
    is.numeric(bayes_opt_iters_per_batch),
    1 < bayes_opt_iters_per_batch,
    1 < n_internal_importance_draws,
    1 <= bayes_opt_batches
  )

  if (is.function(discrepancy)) {
    internal_discrepancy_f <- discrepancy
  } else {
    x <- paste0(discrepancy, "_discrepancy")
    internal_discrepancy_f <- get(x, envir = environment(pbbo))
  }

  if (!is.null(covariate_values)) {
    if (is.vector(covariate_values)) {
      covariate_list <- as.list(covariate_values)
      covariate_entry_lengths <- sapply(covariate_list, length)
      stopifnot(
        length(covariate_list) == length(covariate_values),
        all(covariate_entry_lengths == covariate_entry_lengths[1])
      )
    } else if (is.matrix(covariate_values) | is.data.frame(covariate_values)) {
      covariate_list <- matrixlike_to_rowlist(covariate_values)
      stopifnot(length(covariate_list) == nrow(covariate_values))
    } else {
      stop("Unsure how to handle non vector/matrix-like covariate values.")
    }

    discrep_partial <- build_discrep_covariate(
      target_lcdf = target_lcdf,
      target_lpdf = target_lpdf,
      target_sampler = target_sampler,
      prior_predictive_sampler = prior_predictive_sampler,
      covariate_list = covariate_list,
      internal_discrepancy_f = internal_discrepancy_f,
      n_internal_prior_draws = n_internal_prior_draws,
      importance_method = importance_method,
      importance_args = importance_args,
      n_internal_importance_draws = n_internal_importance_draws
    )
  } else {
    discrep_partial <- build_discrep(
      target_lcdf = target_lcdf,
      target_lpdf = target_lpdf,
      target_sampler = target_sampler,
      prior_predictive_sampler = prior_predictive_sampler,
      internal_discrepancy_f = internal_discrepancy_f,
      n_internal_prior_draws = n_internal_prior_draws,
      importance_method = importance_method,
      importance_args = importance_args,
      n_internal_importance_draws = n_internal_importance_draws
    )
  }

  initial_design <- design_from_crs2(
    discrep = discrep_partial,
    param_set = param_set,
    n_design = bayes_opt_design_points_per_batch,
    n_design_pad = n_design_pad,
    n_crs2_iters = n_crs2_iters
  )

  if (!is.null(extra_objective_term) & is.function(extra_objective_term)) {
    objective_function <- smoof::makeMultiObjectiveFunction(
      name = model_name,
      fn = function(l) c(discrep_partial(l), extra_objective_term(l)),
      n.objectives = 2,
      par.set = param_set,
      minimize = c(TRUE, TRUE),
      noisy = FALSE
    )

    control_obj <- mlrMBO::makeMBOControl(
      n.objectives = 2,
      #final.method = 'best.predicted',
      propose.points = 1,
      #final.evals = 5
    ) %>%
      mlrMBO::setMBOControlInfill(opt = "nsga2") %>%
      mlrMBO::setMBOControlTermination(iters = bayes_opt_iters_per_batch) %>%
      mlrMBO::setMBOControlMultiObj(method = "mspot")

  } else {
    objective_function <- smoof::makeSingleObjectiveFunction(
      name = model_name,
      fn = discrep_partial,
      par.set = param_set,
      minimize = TRUE,
      noisy = TRUE
    )

    control_obj <- mlrMBO::makeMBOControl(final.method = "best.predicted") %>%
      mlrMBO::setMBOControlTermination(iters = bayes_opt_iters_per_batch)
  }

  full_batch_res <- list()
  for (batch_num in seq_len(bayes_opt_batches)) {
    futile.logger::flog.info("Starting Bayes opt batch %d", batch_num)
    if (batch_num == 1) {
      batch_design <- initial_design

      if (!is.null(extra_objective_term)) {
        batch_design <- subset(batch_design, select = -y)
      }

    } else {
      prev_batch_path <- full_batch_res[[batch_num - 1]][["opt.path"]] %>%
        as.data.frame()

      all_names <- colnames(prev_batch_path)
      unneeded_cols <- c(
        "dob", "eol", "error.message", "exec.time", "cb", "cb.y_1", "cb.y_2",
        "error.model", "train.time", "prop.type", "propose.time", "se", "mean",
        "lambda", "hv.contr"
      )

      design_col_names <- setdiff(all_names, unneeded_cols)

      if (!is.null(extra_objective_term)) {
        log_raw_weights <- -exp(prev_batch_path$y_1)

        prev_batch_pf_design <- full_batch_res[[batch_num - 1]][["pareto.set"]] %>%
          lapply(function(x) as.data.frame(x)) %>%
          dplyr::bind_rows()

        bayes_opt_design_points_per_batch <- max(
          1,
          bayes_opt_design_points_per_batch - nrow(prev_batch_pf_design)
        )
      } else {
        log_raw_weights <- -exp(prev_batch_path$y)
      }

      stopifnot(any(!is.na(log_raw_weights)))
      weights <- exp(log_raw_weights - matrixStats::logSumExp(log_raw_weights))
      prev_batch_indices <- sample(
        x = seq_len(nrow(prev_batch_path)),
        size = bayes_opt_design_points_per_batch,
        replace = FALSE,
        prob = weights
      )

      batch_design <- prev_batch_path[prev_batch_indices, design_col_names]
      pad_design <- ParamHelpers::generateDesign(
        n = n_design_pad,
        par.set = param_set,
        fun = lhs::maximinLHS
      )

      batch_design <- dplyr::bind_rows(batch_design, pad_design)

      if (!is.null(extra_objective_term)) {
        batch_design <- subset(batch_design, select = -c(y_1, y_2))
        batch_design <- dplyr::bind_rows(prev_batch_pf_design, batch_design)
      } else {
        batch_design <- subset(batch_design, select = -y)
      }

    }

    batch_mbo_res <- mlrMBO::mbo(
      fun = objective_function,
      design = batch_design,
      control = control_obj,
      show.info = bayes_opt_print
    )

    full_batch_res[[batch_num]] <- batch_mbo_res
  }

  return(full_batch_res)
}
