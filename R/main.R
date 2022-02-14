#' \code{pbbo}: prior by Bayesian optimisation
#' @param model_name String: name of the model/optimisation function
#' @param target_cdf Function: Vectorised target CDF. Takes one argument, a
#'   vector of points at which we wish to evaluate the CDF.
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
#' @param discrepancy String (or Function): One of either \code{'cvm'} or
#'   \code{'ad'} for Cramér–von Mises (default) or Anderson–Darling distance
#'   respectively. The latter is currenly unsupported as it is too numerically
#'   difficult to compute accurately. Alternatively, one can supply a function
#'   with thesignature \code{function(cdf_1, cdf_2, points, weights)}, where the
#'   first two arguments are (E)CDF functions, \code{points} is a vector of
#'   evaluation points, and \code{weights} is the corresponding vector of
#'   importance weights.
#' @param initial_points_to_eval \code{data.frame}: A data.frame of points at which
#'   one might wish to evaluate the discrepancy function at initially. The colums
#'   of this data.farme must match the names as described by \code{param_set}, and
#'   include a column called \code{y}, which may be full of \code{NA}s if the
#'   discrepancy has not been evaluated at the points. Defaults to \code{NULL}.
#' @param n_internal_prior_draws Numeric: Number of draws to generate from the
#'   prior predictive distribution for the given value of \code{lambda} in each
#'   of the optimisation iterations. More draws result in better estimates of the
#'   ECDF.
#' @param importance_method String: Defaults to \code{'uniform'}. Perhaps in the
#'   future this will allow one to adjust the type of importance sampling use to
#'   compute the discrepancy function.
#' @param n_internal_importance_draws Numeric: Number of draws to generate from
#'   the importance distribution, and subsequently used to evalute the discrepancy
#'   intergral.
#' @param bayes_opt_iters Numeric: Number of iterations of Bayesian optimisation
#'   to use. Passed to \code{\link[mlrMBO]{setMBOControlTermination}} as
#'   \code{iters}.
#' @param ... Other options passed to ?. Currently unused.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' library(pbbo)
#' library(mlrMBO)
#'
#' target_cdf <- function(x) {
#'   pnorm(x, mean = 2, sd = 0.5)
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
#'   target_cdf = target_cdf,
#'   target_sampler = target_sampler,
#'   prior_predictive_sampler = prior_predictive_sampler,
#'   discrepancy = 'cvm',
#'   param_set = param_set,
#'   n_internal_prior_draws = 750,
#'   n_internal_importance_draws = 200,
#'   bayes_opt_iters = 50,
#'   bayes_opt_print = FALSE
#' )
#' }
pbbo <- function(
  model_name = 'default',
  target_cdf,
  target_sampler,
  prior_predictive_sampler,
  param_set,
  discrepancy = 'cvm',
  initial_points_to_eval = NULL,
  n_internal_prior_draws = 250,
  importance_method = 'uniform',
  n_internal_importance_draws = 100,
  bayes_opt_iters = 100,
  bayes_opt_print = FALSE,
  ...
) {
  stopifnot(
    is.function(target_cdf),
    is.function(target_sampler),
    is.function(prior_predictive_sampler),
    is.list(param_set),
    discrepancy %in% c('ad', 'cvm') | is.function(discrepancy),
    is.numeric(n_internal_prior_draws),
    is.numeric(n_internal_importance_draws),
    1 < n_internal_prior_draws,
    importance_method %in% c('uniform'),
    1 < n_internal_importance_draws
  )

  if (is.function(discrepancy)) {
    internal_discrepancy_f <- discrepancy
  } else {
    x <- paste0(discrepancy, '_discrepancy')
    internal_discrepancy_f <- match.fun(x)
  }

  discrep <- build_discrep(
    target_cdf = target_cdf,
    target_sampler = target_sampler,
    prior_predictive_sampler = prior_predictive_sampler,
    internal_discrepancy_f = internal_discrepancy_f,
    n_internal_prior_draws = n_internal_prior_draws,
    importance_method = importance_method,
    n_internal_importance_draws = n_internal_importance_draws
  )

  objective_function <- smoof::makeSingleObjectiveFunction(
    name = model_name,
    fn = discrep,
    par.set = param_set,
    minimize = TRUE,
    noisy = TRUE
  )

  control_obj <- mlrMBO::makeMBOControl(final.method = 'best.predicted') %>%
    mlrMBO::setMBOControlTermination(iters = bayes_opt_iters)

  res <- mlrMBO::mbo(
    fun = objective_function,
    design = initial_points_to_eval,
    control = control_obj,
    show.info = bayes_opt_print
  )

  return(res)
}
