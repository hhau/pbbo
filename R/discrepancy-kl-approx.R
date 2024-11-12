build_approx_kl_discrep <- function(
  target_sampler = target_sampler,
  prior_predictive_sampler = prior_predictive_sampler,
  covariate_list = covariate_list,
  n_samples_for_approx = n_internal_prior_draws
) {

  out <- function(lambda) {
    # we could, ofc, fix this to minimise computational
    # effort - later
    target_samples <- matrix(
      data = NA,
      nrow = n_samples_for_approx,
      ncol = length(covariate_list)
    )

    prior_pred_samples <- matrix(
      data = NA,
      nrow = n_samples_for_approx,
      ncol = length(covariate_list)
    )

    # work through the covariate list and draw
    for (ii in seq_len(length(covariate_list))) {
      target_samples[, ii] <- target_sampler(
        n_samples_for_approx,
        covariate_list[[ii]]
      )

      prior_pred_samples[, ii] <- prior_predictive_sampler(
        n_samples_for_approx,
        lambda,
        covariate_list[[ii]]
      )
    }

    # compute the mean/covariance
    target_mean <- apply(target_samples, 2, mean)
    target_cov <- cov(target_samples)
    prior_pred_mean <- apply(prior_pred_samples, 2, mean)
    prior_pred_cov <- cov(prior_pred_samples)

    # compute the KL
    kl_val <- kl_divergence_gaussian(
      mu1 = target_mean,
      sigma1 = target_cov,
      mu2 = prior_pred_mean,
      sigma2 = prior_pred_cov
    )

    return(kl_val)
  }

  # return a function of lambda (l) that computes the
  # KL approx discrep
  return(out)
}

kl_divergence_gaussian <- function(mu1, sigma1, mu2, sigma2) {
  stopifnot(
    length(mu1) == length(mu2),
    nrow(sigma1) == ncol(sigma1),
    nrow(sigma1) == length(mu1),
    nrow(sigma2) == ncol(sigma2),
    nrow(sigma2) == length(mu2)
  )

  # Compute the KL divergence
  d <- length(mu1)
  trace_term <- sum(diag(solve(sigma2, sigma1)))
  det_term <- log(det(sigma2)) - log(det(sigma1))
  mean_term <- t(mu1 - mu2) %*% solve(sigma2, mu1 - mu2)

  kl_div <- 0.5 * (trace_term + det_term + mean_term - d)
  return(log(as.numeric(kl_div)))
}
