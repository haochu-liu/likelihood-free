#' Synthetic Likelihood MCMC
#'
#' Apply Metropolis-Hastings MCMC (MH-MCMC) algorithm for BSL posterior as the target distribution.
#'
#' @param M Number of new data points drawn in each iteration.
#' @param iter Number of iterations.
#' @param obs A vector of the observed statistics.
#' @param init_theta A vector of the initial parameter sampled from the prior.
#' @param prior_func A density function of prior (log density).
#' @param sample_func A function which takes theta and M and return sample mean and variance.
#' @param q_sigma A scale matrix for the Gaussian proposal.
#' @param acc_rate Default acc_rate = FALSE, if TRUE, print acceptance rate and return it.
#' @return A sequence of parameters from the BSL posterior.
SL_MCMC <- function(M, iter, obs, init_theta, prior_func, sample_func, q_sigma,
                    acc_rate=FALSE) {
  # Initial setup
  n_theta <- length(init_theta)
  n_obs <- length(obs)
  theta_matrix <- matrix(NA, nrow=n_theta, ncol=iter)
  accept_num <- 0
  i <- 1

  # Sample and likelihood at i = 1
  theta_old <- init_theta
  stats_old <- sample_func(theta_old, M)
  sl_old <- dmvnorm(x=obs, mean=stats_old$mean, sigma=stats_old$sigma, log=TRUE)
  theta_matrix[, i] <- init_theta

  # M-H MCMC
  for (i in 2:iter) {
    theta_new <- theta_old + as.vector(rmvnorm(n=1, sigma=q_sigma))
    stats_new <- sample_func(theta_new, M)
    sl_new <- dmvnorm(x=obs, mean=stats_new$mean, sigma=stats_new$sigma, log=TRUE)

    log_alpha <- sl_new + prior_func(theta_new) - sl_old - prior_func(theta_old)
    log_alpha <- min(0, log_alpha)
    log_u <- log(runif(1))

    if (log_u < log_alpha) {
      theta_matrix[, i] <- theta_new
      theta_old <- theta_new
      stats_old <- stats_new
      sl_old <- sl_new
      accept_num <- accept_num + 1
    } else {
      theta_matrix[, i] <- theta_old
    }
  }

  if (acc_rate) {
    print(paste0("Acceptance rate: ", accept_num/iter))
    return(list(theta = theta_matrix,
                acc_rate = accept_num/iter))
  } else {
    return(theta_matrix)
  }
}
