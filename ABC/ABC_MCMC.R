#' Approximate Bayesian computation MCMC
#'
#' Apply Metropolis-Hastings MCMC (MH-MCMC) algorithm using approximate Bayesian computation (ABC).
#'
#' @param tol A positive numeric value for the tolerance.
#' @param iter Number of iterations.
#' @param obs A vector of the observed statistics.
#' @param kernel_func A kernel function.
#' @param theta_sigma A covariance matrix for the kernel function.
#' @param init_theta A vector of the initial parameter sampled from the prior.
#' @param prior_func A density function of prior (log density).
#' @param sample_func A function which takes theta and M and return sample mean and variance.
#' @param q_sigma A scale matrix for the Gaussian proposal.
#' @param acc_rate Default acc_rate = FALSE, if TRUE, print acceptance rate and return it.
#' @return A sequence of parameters from the ABC posterior.
ABC_MCMC <- function(tol, iter, obs, kernel_func, theta_sigma, init_theta,
                     prior_func, sample_func, q_sigma,
                     acc_rate=FALSE) {
  # Initial setup
  n_theta <- length(init_theta)
  n_obs <- length(obs)
  theta_matrix <- matrix(NA, nrow=n_theta, ncol=iter)
  if (acc_rate) {accept_num <- 0}
  i <- 1

  # Sample at i = 1
  theta_old <- init_theta
  stats_old <- sample_func(theta_old)
  k_old <- kernel_func(obs, stats_old, tol, theta_sigma)
  theta_matrix[, i] <- theta_old

  # M-H MCMC
  for (i in 2:iter) {
    theta_new <- theta_old + as.vector(rmvnorm(n=1, sigma=q_sigma))
    stats_new <- sample_func(theta_new)
    k_new <- kernel_func(obs, stats_new, tol, theta_sigma)

    log_alpha <- k_new + prior_func(theta_new) - k_old - prior_func(theta_old)
    log_alpha <- min(0, log_alpha)
    log_u <- log(runif(1))

    if (log_u < log_alpha) {
      theta_matrix[, i] <- theta_new
      theta_old <- theta_new
      stats_old <- stats_new
      k_old <- k_new
      if (acc_rate) {accept_num <- accept_num + 1}
    } else {
      theta_matrix[, i] <- theta_old
    }
  }

  result_list <- list(theta=theta_matrix)
  if (acc_rate) {
    print(paste0("Acceptance rate: ", accept_num/iter))
    result_list$acc_rate = accept_num/iter
  }
  return(result_list)
}
