#' Approximate Bayesian computation MCMC
#'
#' Apply Metropolis-Hastings MCMC (MH-MCMC) algorithm using approximate Bayesian computation (ABC).
#'
#' @param tol A positive numberic value for the tolerance.
#' @param iter Number of iterations.
#' @param obs A vector of the observed statistics.
#' @param kernel_func A kernel function.
#' @param init_theta A vector of the initial parameter sampled from the prior.
#' @param prior_func A density function of prior (log density).
#' @param sample_func A function which takes theta and M and return sample mean and variance.
#' @param q_sigma A scale matrix for the Gaussian proposal.
#' @param acc_rate Default acc_rate = FALSE, if TRUE, print acceptance rate and return it.
#' @return A sequence of parameters from the BSL posterior.
#'
#' @param obs The vector of observed data point.
#' @param tol A positive numeric value for the tolerance.
#' @param kernel_func A kernel function.
#' @param p_theta A function to provide sampled theta from proposal.
#' @param d_theta A function to provide the log density of proposal.
#' @param p_s A function to provide the sampled summary statistics by given parameters.
#' @param prior A function to provide the log density of prior.
#' @param theta_0 Initial theta from the prior.
#' @param n_iter Number of iterations.
#' @param sigma The covariance matrix.
#' @return Function value.
#' @export
ABC_MCMC <- function(tol, iter, obs, kernel_func, init_theta, prior_func,
                     sample_func, q_sigma, acc_rate=FALSE) {
  theta_matrix <- matrix(NA, nrow=(n_iter+1), ncol=length(theta_0))
  s_matrix <- matrix(NA, nrow=(n_iter+1), ncol=length(obs))
  accept_vec <- rep(FALSE, n_iter+1)

  s_0 <- as.vector(p_s(theta_0))
  k_0 <- kernel_func(obs, s_0, tol, sigma)

  theta_matrix[1, ] <- theta_0
  s_matrix[1, ] <- s_0

  for (i in 2:(n_iter+1)) {
    theta_1 <- p_theta(theta_0)
    s_1 <- as.vector(p_s(theta_1))
    k_1 <- kernel_func(obs, s_1, tol, sigma)

    log_alpha <- k_1+prior(theta_1)+d_theta(theta_0, theta_1)-
      k_0-prior(theta_0)-d_theta(theta_1, theta_0)
    if (log(runif(1)) < log_alpha) {
      theta_0 <- theta_1
      s_0 <- s_1
      accept_vec[i] <- TRUE
    }

    theta_matrix[i, ] <- theta_0
    s_matrix[i, ] <- s_0
  }

  return(list(theta_matrix=theta_matrix,
              s_matrix=s_matrix,
              accept_vec=accept_vec))
}
