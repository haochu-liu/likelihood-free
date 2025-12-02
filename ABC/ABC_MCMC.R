#' ABC-MCMC
#'
#' Run ABC-MCMC by given proposal and prior.
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
abc_mcmc <- function(obs, tol, kernel_func, p_theta, d_theta, p_s, prior,
                     theta_0, n_iter, sigma=NULL) {
  if (is.null(sigma)) {sigma=diag(rep(1, length(obs)))}

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
