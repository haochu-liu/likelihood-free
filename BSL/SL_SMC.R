#' Synthetic Likelihood SMC
#'
#' Apply sequential Monte Carlo (SMC) algorithm for BSL posterior as the target distribution.
#'
#' @param M Number of new data points drawn in each iteration.
#' @param alpha Number to control the effective sample size in reweight.
#' @param N Number of particles for SMC.
#' @param theta_d Dimension of parameter.
#' @param obs A vector of the observed statistics.
#' @param prior_sampler A function to draw samples form prior.
#' @param prior_func A density function of prior (log density).
#' @param sample_func A function which takes theta and M and return sample mean and variance.
#' @param sigma A scale parameter for the Gaussian proposal in move step.
#' @param theta_history Default theta_history = FALSE, if TRUE, return all particles in history.
#' @param gamma_history Default gamma_history = FALSE, if TRUE, return gamma history.
#' @return A vector of parameters from the BSL posterior.
SL_SMC <- function(M, alpha, N, theta_d, obs, prior_sampler, prior_func,
                   sample_func, sigma,
                   theta_history=FALSE, gamma_history=FALSE) {
  theta_mat <- matrix(NA, nrow=theta_d, ncol=N)
  iter_max <- 50
  if (theta_history) {
    theta_history_mat <- matrix(NA, nrow=(iter_max*theta_d), ncol=N)
  }
  mu_mat <- matrix(NA, nrow=length(obs), ncol=N)
  sigma_array <- array(data = NA, dim = c(length(obs), length(obs), N))
  weight_vec <- rep(-log(N), N)
  ess_flat <- ESS_weight(weight_vec)
  q_sigma <- sigma * diag(theta_d)
  gamma_old <- 0
  iter <- 1
  if (gamma_history) {
    gamma_vec <- c(gamma_old)
  }

  # Initialization
  for (n in 1:N) {
    theta_n <- prior_sampler()
    sample_sta <- sample_func(theta_n, M)

    theta_mat[, n] <- theta_n
    mu_mat[, n] <- sample_sta$mean
    sigma_array[, , n] <- sample_sta$sigma
  }

  while (gamma_old < 1) {
    # Log-likelihood for current parameters
    log_likelihood <- rep(NA, N)
    for (n in 1:N) {
      log_likelihood[n] <- dmvnorm(x=obs,
                                   mean=mu_mat[, n],
                                   sigma=as.matrix(sigma_array[, , n]),
                                   log=TRUE)
    }

    # Binary search 100 times
    gamma_new <- (1 + gamma_old) / 2
    # Reweight
    weight_vec <- (gamma_new - gamma_old) * log_likelihood
    weight_vec <- weight_vec - logSumExp(weight_vec)
    ess_new <- ESS_weight(weight_vec)
    for (i in 1:100) {
      if (abs(ess_new - (log(alpha)+ess_flat)) < 0.01) {
        break
      }

      if (ess_new < (log(alpha)+ess_flat)) {
        # ESS too small
        gamma_new <- (gamma_new + gamma_old) / 2
      } else {
        # ESS too large
        gamma_new <- (gamma_new + 1) / 2
        # Try gamma_new = 1 once
        if (i == 1) {
          weight_vec <- (gamma_new - gamma_old) * log_likelihood
          weight_vec <- weight_vec - logSumExp(weight_vec)
          ess_new <- ESS_weight(weight_vec)
          if (ess_new >= (log(alpha)+ess_flat)) {
            gamma_new <- 1
          }
        }
      }

      # Reweight
      weight_vec <- (gamma_new - gamma_old) * log_likelihood
      weight_vec <- weight_vec - logSumExp(weight_vec)
      ess_new <- ESS_weight(weight_vec)

      if (gamma_new == 1) {break}
    }

    # Resample
    prob_vec <- exp(weight_vec)
    resample_index <- sample(1:N, size=N, replace=TRUE, prob=prob_vec)
    theta_mat <- theta_mat[, resample_index, drop=FALSE]
    mu_mat <- mu_mat[, resample_index, drop=FALSE]
    log_likelihood <- log_likelihood[resample_index]
    sigma_array <- sigma_array[, , resample_index, drop=FALSE]
    weight_vec <- rep(-log(N), N)

    # Move
    for (n in 1:N) {
      theta_new <- theta_mat[, n] + as.vector(rmvnorm(n=1, sigma=q_sigma))
      theta_new <- max(0, theta_new)
      stats_new <- sample_func(theta_new, M)
      sl_new <- dmvnorm(x=obs,
                        mean=stats_new$mean, sigma=stats_new$sigma, log=TRUE)


      log_alpha <- sl_new + prior_func(theta_new) -
        log_likelihood[n] - prior_func(theta_mat[, n] )
      log_alpha <- min(0, log_alpha)
      log_u <- log(runif(1))

      if (log_u < log_alpha) {
        theta_mat[, n] <- theta_new
        mu_mat[, n] <- stats_new$mean
        sigma_array[, , n] <- stats_new$sigma
      }
    }

    # Update
    gamma_old <- gamma_new
    if (gamma_history) {
      gamma_vec <- c(gamma_vec, gamma_old)
    }
    if (theta_history) {
      theta_history_mat[(1+(iter-1)*theta_d):(iter*theta_d), ] <- theta_mat
    }
    iter <- iter + 1
    if (iter >= iter_max) {
      break
    }
  }

  if (theta_history) {
    return(theta_history_mat)
  } else if (gamma_history) {
    return(list(theta = theta_mat,
                gamma = gamma_vec))
  } else {
    return(theta_mat)
  }
}
