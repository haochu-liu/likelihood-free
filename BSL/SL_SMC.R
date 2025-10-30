#' Synthetic Likelihood SMC
#'
#' Apply sequential Monte Carlo (SMC) algorithm for BSL posterior as the target distribution.
#'
#' @param M Number of new data points drawn in each iteration.
#' @param alpha Number to control the effective sample size in reweight.
#' @param N Number of particles for SMC.
#' @param theta_d Dimension of parameter.
#' @param obs A vector of the observed statistics.
#' @param prior_func A function to sample parameters from prior.
#' @param sample_func A function which takes theta and M and return sample mean and variance.
#' @param sigma A scale parameter for the Gaussian proposal in move step.
#' @param theta_history Default theta_history = FALSE, if TRUE, return all particles in history.
#' @return A vector of parameters from the BSL posterior.
SL_MCMC <- function(M, alpha, N, theta_d, obs, prior_func, sample_func, sigma,
                    theta_history=FALSE) {
  theta_mat <- matrix(NA, nrow=theta_d, ncol=N)
  if (theta_history) {
    iter_max <- 100
    theta_history_mat <- matrix(NA, nrow=(iter_max*theta_d), ncol=N)
  }
  mu_mat <- matrix(NA, nrow=length(obs), ncol=N)
  sigma_array <- array(data = NA, dim = c(length(obs), length(obs), N))
  weight_vec <- rep(-log(N), N)
  ess_flat <- ESS_weight(weight_vec)
  gamma_old <- 0

  # Initialization
  for (n in 1:N) {
    theta_n <- prior_func()
    sample_sta <- sample_func(theta_n, M)

    theta_mat[, n] <- theta_n
    mu_mat[, n] <- sample_sta$mean
    sigma_array[, , n] <- sample_sta$sigma
  }

  while (gamma_old < 1) {
    # Binary search 100 times
    gamma_new <- (1 + gamma_old) / 2
    for (i in 1:100) {
      # Reweight
      for (n in 1:N) {
        weight_vec[n] <- dmvnorm(x=obs,
                                 mean=mu_mat[, n],
                                 sigma=as.matrix(sigma_array[, , n]),
                                 log=TRUE) * (gamma_new - gamma_old)
        weight_vec <- weight_vec - logSumExp(weight_vec)
      }

      ess_new <- ESS_weight(weight_vec)

      if (ess_new < (log(alpha)+ess_flat)) {
        # ESS too small
        gamma_new <- (gamma_new + gamma_old) / 2
      } else {
        # ESS too large
        gamma_new <- (gamma_new + 1) / 2
        # Try gamma_new = 1 once
        if (i == 1) {
          for (n in 1:N) {
            weight_vec[n] <- dmvnorm(x=obs,
                                     mean=mu_mat[, n],
                                     sigma=as.matrix(sigma_array[, , n]),
                                     log=TRUE) * (1 - gamma_old)
            weight_vec <- weight_vec - logSumExp(weight_vec)
          }
          ess_new <- ESS_weight(weight_vec)
          if (ess_new >= (log(alpha)+ess_flat)) {
            gamma_new <- 1
            break
          }
        }
      }

      if (abs(ess_new - (log(alpha)+ess_flat)) < 0.01) {
        break
      }
    }


    # Resample
    prob_vec <- exp(weight_vec)
    resample_index <- sample(1:N, size=N, replace=TRUE, prob=prob_vec)
    theta_mat <- theta_mat[, resample_index]

  }
}


