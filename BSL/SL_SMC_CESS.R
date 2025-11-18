#' Synthetic Likelihood SMC (CESS)
#'
#' Apply sequential Monte Carlo (SMC) algorithm for BSL posterior as the target distribution.
#' Use CESS rather than ESS to measure the discrepancy between distributions.
#'
#' @param M Number of new data points drawn in each iteration.
#' @param alpha Number to control the effective sample size in reweight.
#' @param N Number of particles for SMC.
#' @param theta_d Dimension of parameter.
#' @param obs A vector of the observed statistics.
#' @param prior_sampler A function to draw samples form prior.
#' @param prior_func A density function of prior (log density).
#' @param sample_func A function which takes theta and M and return sample mean and variance.
#' @param q_sigma A scale matrix for the Gaussian proposal.
#' @param theta_history Default theta_history = FALSE, if TRUE, return all particles in history.
#' @param gamma_history Default gamma_history = FALSE, if TRUE, return gamma history.
#' @return A vector of parameters from the BSL posterior.
SL_SMC_CESS <- function(M, alpha, N, theta_d, obs, prior_sampler, prior_func,
                        sample_func, q_sigma,
                        theta_history=FALSE, gamma_history=FALSE) {
  theta_mat <- matrix(NA, nrow=theta_d, ncol=N)
  iter_max <- 50
  if (theta_history) {
    theta_history_array <- array(data = NA, dim = c(theta_d, N, iter_max))
  }
  mu_mat <- matrix(NA, nrow=length(obs), ncol=N)
  sigma_array <- array(data = NA, dim = c(length(obs), length(obs), N))
  weight <- rep(-log(N), N)
  incremental_weight <- rep(log(1), N)
  ess_flat <- ESS_weight2(weight, incremental_weight)
  cess_flat <- CESS_weight(weight, incremental_weight)
  gamma_old <- 0
  iter <- 1
  if (gamma_history) {
    gamma_vec <- c(gamma_old)
    ess_vec <- c(ess_flat)
    cess_vec <- c(cess_flat)
  }

  # Initialization
  for (n in 1:N) {
    theta_n <- prior_sampler()
    sample_sta <- sample_func(theta_n, M)

    theta_mat[, n] <- theta_n
    mu_mat[, n] <- sample_sta$mean
    sigma_array[, , n] <- sample_sta$sigma
  }
  if (theta_history) {
    theta_history_array[, , iter] <- theta_mat
  }
  iter <- iter + 1


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
    search_u <- 1
    search_l <- gamma_old
    gamma_new <- (search_u + search_l) / 2
    # Reweight
    incremental_weight <- (gamma_new - gamma_old) * log_likelihood
    # incremental_weight <- incremental_weight - logSumExp(incremental_weight)
    cess_new <- CESS_weight(weight, incremental_weight)
    for (i in 1:100) {
      if (abs(cess_new - (log(alpha)+cess_flat)) < 0.001) {
        break
      }

      if (cess_new < (log(alpha)+cess_flat)) {
        # ESS too small
        search_u <- gamma_new
        gamma_new <- (search_u + search_l) / 2
      } else {
        # ESS too large
        search_l <- gamma_new
        gamma_new <- (search_u + search_l) / 2
        # Try gamma_new = 1 once
        if (i == 1) {
          incremental_weight <- (gamma_new - gamma_old) * log_likelihood
          # incremental_weight <- incremental_weight - logSumExp(incremental_weight)
          cess_new <- CESS_weight(weight, incremental_weight)
          if (cess_new >= (log(alpha)+cess_flat)) {
            gamma_new <- 1
          }
        }
      }

      # Reweight
      incremental_weight <- (gamma_new - gamma_old) * log_likelihood
      # incremental_weight <- incremental_weight - logSumExp(incremental_weight)
      cess_new <- CESS_weight(weight, incremental_weight)

      if (gamma_new == 1) {break}
    }
    ess_new <- ESS_weight2(weight, incremental_weight)
    weight <- weight + incremental_weight
    weight <- weight - logSumExp(weight)

    # Resample
    if (ess_new < (log(0.5)+ess_flat)) {
      prob_vec <- exp(weight)
      resample_index <- sample(1:N, size=N, replace=TRUE, prob=prob_vec)
      theta_mat <- theta_mat[, resample_index, drop=FALSE]
      mu_mat <- mu_mat[, resample_index, drop=FALSE]
      log_likelihood <- log_likelihood[resample_index]
      sigma_array <- sigma_array[, , resample_index, drop=FALSE]
      weight <- rep(-log(N), N)
    }

    # Move
    for (n in 1:N) {
      theta_new <- theta_mat[, n] + as.vector(rmvnorm(n=1, sigma=q_sigma))
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
    cess_old <- cess_new
    if (gamma_history) {
      gamma_vec <- c(gamma_vec, gamma_old)
      ess_vec <- c(ess_vec, ess_new)
      cess_vec <- c(cess_vec, cess_new)
    }
    if (theta_history) {
      theta_history_array[, , iter] <- theta_mat
    }
    iter <- iter + 1
    if (iter >= iter_max) {
      break
    }
  }

  if (theta_history) {
    result_list <- list(theta=theta_history_array)
  } else {
    result_list <- list(theta=theta_mat)
  }

  if (gamma_history) {
    result_list$gamma <- gamma_vec
    result_list$ess <- ess_vec
    result_list$cess <- cess_vec
  }

  return(result_list)
}
