#' Synthetic Likelihood SMC (resampling every iteration)
#'
#' Apply sequential Monte Carlo (SMC) algorithm for BSL posterior target.
#' Resampling at every iteration.
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
SL_SMC_resample <- function(M, alpha, N, theta_d, obs, prior_sampler, prior_func,
                            sample_func, q_sigma, AM=TRUE,
                            theta_history=FALSE, gamma_history=FALSE,
                            acc_history=FALSE) {
  theta_mat <- matrix(NA, nrow=theta_d, ncol=N)
  iter_max <- 50
  if (theta_history) {
    theta_history_array <- array(data = NA, dim = c(theta_d, N, iter_max))
    weight_history_mat <- matrix(NA, nrow = iter_max, ncol = N)
  }
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
  if (acc_history) {
    acc_vec <- c()
    acc_count <- 0
  }

  # Initialization
  log_likelihood <- rep(NA, N)
  for (n in 1:N) {
    theta_n <- prior_sampler()
    sample_sta <- sample_func(theta_n, M)

    theta_mat[, n] <- theta_n
    log_likelihood[n] <- dmvnorm(x=obs,
                                 mean=sample_sta$mean,
                                 sigma=sample_sta$sigma,
                                 log=TRUE)
  }
  if (theta_history) {
    theta_history_array[, , iter] <- theta_mat
    weight_history_mat[iter, ] <- weight
  }
  iter <- iter + 1

  while (gamma_old < 1) {
    # Binary search 100 times
    search_u <- 1
    search_l <- gamma_old
    gamma_new <- (search_u + search_l) / 2
    # Reweight
    incremental_weight <- (gamma_new - gamma_old) * log_likelihood
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
          cess_new <- CESS_weight(weight, incremental_weight)
          if (cess_new >= (log(alpha)+cess_flat)) {
            gamma_new <- 1
          }
        }
      }

      # Reweight
      incremental_weight <- (gamma_new - gamma_old) * log_likelihood
      cess_new <- CESS_weight(weight, incremental_weight)

      if (gamma_new == 1) {break}
    }
    ess_new <- ESS_weight2(weight, incremental_weight)
    weight <- weight + incremental_weight
    weight <- weight - logSumExp(weight)
    if (AM) {
      s_d <- 2.38^2 / theta_d
      q_sigma <- s_d * cov.wt(t(theta_mat), wt=exp(weight))$cov
    }

    if ((gamma_new != 1) & iter < iter_max) {
      # No resample and move in the final iteration
      # Resample
      prob_vec <- exp(weight)
      resample_index <- sample(1:N, size=N, replace=TRUE, prob=prob_vec)
      theta_mat <- theta_mat[, resample_index, drop=FALSE]
      log_likelihood <- log_likelihood[resample_index]
      weight <- rep(-log(N), N)

      # Move
      for (n in 1:N) {
        theta_new <- theta_mat[, n] + as.vector(rmvnorm(n=1, sigma=q_sigma))
        stats_new <- sample_func(theta_new, M)
        sl_new <- dmvnorm(x=obs,
                          mean=stats_new$mean, sigma=stats_new$sigma, log=TRUE)

        log_alpha <- gamma_new * sl_new + prior_func(theta_new) -
          gamma_new * log_likelihood[n] - prior_func(theta_mat[, n])
        log_alpha <- min(0, log_alpha)
        log_u <- log(runif(1))

        if (log_u < log_alpha) {
          theta_mat[, n] <- theta_new
          log_likelihood[n] <- sl_new
          if (acc_history) {acc_count <- acc_count + 1}
        }
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
      weight_history_mat[iter, ] <- weight
    }
    if (acc_history) {
      acc_vec <- c(acc_vec, acc_count/N)
      acc_count <- 0
    }
    if (iter >= iter_max) {
      break
    }
    iter <- iter + 1
  }

  if (theta_history) {
    result_list <- list(theta=theta_history_array)
    result_list$weight <- weight_history_mat
  } else {
    result_list <- list(theta=theta_mat)
    result_list$weight <- weight
  }

  if (gamma_history) {
    result_list$gamma <- gamma_vec
    result_list$ess <- ess_vec
    result_list$cess <- cess_vec
  }

  if (acc_history) {
    result_list$acc <- acc_vec
  }

  return(result_list)
}
