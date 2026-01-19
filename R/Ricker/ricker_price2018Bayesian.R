library(mvtnorm)
library(coda)
library(matrixStats)
library(rmatio)
library(future.apply)


theta <- c(3.8, 10, 0.3)
N_0 <- 1

# Load obs data
github_url <- "https://raw.githubusercontent.com/cdrovandi/Bayesian-Synthetic-Likelihood/master/Ricker/data_ricker.mat"
data_list <- read.mat(github_url)

y_obs <- data_list$y
len_y <- length(y_obs)

# BSL setup
init_theta <- theta

prior_func <- function(theta){
  # independent, uniform, and improper prior
  return(0)
}

sample_func <- function(theta, n) {
  s_mat <- matrix(NA, nrow=13, ncol=n)
  for (i in 1:n) {
    x = simulate_ricker(theta, 1, 50)
    s_mat[, i] <- ricker_summstats(x, y_obs)
  }

  return(list(mean=rowMeans(s_mat),
              sigma=var(t(s_mat))))
}

std_proposal = diag(c(0.15, 0.5, 0.15))
corr_proposal = matrix(c(1, -0.7, -0.6,
                         -0.7, 1, 0.4,
                         -0.6, 0.4, 1), ncol=3, nrow=3)
cov_proposal <- std_proposal %*% corr_proposal %*% std_proposal
proposal <- function(theta_old){
  repeat {
    theta_new <- rmvnorm(1, mean=theta_old, sigma=cov_proposal)
    if (theta_new[3] >= 0) {
      return(theta_new)
    }
  }
}

# BSL-MCMC
SL_MCMC2 <- function(M, iter, obs, init_theta, prior_func, sample_func, proposal,
                     acc_rate=FALSE) {
  # Initial setup
  n_theta <- length(init_theta)
  n_obs <- length(obs)
  theta_matrix <- matrix(NA, nrow=n_theta, ncol=iter)
  if (acc_rate) {accept_num <- 0}
  i <- 1

  # Sample and likelihood at i = 1
  theta_old <- init_theta
  stats_old <- sample_func(theta_old, M)
  sl_old <- dmvnorm(x=obs, mean=stats_old$mean, sigma=stats_old$sigma, log=TRUE)
  theta_matrix[, i] <- init_theta

  # M-H MCMC
  for (i in 2:iter) {
    theta_new <- proposal(theta_old)
    stats_new <- sample_func(theta_new, M)
    sl_new <- dmvnorm(x=obs, mean=stats_new$mean, sigma=stats_new$sigma, log=TRUE)

    log_alpha <- sl_new + prior_func(theta_new) - sl_old - prior_func(theta_old)
    log_alpha <- min(0, log_alpha)
    log_u <- log(runif(1))

    if (log_u < log_alpha & !is.na(log_u < log_alpha)) {
      theta_matrix[, i] <- theta_new
      theta_old <- theta_new
      stats_old <- stats_new
      sl_old <- sl_new
      if (acc_rate) {accept_num <- accept_num + 1}
    } else {
      theta_matrix[, i] <- theta_old
    }
  }

  result_list <- list(theta=theta_matrix)
  if (acc_rate) {
    # print(paste0("Acceptance rate: ", accept_num/iter))
    result_list$acc_rate = accept_num/iter
  }
  return(result_list)
}

# Fix T = 500000, change n between 2 and 120
set.seed(100)
T_iter <- 500000
burn_in <- as.integer(T_iter/2)
n <- c(20, 30, 40, 50, 80, 100, 250)
s_obs <- ricker_summstats(y_obs, y_obs)
ricker_50 <- list(n=n,
                  acc_rate=rep(NA, length(n)),
                  ess=matrix(NA, nrow=3, ncol=length(n)),
                  norm_ess=matrix(NA, nrow=3, ncol=length(n)),
                  var_log_like=rep(NA, length(n)),
                  var_mean_square=rep(NA, length(n)),
                  post_mean=matrix(NA, nrow=3, ncol=length(n)))

plan(multisession, workers = 10)

for (i in 1:length(n)) {
  n_val <- n[i]

  # --- PARALLELIZE FIRST K LOOP ---
  results_k1 <- future_lapply(1:20, function(k) {
    bsl_out <- SL_MCMC2(n_val, T_iter, s_obs, init_theta,
                        prior_func, sample_func,
                        proposal, acc_rate=TRUE)
    ess_vec <- c(effectiveSize(as.mcmc(bsl_out$theta[1, ])),
                 effectiveSize(as.mcmc(bsl_out$theta[2, ])),
                 effectiveSize(as.mcmc(bsl_out$theta[3, ])))

    list(
      acc = bsl_out$acc_rate,
      ess = as.numeric(ess_vec),
      post_mean = rowMeans(bsl_out$theta[, burn_in:T_iter])
    )
  }, future.seed = TRUE)

  ricker_50$acc_rate[i]    <- mean(sapply(results_k1, `[[`, "acc"))
  ricker_50$ess[, i]       <- rowMeans(sapply(results_k1, `[[`, "ess"))
  ricker_50$norm_ess[, i]  <- ricker_50$ess[, i] / n_val
  ricker_50$post_mean[, i] <- rowMeans(sapply(results_k1, `[[`, "post_mean"))

  # --- PARALLELIZE SECOND K LOOP ---
  log_like_vec <- future_sapply(1:100, function(k) {
    stats_n <- sample_func(theta, n_val)
    dmvnorm(x = s_obs,
            mean = stats_n$mean,
            sigma = stats_n$sigma,
            log = TRUE)
  }, future.seed = TRUE)

  ricker_50$var_log_like[i] <- var(log_like_vec)
  like_vec <- exp(log_like_vec)
  ricker_50$var_mean_square[i] <- var(like_vec) / (mean(like_vec)^2)

  print(paste0("n = ", n_val, " finish."))
}

plan(sequential)

save(ricker_50, file="data/ricker_50.RData")


