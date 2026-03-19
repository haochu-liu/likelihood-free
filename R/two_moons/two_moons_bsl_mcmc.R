library(mvtnorm)
library(coda)
library(matrixStats)
library(rmatio)


set.seed(100)
x_obs <- c(0, 0)

# BSL setup
init_theta <- c(0, 0)

prior_func <- function(theta){
  return(0)
}

sample_func <- function(theta, n) {
  x_mat <- matrix(NA, nrow=2, ncol=n)
  for (i in 1:n) {
    x_mat[, i] <- simulate_two_moons(theta)
  }

  return(list(mean=rowMeans(x_mat),
              sigma=var(t(x_mat))))
}

cov_proposal <- matrix(c(0.05, 0, 0, 0.05), ncol=2, nrow=2)
proposal <- function(theta_old){
  theta_new <- rmvnorm(1, mean=theta_old, sigma=cov_proposal)
  return(theta_new)
}

# Find optimal n
M_seq <- seq(5, 20, by=1)
log_likelihood <- rep(NA, length(M_seq))
names(log_likelihood) <- M_seq
for (j in 1:length(M_seq)) {
  M <- M_seq[j]
  sl_vec <- rep(NA, 100)
  for (i in 1:100) {
    stats_M <- sample_func(c(0.22, 0.22), M)
    sl_vec[i] <- dmvnorm(x=x_obs,
                         mean=stats_M$mean,
                         sigma=stats_M$sigma,
                         log=TRUE)
  }
  log_likelihood[j] <- var(sl_vec)
  print(paste0("Finish M = ", M))
}

print(log_likelihood)
M <- 7

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


for (i in 1:10) {
  set.seed(i)
  result <- SL_MCMC2(M, 10000, x_obs, init_theta, prior_func, sample_func, proposal)
  path_str <- paste0("output/two_moons/bsl_mcmc_post_sims10000_seed", i, ".csv")
  post_mat <- unique(t(result$theta))
  write.table(post_mat[101:nrow(post_mat), ], file=path_str, sep=",",
              row.names=FALSE, col.names=FALSE)
  print(paste0("Finish i = ", i))
}
