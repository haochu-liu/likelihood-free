library(mvtnorm)
library(coda)
library(matrixStats)
library(rmatio)
library(future.apply)


lambda <- 30
alpha <- 0.001
beta <- 0.001

# Load obs data
github_url <- "https://raw.githubusercontent.com/cdrovandi/Bayesian-Synthetic-Likelihood/master/Simple/data_poisson.mat"
data_list <- read.mat(github_url)

# Sampling and densities
y_obs <- data_list$y

# BSL setup
init_theta <- lambda

prior_func <- function(theta){
  (alpha - 1) * log(theta) - beta * theta
}

sample_func <- function(theta, n) {
  x_mat <- matrix(rpois(n*N_val, theta), nrow=N_val, ncol=n)
  return(list(mean=rowMeans(x_mat),
              sigma=var(t(x_mat))))
}

proposal <- function(theta_old){
  sd_proposal <- sqrt(alpha + sum(y_obs)) / (beta + N_val)
  repeat {
    theta_new <- rnorm(1, mean=theta_old, sd=sd_proposal)
    if (theta_new > 0) {
      return(c(theta_new))
    }
  }
}

posterior_mean <- function(y, N_val) {
  return((alpha + sum(y)) / (beta + N_val))
}

posterior_var <- function(y, N_val) {
  return((alpha + sum(y)) / (beta + N_val)^2)
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


# Fix T = 10000, change n
set.seed(100)
T_iter <- 10000
burn_in <- as.integer(T_iter/2)
n <- seq(from=5, to=20, by=5) # simulations
N <- seq(from=2, to=10, by=2)  # dimensions
# n <- seq(from=5, to=200, by=5) # simulations
# N <- seq(from=2, to=30, by=2)  # dimensions

acc_rate <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(acc_rate) <- N
rownames(acc_rate) <- n

ess <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(ess) <- N
rownames(ess) <- n

err_mean <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(err_mean) <- N
rownames(err_mean) <- n

err_var <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(err_var) <- N
rownames(err_var) <- n

var_log <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(var_log) <- N
rownames(var_log) <- n

relative_var <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(relative_var) <- N
rownames(relative_var) <- n


plan(multisession, workers = 10)

for (i in 1:length(N)) {
  N_val <- N[i]
  for (j in 1:length(n)) {
    n_val <- n[j]

    # BSL-MCMC
    results_k <- future_lapply(1:10, function(k) {
      bsl_out <- SL_MCMC2(n_val, T_iter, y_obs[1:N_val], init_theta,
                          prior_func, sample_func,
                          proposal, acc_rate=TRUE)
      kde_res <- density(bsl_out$theta[1, burn_in:T_iter])

      list(
        acc = bsl_out$acc_rate,
        ess = as.numeric(effectiveSize(as.mcmc(bsl_out$theta[1, ]))),
        err_mean = abs(posterior_mean(y_obs[1:N_val], N_val) - mean(bsl_out$theta[1, burn_in:T_iter])),
        err_var = abs(posterior_var(y_obs[1:N_val], N_val) - var(bsl_out$theta[1, burn_in:T_iter]))
      )
    }, future.seed = TRUE)

    acc_rate[j, i] <- mean(sapply(results_k, `[[`, "acc"))
    ess[j, i]      <- mean(sapply(results_k, `[[`, "ess"))
    err_mean[j, i] <- mean(sapply(results_k, `[[`, "err_mean"))
    err_var[j, i] <- mean(sapply(results_k, `[[`, "err_var"))

    # Compute var log-likelihood
    log_like_vec <- future_sapply(1:100, function(k) {
      stats_n <- sample_func(lambda, n_val)
      dmvnorm(x = y_obs[1:N_val],
              mean = stats_n$mean,
              sigma = stats_n$sigma,
              log = TRUE)
    }, future.seed = TRUE)

    like_vec <- exp(log_like_vec)
    var_log[j, i] <- var(log_like_vec)
    relative_var[j, i] <- var(like_vec) / (mean(like_vec)^2)

    print(paste0("N = ", N_val, ", n = ", n_val, " finish."))
  }
}

plan(sequential)


norm_ess <- ess
for (j in 1:length(n)) {
  norm_ess[j, ] <- ess[j, ] / n[j]
}

poisson_toy_table <- list(acc_rate=acc_rate,
                           ess=ess,
                           norm_ess=norm_ess,
                           err_mean=err_mean,
                           err_var=err_var,
                           var_log=var_log,
                           relative_var=relative_var)

save(poisson_toy_table, file="data/poisson_toy_table.RData")
