library(mvtnorm)
library(coda)
library(matrixStats)
library(rmatio)
source("R/BSL/SL_MCMC2.R")


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
  repeat {
    x_mat <- matrix(rpois(n*N_val, theta), nrow=N_val, ncol=n)
    return(list(mean=rowMeans(x_mat),
                sigma=var(t(x_mat))))
  }
}
sample_func_fix_sigma <- function(theta, n) {
  x_mat <- matrix(rpois(n*N_val, theta), nrow=N_val, ncol=n)
  return(list(mean=rowMeans(x_mat),
              sigma=diag(30, nrow=N_val)))
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


# Fix T = 10000, change n
set.seed(100)
T_iter <- 10000
burn_in <- as.integer(T_iter/2)
n <- c(1, 2, 5, 6, 7, 10, 15, 20)
N <- c(1, 2, 5, 6, 7, 10, 15, 20)

acc_rate.fix_var <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(acc_rate.fix_var) <- N
rownames(acc_rate.fix_var) <- n

# acc_rate.est_var <- matrix(NA, nrow=length(n), ncol=length(N))
# colnames(acc_rate.est_var) <- N
# rownames(acc_rate.est_var) <- n

ess.fix_var <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(ess.fix_var) <- N
rownames(ess.fix_var) <- n

# ess.est_var <- matrix(NA, nrow=length(n), ncol=length(N))
# colnames(ess.est_var) <- N
# rownames(ess.est_var) <- n

var_log.fix_var <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(var_log.fix_var) <- N
rownames(var_log.fix_var) <- n

# var_log.est_var <- matrix(NA, nrow=length(n), ncol=length(N))
# colnames(var_log.est_var) <- N
# rownames(var_log.est_var) <- n

err_mean.fix_var <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(err_mean.fix_var) <- N
rownames(err_mean.fix_var) <- n


for (i in 1:length(N)) {
  N_val <- N[i]
  for (j in 1:length(n)) {
    n_val <- n[j]

    acc_rate_vec <- rep(NA, 10)
    ess_vec <- rep(NA, 10)
    err_mean_vec <- rep(NA, 10)
    for (k in 1:10) {
      bsl_out <- SL_MCMC2(n_val, T_iter, y_obs[1:N_val], init_theta, prior_func, sample_func_fix_sigma,
                          proposal, acc_rate=TRUE)
      acc_rate_vec[k] <- bsl_out$acc_rate
      ess_vec[k] <- as.numeric(effectiveSize(as.mcmc(bsl_out$theta[1, ])))
      err_mean_vec[k] = abs(posterior_mean(y_obs[1:N_val], N_val) - mean(bsl_out$theta[1, burn_in:T_iter]))
    }

    acc_rate.fix_var[j, i] <- mean(acc_rate_vec)
    ess.fix_var[j, i] <- mean(ess_vec)
    err_mean.fix_var[j, i] <- mean(err_mean_vec)

    log_like_vec <- rep(NA, 100)
    for (k in 1:100) {
      stats_n <- sample_func_fix_sigma(lambda, n_val)
      log_like_vec[k] <- dmvnorm(x=y_obs[1:N_val],
                                 mean=stats_n$mean,
                                 sigma=stats_n$sigma,
                                 log=TRUE)
    }
    var_log.fix_var[j, i] <- var(log_like_vec)

    # if (n_val > 1) {
    #   acc_rate_vec <- rep(NA, 10)
    #   ess_vec <- rep(NA, 10)
    #   for (k in 1:10) {
    #     bsl_out <- SL_MCMC2(n_val, T_iter, y_obs[1:N_val], init_theta, prior_func, sample_func,
    #                         proposal, acc_rate=TRUE)
    #     acc_rate_vec[k] <- bsl_out$acc_rate
    #     ess_vec[k] <- as.numeric(effectiveSize(as.mcmc(bsl_out$theta[1, ])))
    #   }
    #
    #   acc_rate.est_var[j, i] <- mean(acc_rate_vec)
    #   ess.est_var[j, i] <- mean(ess_vec)
    #
    #   log_like_vec <- rep(NA, 100)
    #   for (k in 1:100) {
    #     stats_n <- sample_func(lambda, n_val)
    #     log_like_vec[k] <- dmvnorm(x=y_obs[1:N_val],
    #                                mean=stats_n$mean,
    #                                sigma=stats_n$sigma,
    #                                log=TRUE)
    #   }
    #  var_log.est_var[j, i] <- var(log_like_vec)
    # }

    print(paste0("N = ", N_val, ", n = ", n_val, " finish."))
  }
}


norm_ess.fix_var <- ess.fix_var
for (i in 1:length(N)) {
  N_val <- N[i]
  for (j in 1:length(n)) {
    n_val <- n[j]
    norm_ess.fix_var[j, i] <- ess.fix_var[j, i] / (N_val * n_val)
  }
}

toy_fix_var <- list(acc_rate=acc_rate.fix_var,
                    ess=ess.fix_var,
                    norm_ess=norm_ess.fix_var,
                    var_log=var_log.fix_var,
                    err_mean=err_mean.fix_var)

save(toy_fix_var, file="data/toy_fix_var.RData")
