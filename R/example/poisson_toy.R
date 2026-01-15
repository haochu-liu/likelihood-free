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
    x_mat <- matrix(rpois(n*N, theta), nrow=N, ncol=n)
    if (var(s) > 0) {
      return(list(mean=rowMeans(x_mat),
                  sigma=var(t(x_mat))))
    }
  }
}
sample_func_fix_sigma <- function(theta, n) {
  x_mat <- matrix(rpois(n*N, theta), nrow=N, ncol=n)
  return(list(mean=rowMeans(x_mat),
              sigma=diag(30, nrow=N)))
}

proposal <- function(theta_old){
  sd_proposal <- sqrt(alpha + sum(y_obs)) / (beta + N)
  repeat {
    theta_new <- rnorm(1, mean=theta_old, sd=sd_proposal)
    if (theta_new > 0) {
      return(c(theta_new))
    }
  }
}

# Fix T = 100000, change n
set.seed(100)
T_iter <- 100000
mcmc_quality1 <- data.frame(n=c(1, 2, 5, 6, 7, 10, 15, 20),
                            acc_rate.fix_var=NA,
                            ess.fix_var=NA,
                            var_log.fix_var=NA,
                            acc_rate.est_var=NA,
                            ess.est_var=NA,
                            var_log.est_var=NA)
for (i in 1:nrow(mcmc_quality1)) {
  n <- mcmc_quality1$n[i]

  acc_rate_vec <- rep(NA, 10)
  ess_vec <- rep(NA, 10)
  for (j in 1:10) {
    bsl_out <- SL_MCMC2(n, T_iter, s_obs, init_theta, prior_func, sample_func_fix_sigma,
                        proposal, acc_rate=TRUE)
    acc_rate_vec[j] <- bsl_out$acc_rate
    ess_vec[j] <- as.numeric(effectiveSize(as.mcmc(bsl_out$theta[1, ])))
  }

  mcmc_quality1$acc_rate.fix_var[i] <- mean(acc_rate_vec)
  mcmc_quality1$ess.fix_var[i] <- mean(ess_vec)

  log_like_vec <- rep(NA, 100)
  for (k in 1:100) {
    stats_n <- sample_func_fix_sigma(lambda, n)
    log_like_vec[k] <- dmvnorm(x=s_obs,
                               mean=stats_n$mean,
                               sigma=stats_n$sigma,
                               log=TRUE)
  }
  mcmc_quality1$var_log.fix_var[i] <- var(log_like_vec)

  if (n > 1) {
    acc_rate_vec <- rep(NA, 10)
    ess_vec <- rep(NA, 10)
    for (j in 1:10) {
      bsl_out <- SL_MCMC2(n, T_iter, s_obs, init_theta, prior_func, sample_func,
                          proposal, acc_rate=TRUE)
      acc_rate_vec[j] <- bsl_out$acc_rate
      ess_vec[j] <- as.numeric(effectiveSize(as.mcmc(bsl_out$theta[1, ])))
    }

    mcmc_quality1$acc_rate.est_var[i] <- mean(acc_rate_vec)
    mcmc_quality1$ess.est_var[i] <- mean(ess_vec)

    log_like_vec <- rep(NA, 100)
    for (k in 1:100) {
      stats_n <- sample_func(lambda, n)
      log_like_vec[k] <- dmvnorm(x=s_obs,
                                 mean=stats_n$mean,
                                 sigma=stats_n$sigma,
                                 log=TRUE)
    }
    mcmc_quality1$var_log.est_var[i] <- var(log_like_vec)
  }

  print(paste0("n = ", n, " finish."))
}











