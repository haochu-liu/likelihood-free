library(mvtnorm)
library(coda)
library(matrixStats)
library(rmatio)
source("R/BSL/SL_MCMC2.R")


lambda <- 30
alpha <- 0.001
beta <- 0.001
N <- 100

# Load obs data
github_url <- "https://raw.githubusercontent.com/cdrovandi/Bayesian-Synthetic-Likelihood/master/Simple/data_poisson.mat"
data_list <- read.mat(github_url)

# Sampling and densities
y_obs <- data_list$y
s_obs <- mean(y_obs)
lambda_vec <- seq(0, 50, length.out = 50001)
prior_density <- dgamma(lambda_vec, shape=alpha, rate=beta)
posterior_density <- dgamma(lambda_vec, shape=alpha+sum(y_obs), rate=beta+N)
posterior_mean <- (alpha + sum(y_obs)) / (beta + N)

plot(lambda_vec, prior_density,
     ylim = c(0, max(posterior_density)),
     xlab = expression(lambda),
     ylab = "Density",
     type = "l",
     main="Prior and posterior densities",
     lwd = 2,
     col = "green")
lines(lambda_vec, posterior_density,
      col = "darkblue", lwd = 2)
abline(v = posterior_mean, col = "red", lwd = 2)
legend("topright",
       c("prior", "posterior", "posterior mean"),
       col = c("green", "darkblue", "red"),
       lwd = 2)

# BSL setup
init_theta <- lambda

prior_func <- function(theta){
  (alpha - 1) * log(theta) - beta * theta
}

sample_func <- function(theta, n) {
  repeat {
    x_mat <- matrix(rpois(n*N, theta), nrow=N, ncol=n)
    s <- colMeans(x_mat)
    if (var(s) > 0) {
      return(list(mean=c(mean(s)),
                  sigma=matrix(var(s), ncol=1, nrow=1)))
    }
  }
}
sample_func_fix_sigma <- function(theta, n) {
  x_mat <- matrix(rpois(n*N, theta), nrow=N, ncol=n)
  s <- colMeans(x_mat)
  return(list(mean=c(mean(s)),
              sigma=matrix(theta/N, ncol=1, nrow=1)))
}

sd_proposal <- sqrt(alpha + sum(y_obs)) / (beta + N)
proposal <- function(theta_old){
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

save(mcmc_quality1, file="mcmc_quality1.RData")

# Fix n*T = 100000, change n and T
set.seed(100)
mcmc_quality2 <- data.frame(n=c(1, 2, 5, 6, 7, 10, 15, 20),
                            acc_rate.fix_var=NA,
                            ess.fix_var=NA,
                            var_log.fix_var=NA,
                            acc_rate.est_var=NA,
                            ess.est_var=NA,
                            var_log.est_var=NA)
for (i in 1:nrow(mcmc_quality2)) {
  n <- mcmc_quality2$n[i]
  T_iter <- as.integer(100000 / n)

  acc_rate_vec <- rep(NA, 10)
  ess_vec <- rep(NA, 10)
  for (j in 1:10) {
    bsl_out <- SL_MCMC2(n, T_iter, s_obs, init_theta, prior_func, sample_func_fix_sigma,
                        proposal, acc_rate=TRUE)
    acc_rate_vec[j] <- bsl_out$acc_rate
    ess_vec[j] <- as.numeric(effectiveSize(as.mcmc(bsl_out$theta[1, ])))
  }

  mcmc_quality2$acc_rate.fix_var[i] <- mean(acc_rate_vec)
  mcmc_quality2$ess.fix_var[i] <- mean(ess_vec)

  log_like_vec <- rep(NA, 100)
  for (k in 1:100) {
    stats_n <- sample_func_fix_sigma(lambda, n)
    log_like_vec[k] <- dmvnorm(x=s_obs,
                               mean=stats_n$mean,
                               sigma=stats_n$sigma,
                               log=TRUE)
  }
  mcmc_quality2$var_log.fix_var[i] <- var(log_like_vec)

  if (n > 1) {
    acc_rate_vec <- rep(NA, 10)
    ess_vec <- rep(NA, 10)
    for (j in 1:10) {
      bsl_out <- SL_MCMC2(n, T_iter, s_obs, init_theta, prior_func, sample_func,
                          proposal, acc_rate=TRUE)
      acc_rate_vec[j] <- bsl_out$acc_rate
      ess_vec[j] <- as.numeric(effectiveSize(as.mcmc(bsl_out$theta[1, ])))
    }

    mcmc_quality2$acc_rate.est_var[i] <- mean(acc_rate_vec)
    mcmc_quality2$ess.est_var[i] <- mean(ess_vec)

    log_like_vec <- rep(NA, 100)
    for (k in 1:100) {
      stats_n <- sample_func(lambda, n)
      log_like_vec[k] <- dmvnorm(x=s_obs,
                                 mean=stats_n$mean,
                                 sigma=stats_n$sigma,
                                 log=TRUE)
    }
    mcmc_quality2$var_log.est_var[i] <- var(log_like_vec)
  }

  print(paste0("n = ", n, " finish."))
}

save(mcmc_quality2, file="mcmc_quality2.RData")


