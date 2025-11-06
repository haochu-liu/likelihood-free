library(moments)
library(mvtnorm)
library(coda)
library(matrixStats)
setwd("C:/Users/u2008181/likelihood-free/example")
source("../BSL/SL_MCMC.R")
source("../BSL/SL_AM.R")
source("../BSL/SL_SMC.R")
source("../BSL/ESS_weight.R")
set.seed(100)


sample_mixture_cauchy <- function(theta = 5, gamma = 1, n = 20) {
  u <- runif(n)
  locs <- ifelse(u < 0.5, theta, -theta)
  rcauchy(n, location = locs, scale = gamma)
}

summary_stats <- function(y) {
  mean_y <- mean(y)
  # median_y <- median(y)
  sd_y <- sd(y)
  mad_y <- mad(y)
  skew_y <- moments::skewness(y)
  kurt_y <- moments::kurtosis(y) + 3
  bc <- (skew_y^2 + 1) / kurt_y

  stats_vec <- c(mean_y, sd_y, mad_y, skew_y, kurt_y, bc)
  return(stats_vec)
}

n <- 100
y_obs <- sample_mixture_cauchy(theta=5, gamma=1, n=n)
s_obs <- summary_stats(y_obs)

iter <- 50000
init_theta <- rcauchy(1, location=0, scale=1)
prior_func <- function(theta){
  dcauchy(theta, location=0, scale=1, log=TRUE)
}
sample_func <- function(theta, M) {
  s_mat <- matrix(NA, nrow=6, ncol=M)
  for (i in 1:M) {
    z <- sample_mixture_cauchy(theta=theta, gamma=1, n=n)
    s_mat[, i] <- summary_stats(z)
  }

  return(list(mean=rowMeans(s_mat),
              sigma=var(t(s_mat))))
}

M_seq <- seq(5, 100, by=5)
log_likelihood <- rep(NA, length(M_seq))
for (j in 1:length(M_seq)) {
  M <- M_seq[j]
  sl_vec <- rep(NA, 100)
  for (i in 1:100) {
    stats_M <- sample_func(5, M)
    sl_vec[i] <- dmvnorm(x=s_obs,
                         mean=stats_M$mean,
                         sigma=stats_M$sigma,
                         log=TRUE)
  }
  log_likelihood[j] <- var(sl_vec)
  print(paste0("Finish M = ", M))
}

print(log_likelihood)
M <- M_seq[9]

q_sigma <- matrix(1, nrow=1, ncol=1)
theta_seq <- SL_MCMC(M, iter, s_obs, init_theta, prior_func, sample_func, q_sigma,
                     acc_rate=TRUE)

plot(1:iter, theta_seq$theta[1, ], type = "l", xlab="Iterations", ylab="Theta")

theta_density <- density(theta_seq$theta[1, 10000:iter])

plot(theta_density,
     main = "Density functions",
     xlab = "Theta",
     ylab = "Density",
     col = "black",
     lwd = 2)


prior_sampler <- function(){
  return(rcauchy(1, location=0, scale=1))
}
prior_func <- function(theta){
  dcauchy(theta, location=0, scale=1, log=TRUE)
}
sample_func <- function(theta, M) {
  s_mat <- matrix(NA, nrow=6, ncol=M)
  for (i in 1:M) {
    z <- sample_mixture_cauchy(theta=theta, gamma=1, n=n)
    s_mat[, i] <- summary_stats(z)
  }

  return(list(mean=rowMeans(s_mat),
              sigma=var(t(s_mat))))
}

alpha <- 0.9
N <- 10000
theta_d <- 1
q_sigma <- matrix(1, nrow=1, ncol=1)

theta_smc <- SL_SMC(M, alpha, N, theta_d, s_obs, prior_sampler, prior_func,
                    sample_func, q_sigma, gamma_history=TRUE)

theta_density <- density(theta_smc$theta[1, ])

plot(theta_density,
     main = "Density functions",
     xlab = "Theta",
     ylab = "Density",
     col = "black",
     lwd = 2)
