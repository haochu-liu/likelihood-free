library(mvtnorm)
library(coda)
library(matrixStats)
library(rmatio)
source("R/BSL/SL_MCMC2.R")


theta <- c(3.8, 10, 0.3)
N_0 <- 1

# Load obs data
github_url <- "https://raw.githubusercontent.com/cdrovandi/Bayesian-Synthetic-Likelihood/master/Ricker/data_ricker.mat"
data_list <- read.mat(github_url)

y_obs <- data_list$y
T_iter <- length(y_obs)

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



