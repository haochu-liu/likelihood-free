library(moments)
library(mvtnorm)
library(coda)
library(matrixStats)
library(future.apply)
library(ggplot2)
source("BSL/SL_MCMC.R")
source("BSL/SL_SMC.R")
source("BSL/SL_SMC_resample.R")
source("BSL/SL_WF_SMC.R")
source("BSL/SL_WF_SMC_par.R")
source("ABC/ABC_MCMC.R")
source("ABC/ABC_SMC.R")
source("ABC/gaussian_kernel.R")
source("ABC/gaussian_kernel2.R")
source("BSL/ESS_weight2.R")
source("BSL/CESS_weight.R")
source("wasserstein1.R")
source("MMD.R")
set.seed(100)


# Get obs data
sample_normal <- function(theta = 5, n = 20) {
  rnorm(n, mean = theta, sd = 1)
}

summary_stats <- function(y) {
  mean_y <- mean(y)
  stats_vec <- c(mean_y)
  return(stats_vec)
}

y_obs <- sample_normal()
s_obs <- summary_stats(y_obs)

# Theoretical posterior
# posterior: N(9*s_obs/(0.05+9), 1/(1/0.05+1/9))












