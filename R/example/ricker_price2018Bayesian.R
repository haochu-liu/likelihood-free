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


