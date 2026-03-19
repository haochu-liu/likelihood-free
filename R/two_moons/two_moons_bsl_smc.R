library(mvtnorm)
library(coda)
library(matrixStats)
library(rmatio)
source("R/BSL/SL_SMC.R")
source("R/BSL/ESS_weight2.R")
source("R/BSL/CESS_weight.R")


set.seed(100)
x_obs <- c(0, 0)

# BSL setup
prior_sampler <- function() {
  return(runif(2, min=-1, max=1))
}

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

M <- 7
alpha <- 0.9
N <- 1000

for (i in 1:10) {
  set.seed(i)
  result <- SL_SMC(M, alpha, N, 2, x_obs, prior_sampler, prior_func,
                   sample_func, cov_proposal)
  path_str <- paste0("output/two_moons/bsl_smc_post_sims_seed", i, ".csv")
  post_mat <- t(result$theta)
  write.table(post_mat, file=path_str, sep=",",
              row.names=FALSE, col.names=FALSE)
  print(paste0("Finish i = ", i))
}
