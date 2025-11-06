library(moments)
library(mvtnorm)
library(coda)
library(matrixStats)
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
  median_y <- median(y)
  sd_y <- sd(y)
  mad_y <- mad(y)
  skew_y <- moments::skewness(y)
  kurt_y <- moments::kurtosis(y) + 3  # convert to Pearson kurtosis
  bc <- (skew_y^2 + 1) / kurt_y       # bimodality coefficient

  data.frame(
    mean = mean_y, median = median_y, sd = sd_y, mad = mad_y,
    skew = skew_y, kurtosis = kurt_y, bimodality_coeff = bc
  )
}

n <- 20
obs <- rnbinom(1, size=5*n, prob=0.5)

iter <- 50000
init_theta <- c(rgamma(1, shape=2, rate=0.5))
prior_func <- function(theta){
  dgamma(theta, shape=2, rate=0.5, log=TRUE)
}
sample_func <- function(theta, M) {
  s <- rpois(M, theta*n)
  return(list(mean=c(mean(s)),
              sigma=matrix(var(s), ncol=1, nrow=1)))
}



