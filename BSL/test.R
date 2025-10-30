library(mvtnorm)
library(matrixStats)
set.seed(100)


obs <- rnbinom(1, size=5*n, prob=0.5)

prior_sampler <- function(){
  return(rgamma(1, shape=2, rate=0.5))
}

prior_func <- function(theta){
  dgamma(theta, shape=2, rate=0.5, log=TRUE)
}

sample_func <- function(theta, M) {
  s <- rpois(M, theta*20)
  return(list(mean=c(mean(s)),
              sigma=matrix(var(s), ncol=1, nrow=1)))
}

M <- 25
alpha <- 0.9
N <- 100
theta_d <- 1
obs
sigma <- 0.1
theta_history=FALSE

results <- SL_SMC(M, alpha, N, theta_d, obs, prior_sampler, prior_func,
                  sample_func, sigma)
