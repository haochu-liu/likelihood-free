library(mvtnorm)
library(matrixStats)
set.seed(100)


obs <- rnbinom(1, size=5*n, prob=0.5)

prior_func <- function(){
  return(rgamma(1, shape=2, rate=0.5))
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
theta_history=FALSE

