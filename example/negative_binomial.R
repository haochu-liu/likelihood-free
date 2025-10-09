library(mvtnorm)
source("BSL/SL_MCMC.R")


set.seed(100)
n <- 20
y <- rnbinom(n, size=5, prob=0.5)
obs <- c(mean(y))

M <- 20
iter <- 10000
init_theta <- c(rgamma(1, shape=2, rate=0.5))
prior_func <- function(theta){
  dgamma(theta, shape=2, rate=0.5, log=TRUE)
  }
sample_func <- function(theta, M) {
  s_z <- rpois(M, theta)
  return(list(mean=c(mean(s_z)),
              sigma=matrix(var(s_z), ncol=1, nrow=1)))
}

theta_seq <- SL_MCMC(M, iter, obs, init_theta, prior_func, sample_func, 1)


# Trace plot
plot(1:iter, theta_seq[1, ], type = "l", xlab="Iterations", ylab="Theta")

# Density plot




