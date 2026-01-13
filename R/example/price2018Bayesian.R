library(mvtnorm)
library(rjags)
library(coda)
library(matrixStats)
source("R/BSL/SL_MCMC.R")
set.seed(100)


lambda <- 30
alpha <- 0.001
beta <- 0.001
N <- 100

# Sampling and densities
obs <- rpois(N, lambda)
lambda_vec <- seq(0, 50, length.out = 50001)
prior_density <- dgamma(lambda_vec, shape=alpha, rate=beta)
posterior_density <- dgamma(lambda_vec, shape=alpha+sum(obs), rate=beta+N)
posterior_mean <- (alpha + sum(obs)) / (beta + N)

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
init_theta <- c(rgamma(1, shape=alpha, rate=beta))
prior_func <- function(theta){
  dgamma(theta, shape=alpha, rate=beta, log=TRUE)
}
sample_func <- function(theta, n) {
  x <- rpois(n, theta)
  return(list(mean=c(mean(x)),
              sigma=matrix(var(s), ncol=1, nrow=1)))
}

















