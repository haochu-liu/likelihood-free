library(mvtnorm)
source("BSL/SL_MCMC.R")


set.seed(100)
n <- 20
y <- rnbinom(n, size=5, prob=0.5)
obs <- c(mean(y))

M <- 20
iter <- 50000
init_theta <- c(rgamma(1, shape=2, rate=0.5))
prior_func <- function(theta){
  dgamma(theta, shape=2, rate=0.5, log=TRUE)
  }
sample_func <- function(theta, M) {
  s_z <- rpois(M, theta)
  return(list(mean=c(mean(s_z)),
              sigma=matrix(var(s_z), ncol=1, nrow=1)))
}

theta_seq <- SL_MCMC(M, iter, obs, init_theta, prior_func, sample_func, 0.1)


# Trace plot
plot(1:iter, theta_seq[1, ], type = "l", xlab="Iterations", ylab="Theta")

# Density plot
theta_density <- density(theta_seq[1, 10000:iter])
x_values <- seq(min(theta_seq[1, 10000:iter]),
                max(theta_seq[1, 10000:iter]),
                length.out = 1000)
y_values <- dgamma(x_values,
                   shape=(2 + obs),
                   rate=(0.5 + 1))
plot(theta_density,
     main = "BSL posterior",
     xlab = "Theta",
     ylab = "Density",
     col = "black",
     lwd = 2)
lines(x_values, y_values, col = "blue", lwd = 2, lty = 2)
legend("topright",
       c("Synthetic posterior", "Exact posterior"),
       col = c("black", "blue"),
       lty = c(1, 2),
       lwd = 2)


