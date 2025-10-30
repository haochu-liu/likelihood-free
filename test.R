M <- 25
alpha <- 0.9
N <- 100
theta_d <- 1
obs

prior_func <- function(){
  return(rgamma(1, shape=2, rate=0.5))
}

sample_func <- function(theta, M) {
  s <- rpois(M, theta*n)
  return(list(mean=c(mean(s)),
              sigma=matrix(var(s), ncol=1, nrow=1)))
}

