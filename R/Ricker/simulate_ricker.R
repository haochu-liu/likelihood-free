#' Simulate Ricker's model
#'
#' simulate_ricker simulates one data set from the model.
#'
#' @param theta A three dimension vector for the parameters.
#' @param N_t The starting population (equal to 1 in our application)
#' @param T_iter The length of the data set
#' @return The simulated data set of length T_iter.
simulate_ricker <- function(theta, N_t, T_iter) {
  y <- rep(NA, T_iter)
  r <- exp(theta[1])
  phi <- theta[2]
  sigmae <- theta[3]

  for (i in 1:T_iter) {
    e_t <- rnorm(1, mean=0, sd=sigmae)
    y[i] <- rpois(1, phi*r*N_t*exp(-N_t+e_t))
  }

  return(y)
}
