#' Lotka-Volterra Model simulation
#'
#' Simulate Lotka-Volterra model with Gillespie's method.
#'
#' @param theta A three dimensional vector for the parameters.
#' @param state A two dimensional integer vector for the initial state of prey and predator.
#' @param time_points A vector of time points (in ascending order) for measurement.
#' @param S A transition matrix for changing states.
#' @return A dataframe which records the numbers of predator and prey at every time point.
Lotka_Volterra <- function(theta, state, time_points, S) {
  time <- 0
  num_time <- length(time_points)
  time_final <- time_points[num_time]
  results <- data.frame(Time = time_points,
                        Prey = c(state[1], rep(0, num_time-1)),
                        Predator = c(state[2], rep(0, num_time-1)))
  measurement_index <- 2
  i <- 1

  while (time < time_final) {
    h <- c(theta[1]*state[1],
           theta[2]*state[1]*state[2],
           theta[3]*state[2])
    h0 <- sum(h)

    if (h0 == 0) {
      time <- time_final
      break
    }

    tau <- rexp(1, rate = h0) # Time step tau ~ Exp(h0)
    j <- sample.int(nrow(S), size = 1, prob = h/h0) # Select situation j
    time_new <- time + tau

    # Record measurements within the time step
    while ((measurement_index <= num_time) &
           (time_points[measurement_index] <= time_new)) {
      results[measurement_index, 2:3] <- state
      measurement_index <- measurement_index + 1
    }

    # Update state
    state <- state + S[j, ]
    time <- time_new
    i <- i + 1
  }
  return(results)
}
