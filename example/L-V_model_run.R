theta <- c(1.0, 0.01, 0.6)
state <- c(100, 50)
time_points <- c(0, Measurement_times)
S <- matrix(c(
  1, 0,
  -1, 1,
  0, -1
), nrow = 3, byrow = TRUE)

set.seed(100)
results <- Lotka_Volterra(theta, state, time_points, S)
