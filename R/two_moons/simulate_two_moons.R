simulate_two_moons <- function(theta) {
  mean_radius <- 0.1
  sd_radius <- 0.01
  baseoffset <- 0.25

  a <- pi * (runif(1) - 0.5)
  r <- mean_radius + rnorm(1, mean = 0, sd = 1) * sd_radius
  p <- c(r * cos(a) + baseoffset, r * sin(a))

  z0 <- -abs(theta[1] + theta[2]) / sqrt(2.0)
  z1 <- (-theta[1] + theta[2]) / sqrt(2.0)

  return(p + c(z0, z1))
}
