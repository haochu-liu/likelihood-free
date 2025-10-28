library(ggplot2)


theta <- c(1.0, 0.005, 0.6)
state <- c(100, 50)
time_points <- c(0, Measurement_times)
S <- matrix(c(
  1, 0,
  -1, 1,
  0, -1
), nrow = 3, byrow = TRUE)

set.seed(100)
results <- Lotka_Volterra(theta, state, time_points, S)

ggplot(results, aes(x = Time)) +
  geom_line(aes(y = Prey, color = "Prey")) +
  geom_line(aes(y = Predator, color = "Predator")) +
  geom_point(aes(y = Prey, color = "Prey")) +
  geom_point(aes(y = Predator, color = "Predator")) +
  labs(title = "Lotka-Volterra Simulation",
       y = "Population Size", x = "Time") +
  scale_color_manual(values = c("Prey" = "darkred", "Predator" = "darkblue"), name = "Species")

stats <- L_V_summary_stats(results)
print(stats)


