library(ggplot2)
library(GGally)


theta <- c(1.0, 0.01, 0.6)
state <- c(100, 50)
time_points <- c(0, seq(2, 64, by=2))
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

stats_mat <- matrix(NA, nrow=9, ncol=100)
for (i in 1:100) {
  set.seed(i)
  results <- Lotka_Volterra(theta, state, time_points, S)
  stats_mat[, i] <- L_V_summary_stats(results)
  print(paste0("Finish iteration ", i))
}

ggpairs(as.data.frame(stats_mat[, 1:25]))
