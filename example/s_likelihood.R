s_likelihood <- rep(NA, 481)
theta_vec <- seq(-12, 12, length.out = 481)
for (j in 1:481) {
  theta <- theta_vec[j]
  sl_vec <- rep(NA, 100)
  for (k in 1:100) {
    stats_M <- sample_func(theta, 20)
    sl_vec[k] <- dmvnorm(x=s_obs,
                         mean=stats_M$mean,
                         sigma=stats_M$sigma)
  }
  s_likelihood[j] <- mean(sl_vec)
}

prior <- dcauchy(theta_vec, location = 0, scale = 1)
prior <- prior / sum(prior) /
  (theta_vec[2] - theta_vec[1])
s_likelihood <- s_likelihood / sum(s_likelihood) /
  (theta_vec[2] - theta_vec[1])
posterior_unnorm <- prior * s_likelihood
posterior <- posterior_unnorm / sum(posterior_unnorm) /
  (theta_vec[2] - theta_vec[1])

plot(theta_vec, prior, type = "l", lwd = 2, col = "blue",
     ylim = c(0, 1.5),
     ylab = "Density",
     xlab = expression(theta),
     main = "Density functions")
lines(theta_vec, s_likelihood, col = "red", lwd = 2, lty = 2)
lines(theta_vec, posterior, col = "darkgreen", lwd = 3, lty = 1)
legend("topright",
       legend = c("Prior",
                  "N(s_obs; mu(theta), Sigma(theta))",
                  "pi_BSL(theta|s_obs)"),
       col = c("blue", "red", "darkgreen"),
       lwd = c(2, 2, 3),
       lty = c(1, 2, 1),
       bg = "white")

