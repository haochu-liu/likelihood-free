library(ggplot2)


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



theta_candidate <- c(0, 1, 5, -5)
s_df <- data.frame(theta=rep(theta_candidate, each=1000),
                   s1=NA,
                   s2=NA,
                   s3=NA,
                   s4=NA)
s_df_index <- 1
for (i in 1:length(theta_candidate)) {
  theta <- theta_candidate[i]
  for (j in 1:1000) {
    stats_M <- sample_func(theta, 20)
    s_df[s_df_index, 2:5] <- stats_M$mean
    s_df_index <- s_df_index + 1
  }
}

s_df$theta <- as.character(s_df$theta)

ggplot(s_df, aes(x = s1, fill = theta)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = s_obs[1],
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(
    title = "mu(theta) for theta = -5, 0, 1, 5",
    x = "mu_1(theta)",
    y = "Density",
    fill = "theta"
  ) +
  theme_minimal()

ggplot(s_df, aes(x = s2, fill = theta)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = s_obs[2],
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(
    title = "mu(theta) for theta = -5, 0, 1, 5",
    x = "mu_2(theta)",
    y = "Density",
    fill = "theta"
  ) +
  theme_minimal()

ggplot(s_df, aes(x = s3, fill = theta)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = s_obs[3],
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  coord_cartesian(xlim = c(0, 12)) +
  labs(
    title = "mu(theta) for theta = -5, 0, 1, 5",
    x = "mu_3(theta)",
    y = "Density",
    fill = "theta"
  ) +
  theme_minimal()

ggplot(s_df, aes(x = s4, fill = theta)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = s_obs[4],
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  coord_cartesian(xlim = c(0, 6)) +
  labs(
    title = "mu(theta) for theta = -5, 0, 1, 5",
    x = "mu_4(theta)",
    y = "Density",
    fill = "theta"
  ) +
  theme_minimal()



theta_candidate <- c(0, 1, 5, -5)
s_df <- data.frame(theta=rep(theta_candidate, each=1000),
                   s1=NA,
                   s2=NA,
                   s3=NA,
                   s4=NA)
s_df_index <- 1
for (i in 1:length(theta_candidate)) {
  theta <- theta_candidate[i]
  for (j in 1:1000) {
    y_sample <- sample_mix_cauchy(theta=theta)
    s_sample <- summary_stats(y_sample)
    s_df[s_df_index, 2:5] <- s_sample
    s_df_index <- s_df_index + 1
  }
}

s_df$theta <- as.character(s_df$theta)

ggplot(s_df, aes(x = s1, fill = theta)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = s_obs[1],
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(
    title = "s1(theta) for theta = -5, 0, 1, 5",
    x = "s1(theta)",
    y = "Density",
    fill = "theta"
  ) +
  theme_minimal()

ggplot(s_df, aes(x = s2, fill = theta)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = s_obs[2],
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(
    title = "s2(theta) for theta = -5, 0, 1, 5",
    x = "s2(theta)",
    y = "Density",
    fill = "theta"
  ) +
  theme_minimal()

ggplot(s_df, aes(x = s3, fill = theta)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = s_obs[3],
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  coord_cartesian(xlim = c(0, 12)) +
  labs(
    title = "s3(theta) for theta = -5, 0, 1, 5",
    x = "s3(theta)",
    y = "Density",
    fill = "theta"
  ) +
  theme_minimal()

ggplot(s_df, aes(x = s4, fill = theta)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = s_obs[4],
             linetype = "dashed",
             color = "red",
             linewidth = 1) +
  coord_cartesian(xlim = c(0, 6)) +
  labs(
    title = "s4(theta) for theta = -5, 0, 1, 5",
    x = "s4(theta)",
    y = "Density",
    fill = "theta"
  ) +
  theme_minimal()

