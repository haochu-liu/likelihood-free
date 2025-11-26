theta_seq <- theta_wf_ibis
par(mfrow = c(5, 1), mar = c(2, 4, 2, 1))
for (i in 1:5) {
  theta_density <- density(theta_seq$theta[1, ,i], bw = 0.5,
                           weights = exp(theta_seq$weight[i, ]))
  plot(theta_density,
       main = paste0("Iteration ", i),
       xlab = "Theta",
       ylab = "Density",
       xlim = c(-20, 20),
       col = "black",
       lwd = 2)
}
par(mfrow = c(1, 1))

par(mfrow = c(gamma_n-5, 1), mar = c(2, 4, 2, 1))
for (i in 6:gamma_n) {
  theta_density <- density(theta_seq$theta[1, ,i], bw = 0.5,
                           weights = exp(theta_seq$weight[i, ]))
  plot(theta_density,
       main = paste0("Iteration ", i),
       xlab = "Theta",
       ylab = "Density",
       xlim = c(-20, 20),
       col = "black",
       lwd = 2)
}
par(mfrow = c(1, 1))
