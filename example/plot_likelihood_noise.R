x_value <- seq(-5, 5, length.out = 1001)
noise_vec <- rnorm(1001, mean=0, sd=0.01)
y_value <- dnorm(x_value, mean=0, sd=1) + noise_vec

plot(x_value, y_value,
     type = "l",
     col = "black",
     xlab = "theta",
     ylab = "f(x|theta)",
     lwd = 2)

plot(x_value, log(y_value),
     type = "l",
     col = "black",
     xlab = "theta",
     ylab = "log(f(x|theta))",
     lwd = 2)
