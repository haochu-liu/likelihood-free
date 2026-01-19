#' Summary statistics for Ricker model
#'
#' ricker_summstats computes the summary statistics used in the Ricker model (Wood, 2010).
#'
#' @param x A vector of simulated data.
#' @param y A vector of observed data.
#' @return Summary statistics for the simulated data x, a vector of dim 13.
ricker_summstats <- function(x, y) {
  len_x <- length(x)
  ss1x <- acf(x, lag.max = 5, type = "covariance", plot = FALSE)$acf
  lags <- 0:5
  ss1x <- ss1x * (len_x / (len_x - lags))

  order_diff <- 1
  x_diff <- diff(x, differences = order_diff)
  y_diff <- diff(y, differences = order_diff)
  x_diff <- sort(x_diff)
  y_diff <- sort(y_diff)
  x_diff2 <- x_diff - mean(x_diff)
  y_diff2 <- y_diff - mean(y_diff)

  fit_cubic <- lm(x_diff2 ~ 0 + y_diff2 + I(y_diff2^2) + I(y_diff2^3))
  ss2x <- as.vector(coef(fit_cubic))

  x_mod <- x^0.3
  x_mod <- x_mod - mean(x_mod)
  x_mod2 <- x_mod[2:length(x_mod)]
  x_pred2 <- x_mod[1:(length(x_mod) - 1)]

  fit_ar <- lm(x_mod2 ~ 0 + x_pred2 + I(x_pred2^2))
  ss3x <- as.vector(coef(fit_ar))

  ss4x <- mean(x)

  ss5x <- sum(x == 0)

  return(c(ss1x, ss2x, ss3x, ss4x, ss5x))
}
