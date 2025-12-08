#' Wasserstein distance (1D)
#'
#' Compute the Wasserstein distance between one dimensional samples and quantiles.
#'
#' @param X One vector, the samples from distribution P.
#' @param q One vector, the queantiles from distribution Q.
#' @return W_1(P, Q)
wasserstein1 <- function(X, q) {
  X_sorted <- sort(X)
  # Example: q <- qnorm((1:n)/(n+1), mean = 0, sd = 1).
  W1 <- mean(abs(X_sorted - q))
  return(W1)
}
