#' Maximum mean discrepancy
#'
#' Compute maximum mean discrepancy (MMD) between two vector of samples with RBF kernel.
#'
#' @param X One matrix which each column is a sample.
#' @param Y One matrix which each column is a sample.
#' @param sigma One numerical value for the RBF kernel.
#' @return Estimated MMD
MMD <- function(X, Y, sigma=1) {
  n <- ncol(X)
  m <- ncol(Y)

  sum_kxx <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      sum_kxx <- sum_kxx + exp(-sum((X[, i] - X[, j])^2) / (2 * sigma^2))
    }
  }

  sum_kyy <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      sum_kyy <- sum_kyy + exp(-sum((Y[, i] - Y[, j])^2) / (2 * sigma^2))
    }
  }

  sum_kxy <- 0
  for (i in 1:n) {
    for (j in 1:m) {
      sum_kxy <- sum_kxy + exp(-sum((X[, i] - Y[, j])^2) / (2 * sigma^2))
    }
  }

  mmd2_value <- sum_kxx / n^2 + sum_kyy / m^2 - sum_kxy *2 / (m * n)
  return(sqrt(mmd2_value))
}
