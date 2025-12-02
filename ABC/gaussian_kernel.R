#' Gaussian kernel function
#'
#' Give the function value for a Gaussian kernel.
#'
#' @param y One vector.
#' @param z One vector.
#' @param tol The tolerance epsilon.
#' @param sigma The covariance matrix.
#' @param log.kernel If TRUE, return value in log scale.
#' @return Function value.
#' @export
gaussian_kernel <- function(y, z, tol=1, sigma, log.kernel=TRUE) {
  if (length(y)!=length(z)) {
    stop("Two vectors should have same length.")
  }

  d <- as.matrix(y - z)
  u <- sqrt(t(d) %*% solve(sigma) %*% d) / tol
  if (log.kernel) {
    return(dnorm(u, log=TRUE))
  } else {
    return(dnorm(u, log=FALSE))
  }
}
