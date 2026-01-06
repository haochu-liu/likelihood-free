#' Gaussian kernel function 2
#'
#' Give the function value for a Gaussian kernel.
#' Take normalised distance as input.
#'
#' @param u A numerical value.
#' @param tol The tolerance epsilon.
#' @param log.kernel If TRUE, return value in log scale.
#' @return Function value.
gaussian_kernel2 <- function(u, tol=1, log.kernel=TRUE) {
  v <- u / tol
  if (log.kernel) {
    return(dnorm(v, log=TRUE))
  } else {
    return(dnorm(v, log=FALSE))
  }
}
