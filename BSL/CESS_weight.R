#' Conditional Effective sample size
#'
#' Compute the conditional effective sample size for log normalized weights using a more complex formula.
#'
#' @param W A vector of normalized log-weights of last iteration.
#' @param w A vector of un-normalized log-weights of current iteration.
#' @return log of effective sample size.
CESS_weight <- function(W, w) {
  N <- length(W)
  log_sess <- log(N) + 2 * logSumExp(W + w) - logSumExp(W + 2 * w)
  return(log_sess)
}
