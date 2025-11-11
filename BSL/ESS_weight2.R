#' Effective sample size
#'
#' Compute the effective sample size for log normalized weights using a more complex formula.
#'
#' @param W A vector of normalized log-weights of last iteration.
#' @param w A vector of un-normalized log-weights of current iteration.
#' @return log of effective sample size.
ESS_weight2 <- function(W, w) {
  log_add <- W + w
  log_ess <- 2 * logSumExp(log_add) - logSumExp(2 * log_add)
  return(log_ess)
}
