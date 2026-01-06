#' Effective sample size
#'
#' Compute the effective sample size for log normalized weights.
#'
#' @param W A vector of log-weights.
#' @return log of effective sample size.
ESS_weight <- function(W) {
  W_squared <- 2 * W
  log_ess <- -logSumExp(W_squared)
  return(log_ess)
}
