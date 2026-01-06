#' Summary statistics for Lotka-Volterra Model results
#'
#' Provide a 9-dimensional vector composed of the mean, log variance,
#' first two auto-correlations of each time series,
#' and the cross-correlation between them.
#'
#' @param results A dataframe coming from `Lotka_Volterra()` function.
#' @return A 9-dimensional vector for the summary statistics.
L_V_summary_stats <- function(results) {
  summary_stats <- rep(NA, 9)

  # Summary statistics first 4 dimensions
  summary_stats[1] <- mean(results$Prey)
  summary_stats[2] <- mean(results$Predator)
  summary_stats[3] <- log(var(results$Prey))
  summary_stats[4] <- log(var(results$Predator))

  # Summary statistics 5-8 dimensions
  prey_acf <- acf(results$Prey, plot = FALSE, lag.max = 2)
  predator_acf <- acf(results$Predator, plot = FALSE, lag.max = 2)
  summary_stats[5] <- prey_acf$acf[2]
  summary_stats[6] <- prey_acf$acf[3]
  summary_stats[7] <- predator_acf$acf[2]
  summary_stats[8] <- predator_acf$acf[3]

  # Summary statistics last dimension
  series_ccf <- ccf(results$Prey, results$Predator, plot = FALSE, lag.max = 0)
  summary_stats[9] <- series_ccf$acf[1]

  return(summary_stats)
}
