#' Synthetic Likelihood SMC
#'
#' Apply sequential Monte Carlo (SMC) algorithm for BSL posterior as the target distribution.
#'
#' @param M Number of new data points drawn in each iteration.
#' @param alpha Number to control the effective sample size in reweight.
#' @param N Number of particles for SMC.
#' @param obs A vector of the observed statistics.
#' @param prior_func A function to sample parameters from prior.
#' @param sample_func A function which takes theta and M and return sample mean and variance.
#' @param sigma A scale parameter for the Gaussian proposal in move step.
#' @param history Default history = FALSE, if TRUE, return all particles in history.
#' @return A vector of parameters from the BSL posterior.
SL_MCMC <- function() {

}
