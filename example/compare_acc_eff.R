library(moments)
library(mvtnorm)
library(coda)
library(matrixStats)
library(future.apply)
library(ggplot2)
source("BSL/SL_MCMC.R")
source("BSL/SL_SMC.R")
source("BSL/SL_SMC_resample.R")
source("BSL/SL_WF_SMC.R")
source("BSL/SL_WF_SMC_par.R")
source("ABC/ABC_MCMC.R")
source("ABC/ABC_SMC.R")
source("ABC/gaussian_kernel.R")
source("ABC/gaussian_kernel2.R")
source("BSL/ESS_weight2.R")
source("BSL/CESS_weight.R")
source("wasserstein1.R")
source("MMD.R")
set.seed(100)


# Get obs data
sample_normal <- function(theta = 5, n = 20) {
  rnorm(n, mean = theta, sd = 1)
}

summary_stats <- function(y) {
  mean_y <- mean(y)
  stats_vec <- c(mean_y)
  return(stats_vec)
}

y_obs <- sample_normal()
s_obs <- summary_stats(y_obs)

# Theoretical posterior
# posterior: N(9*s_obs/(0.05+9), 1/(1/0.05+1/9))
post_mean <- 9*s_obs/(0.05+9)
post_sd <- sqrt(1/(1/0.05+1/9))

# Setup
M <- 10
init_theta <- rnorm(1, mean=0, sd=3)
prior_sampler <- function(){
  rnorm(1, mean=0, sd=3)
}
prior_func <- function(theta){
  dnorm(theta, mean=0, sd=3, log=TRUE)
}
sample_bsl <- function(theta, M) {
  s_mat <- matrix(NA, nrow=1, ncol=M)
  for (i in 1:M) {
    z <- sample_normal(theta=theta)
    s_mat[, i] <- summary_stats(z)
  }

  return(list(mean=rowMeans(s_mat),
              sigma=var(t(s_mat))))
}
sample_abc_mcmc <- function(theta) {
  z <- sample_normal(theta=theta)
  return(summary_stats(z))
}
sample_abc_smc <- function(theta, M) {
  s_mat <- matrix(NA, nrow=1, ncol=M)
  for (i in 1:M) {
    z <- sample_normal(theta=theta)
    s_mat[, i] <- summary_stats(z)
  }

  return(s_mat)
}

q_wd <- qnorm((1:1000)/1001, mean = post_mean, sd = post_sd)

y_mmd <- matrix(rnorm(1000, mean = post_mean, sd = post_sd), ncol=1000)

N_sample = 10
plan(multisession, workers = N_sample)

# Construct df
df <- data.frame(method=c("ABC-MCMC", "ABC-SMC", "SL-MCMC", "SL-SMC",
                          "SL-SMC(R)", "SL-WF-SMC", "SL-WF-SMC.par"),
                 Wasserstein_dist = NA,
                 MMD = NA,
                 bias.post_mean = NA,
                 mse.post_mean = NA,
                 bias.post_sd = NA,
                 mse.post_sd = NA,
                 time = NA,
                 VU.time = NA,
                 VL.time = NA,
                 VE.time = NA)

# Run
num_runs <- 20
measure_mat <- matrix(NA, nrow=num_runs, ncol=7)
colnames(measure_mat) <- c("Wasserstein_dist", "MMD", "est_mean", "est_sd",
                           "time", "QU", "QL")

# ABC-MCMC
set.seed(100)
iter <- 1
theta_sigma <- matrix(0.1, nrow=1, ncol=1)
q_sigma <- matrix(1, nrow=1, ncol=1)
for (i in 1:num_runs) {
  time_result <- system.time(
    theta <- ABC_MCMC(0.1, 50000, s_obs, gaussian_kernel, theta_sigma,
                      init_theta, prior_func, sample_abc_mcmc, q_sigma,
                      acc_rate=FALSE)
  )

  X_vec <- theta$theta[1, 49001:50000]
  measure_mat[i, 1] <- wasserstein1(X_vec, q_wd)
  measure_mat[i, 2] <- MMD(matrix(X_vec, ncol=1000), y_mmd)
  measure_mat[i, 3] <- mean(X_vec)
  measure_mat[i, 4] <- sd(X_vec)
  measure_mat[i, 5] <- as.numeric(time_result["elapsed"])
  measure_mat[i, 6] <- quantile(X_vec, probs=0.95)
  measure_mat[i, 7] <- quantile(X_vec, probs=0.05)

  if (i %% 5 == 0) {print(paste0("Finish iteration ", i))}
}
df$Wasserstein_dist[iter] <- mean(measure_mat[, 1])
df$MMD[iter] <- mean(measure_mat[, 2])
df$bias.post_mean[iter] <- abs(mean(measure_mat[, 3]) - post_mean)
df$mse.post_mean[iter] <- mean((measure_mat[, 3] - post_mean)^2)
df$bias.post_sd[iter] <- abs(mean(measure_mat[, 4]) - post_sd)
df$mse.post_sd[iter] <- mean((measure_mat[, 4] - post_sd)^2)
df$time[iter] <- mean(measure_mat[, 5])
df$VU.time[iter] <- var(measure_mat[, 6]) * mean(measure_mat[, 5])
df$VL.time[iter] <- var(measure_mat[, 7]) * mean(measure_mat[, 5])
df$VE.time[iter] <- var(measure_mat[, 3]) * mean(measure_mat[, 5])

# ABC-SMC
set.seed(100)
iter <- 2
theta_sigma <- matrix(0.1, nrow=1, ncol=1)
q_sigma <- matrix(1, nrow=1, ncol=1)
for (i in 1:num_runs) {
  time_result <- system.time(
    theta <- ABC_SMC(1, 0.01, 0.9, M, 1000, 1, s_obs,
                     gaussian_kernel2, theta_sigma, prior_sampler, prior_func,
                     sample_abc_smc, q_sigma, AM=TRUE)
  )

  X_vec <- theta$theta[1, ]
  measure_mat[i, 1] <- wasserstein1(X_vec, q_wd)
  measure_mat[i, 2] <- MMD(matrix(X_vec, ncol=1000), y_mmd)
  measure_mat[i, 3] <- mean(X_vec)
  measure_mat[i, 4] <- sd(X_vec)
  measure_mat[i, 5] <- as.numeric(time_result["elapsed"])
  measure_mat[i, 6] <- quantile(X_vec, probs=0.95)
  measure_mat[i, 7] <- quantile(X_vec, probs=0.05)

  if (i %% 5 == 0) {print(paste0("Finish iteration ", i))}
}
df$Wasserstein_dist[iter] <- mean(measure_mat[, 1])
df$MMD[iter] <- mean(measure_mat[, 2])
df$bias.post_mean[iter] <- abs(mean(measure_mat[, 3]) - post_mean)
df$mse.post_mean[iter] <- mean((measure_mat[, 3] - post_mean)^2)
df$bias.post_sd[iter] <- abs(mean(measure_mat[, 4]) - post_sd)
df$mse.post_sd[iter] <- mean((measure_mat[, 4] - post_sd)^2)
df$time[iter] <- mean(measure_mat[, 5])
df$VU.time[iter] <- var(measure_mat[, 6]) * mean(measure_mat[, 5])
df$VL.time[iter] <- var(measure_mat[, 7]) * mean(measure_mat[, 5])
df$VE.time[iter] <- var(measure_mat[, 3]) * mean(measure_mat[, 5])

# BSL-MCMC
set.seed(100)
iter <- 3
q_sigma <- matrix(1, nrow=1, ncol=1)
for (i in 1:num_runs) {
  time_result <- system.time(
    theta <- SL_MCMC(M, 50000, s_obs, init_theta, prior_func, sample_bsl,
                     q_sigma, acc_rate=FALSE)
  )

  X_vec <- theta$theta[1, 49001:50000]
  measure_mat[i, 1] <- wasserstein1(X_vec, q_wd)
  measure_mat[i, 2] <- MMD(matrix(X_vec, ncol=1000), y_mmd)
  measure_mat[i, 3] <- mean(X_vec)
  measure_mat[i, 4] <- sd(X_vec)
  measure_mat[i, 5] <- as.numeric(time_result["elapsed"])
  measure_mat[i, 6] <- quantile(X_vec, probs=0.95)
  measure_mat[i, 7] <- quantile(X_vec, probs=0.05)

  if (i %% 5 == 0) {print(paste0("Finish iteration ", i))}
}
df$Wasserstein_dist[iter] <- mean(measure_mat[, 1])
df$MMD[iter] <- mean(measure_mat[, 2])
df$bias.post_mean[iter] <- abs(mean(measure_mat[, 3]) - post_mean)
df$mse.post_mean[iter] <- mean((measure_mat[, 3] - post_mean)^2)
df$bias.post_sd[iter] <- abs(mean(measure_mat[, 4]) - post_sd)
df$mse.post_sd[iter] <- mean((measure_mat[, 4] - post_sd)^2)
df$time[iter] <- mean(measure_mat[, 5])
df$VU.time[iter] <- var(measure_mat[, 6]) * mean(measure_mat[, 5])
df$VL.time[iter] <- var(measure_mat[, 7]) * mean(measure_mat[, 5])
df$VE.time[iter] <- var(measure_mat[, 3]) * mean(measure_mat[, 5])

# BSL-SMC
set.seed(100)
iter <- 4
q_sigma <- matrix(1, nrow=1, ncol=1)
for (i in 1:num_runs) {
  time_result <- system.time(
    theta <- SL_SMC(M, 0.9, 1000, 1, s_obs, prior_sampler,
                    prior_func, sample_bsl, q_sigma, AM=TRUE)
  )

  X_vec <- theta$theta[1, ]
  measure_mat[i, 1] <- wasserstein1(X_vec, q_wd)
  measure_mat[i, 2] <- MMD(matrix(X_vec, ncol=1000), y_mmd)
  measure_mat[i, 3] <- mean(X_vec)
  measure_mat[i, 4] <- sd(X_vec)
  measure_mat[i, 5] <- as.numeric(time_result["elapsed"])
  measure_mat[i, 6] <- quantile(X_vec, probs=0.95)
  measure_mat[i, 7] <- quantile(X_vec, probs=0.05)

  if (i %% 5 == 0) {print(paste0("Finish iteration ", i))}
}
df$Wasserstein_dist[iter] <- mean(measure_mat[, 1])
df$MMD[iter] <- mean(measure_mat[, 2])
df$bias.post_mean[iter] <- abs(mean(measure_mat[, 3]) - post_mean)
df$mse.post_mean[iter] <- mean((measure_mat[, 3] - post_mean)^2)
df$bias.post_sd[iter] <- abs(mean(measure_mat[, 4]) - post_sd)
df$mse.post_sd[iter] <- mean((measure_mat[, 4] - post_sd)^2)
df$time[iter] <- mean(measure_mat[, 5])
df$VU.time[iter] <- var(measure_mat[, 6]) * mean(measure_mat[, 5])
df$VL.time[iter] <- var(measure_mat[, 7]) * mean(measure_mat[, 5])
df$VE.time[iter] <- var(measure_mat[, 3]) * mean(measure_mat[, 5])

# BSL-SMC(R)
set.seed(100)
iter <- 5
q_sigma <- matrix(1, nrow=1, ncol=1)
for (i in 1:num_runs) {
  time_result <- system.time(
    theta <- SL_SMC_resample(M, 0.9, 1000, 1, s_obs, prior_sampler,
                             prior_func, sample_bsl, q_sigma, AM=TRUE)
  )

  X_vec <- theta$theta[1, ]
  measure_mat[i, 1] <- wasserstein1(X_vec, q_wd)
  measure_mat[i, 2] <- MMD(matrix(X_vec, ncol=1000), y_mmd)
  measure_mat[i, 3] <- mean(X_vec)
  measure_mat[i, 4] <- sd(X_vec)
  measure_mat[i, 5] <- as.numeric(time_result["elapsed"])
  measure_mat[i, 6] <- quantile(X_vec, probs=0.95)
  measure_mat[i, 7] <- quantile(X_vec, probs=0.05)

  if (i %% 5 == 0) {print(paste0("Finish iteration ", i))}
}
df$Wasserstein_dist[iter] <- mean(measure_mat[, 1])
df$MMD[iter] <- mean(measure_mat[, 2])
df$bias.post_mean[iter] <- abs(mean(measure_mat[, 3]) - post_mean)
df$mse.post_mean[iter] <- mean((measure_mat[, 3] - post_mean)^2)
df$bias.post_sd[iter] <- abs(mean(measure_mat[, 4]) - post_sd)
df$mse.post_sd[iter] <- mean((measure_mat[, 4] - post_sd)^2)
df$time[iter] <- mean(measure_mat[, 5])
df$VU.time[iter] <- var(measure_mat[, 6]) * mean(measure_mat[, 5])
df$VL.time[iter] <- var(measure_mat[, 7]) * mean(measure_mat[, 5])
df$VE.time[iter] <- var(measure_mat[, 3]) * mean(measure_mat[, 5])

# BSL-WF-SMC(R)
set.seed(100)
iter <- 6
q_sigma <- matrix(1, nrow=1, ncol=1)
for (i in 1:num_runs) {
  time_result <- system.time(
    theta <- SL_WF_SMC(M, 0.9, 1000, N_sample, 1, s_obs, prior_sampler,
                       prior_func, sample_bsl, q_sigma, AM=TRUE)
  )

  X_vec <- theta$theta[1, ]
  measure_mat[i, 1] <- wasserstein1(X_vec, q_wd)
  measure_mat[i, 2] <- MMD(matrix(X_vec, ncol=1000), y_mmd)
  measure_mat[i, 3] <- mean(X_vec)
  measure_mat[i, 4] <- sd(X_vec)
  measure_mat[i, 5] <- as.numeric(time_result["elapsed"])
  measure_mat[i, 6] <- quantile(X_vec, probs=0.95)
  measure_mat[i, 7] <- quantile(X_vec, probs=0.05)

  if (i %% 5 == 0) {print(paste0("Finish iteration ", i))}
}
df$Wasserstein_dist[iter] <- mean(measure_mat[, 1])
df$MMD[iter] <- mean(measure_mat[, 2])
df$bias.post_mean[iter] <- abs(mean(measure_mat[, 3]) - post_mean)
df$mse.post_mean[iter] <- mean((measure_mat[, 3] - post_mean)^2)
df$bias.post_sd[iter] <- abs(mean(measure_mat[, 4]) - post_sd)
df$mse.post_sd[iter] <- mean((measure_mat[, 4] - post_sd)^2)
df$time[iter] <- mean(measure_mat[, 5])
df$VU.time[iter] <- var(measure_mat[, 6]) * mean(measure_mat[, 5])
df$VL.time[iter] <- var(measure_mat[, 7]) * mean(measure_mat[, 5])
df$VE.time[iter] <- var(measure_mat[, 3]) * mean(measure_mat[, 5])

# BSL-WF-SMC(R)
set.seed(100)
iter <- 7
q_sigma <- matrix(1, nrow=1, ncol=1)
for (i in 1:num_runs) {
  time_result <- system.time(
    theta <- SL_WF_SMC.par(M, 0.9, 1000, N_sample, 1, s_obs,
                           prior_sampler, prior_func, sample_bsl,
                           q_sigma, AM=TRUE, future_seed=NULL)
  )

  X_vec <- theta$theta[1, ]
  measure_mat[i, 1] <- wasserstein1(X_vec, q_wd)
  measure_mat[i, 2] <- MMD(matrix(X_vec, ncol=1000), y_mmd)
  measure_mat[i, 3] <- mean(X_vec)
  measure_mat[i, 4] <- sd(X_vec)
  measure_mat[i, 5] <- as.numeric(time_result["elapsed"])
  measure_mat[i, 6] <- quantile(X_vec, probs=0.95)
  measure_mat[i, 7] <- quantile(X_vec, probs=0.05)

  if (i %% 5 == 0) {print(paste0("Finish iteration ", i))}
}
df$Wasserstein_dist[iter] <- mean(measure_mat[, 1])
df$MMD[iter] <- mean(measure_mat[, 2])
df$bias.post_mean[iter] <- abs(mean(measure_mat[, 3]) - post_mean)
df$mse.post_mean[iter] <- mean((measure_mat[, 3] - post_mean)^2)
df$bias.post_sd[iter] <- abs(mean(measure_mat[, 4]) - post_sd)
df$mse.post_sd[iter] <- mean((measure_mat[, 4] - post_sd)^2)
df$time[iter] <- mean(measure_mat[, 5])
df$VU.time[iter] <- var(measure_mat[, 6]) * mean(measure_mat[, 5])
df$VL.time[iter] <- var(measure_mat[, 7]) * mean(measure_mat[, 5])
df$VE.time[iter] <- var(measure_mat[, 3]) * mean(measure_mat[, 5])

# Save df
saveRDS(df, file="df_acc_eff.rds")



