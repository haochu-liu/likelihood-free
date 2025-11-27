library(moments)
library(mvtnorm)
library(coda)
library(matrixStats)
library(future.apply)
source("../BSL/SL_SMC.R")
source("../BSL/SL_SMC_resample.R")
source("../BSL/SL_IBIS.R")
source("../BSL/SL_WF_SMC.R")
source("../BSL/SL_WF_IBIS.R")
source("../BSL/SL_WF_SMC_par.R")
source("../BSL/SL_WF_IBIS_par.R")
source("../BSL/ESS_weight2.R")
source("../BSL/CESS_weight.R")
set.seed(100)


sample_mix_cauchy <- function(theta = 5, n = 20) {
  u <- runif(n)
  locs <- ifelse(u < 0.5, theta, -theta)
  rcauchy(n, location = locs, scale = 1)
}

summary_stats <- function(y) {
  median_y <- median(y)
  mad_y <- mad(y) # median absolute deviation
  quantile_y <- quantile(y, probs = c(0.25, 0.75))
  median_abs_y <- median(abs(y))
  stats_vec <- c(median_y, mad_y, quantile_y[2] - quantile_y[1], median_abs_y)
  return(stats_vec)
}

y_obs <- sample_mix_cauchy()
s_obs <- summary_stats(y_obs)

prior_sampler <- function(){
  rcauchy(1, location=0, scale=1)
}
prior_func <- function(theta){
  dcauchy(theta, location=0, scale=1, log=TRUE)
}
sample_func <- function(theta, M) {
  s_mat <- matrix(NA, nrow=4, ncol=M)
  for (i in 1:M) {
    z <- sample_mix_cauchy(theta=theta)
    s_mat[, i] <- summary_stats(z)
  }

  return(list(mean=rowMeans(s_mat),
              sigma=var(t(s_mat))))
}

M <- 20
alpha <- 0.9
N <- 1000
N_sample <- 10
theta_d <- 1
q_sigma <- matrix(1, nrow=1, ncol=1)
plan(multisession, workers = N_sample)

df1 <- data.frame(type=rep(c("SMC", "SMC(R)", "IBIS"), each=50),
                  move=rep("std", 150),
                  run=rep("sequential", 150),
                  time=rep(NA, 150),
                  QU=rep(NA, 150),
                  QL=rep(NA, 150),
                  E=rep(NA, 150))

df2 <- data.frame(type=rep(c("SMC(R)", "IBIS"), each=50),
                  move=rep("waste-free", 100),
                  run=rep("sequential", 100),
                  time=rep(NA, 100),
                  QU=rep(NA, 100),
                  QL=rep(NA, 100),
                  E=rep(NA, 100))

df3 <- data.frame(type=rep(c("SMC(R)", "IBIS"), each=50),
                  move=rep("waste-free", 100),
                  run=rep("parallel", 100),
                  time=rep(NA, 100),
                  QU=rep(NA, 100),
                  QL=rep(NA, 100),
                  E=rep(NA, 100))

df_eff <- rbind(df1, df2)
df_eff <- rbind(df_eff, df3)

for (i in 1:350) {
  func_type <- df_eff[i, 1:3]
  if (func_type$type == "SMC") {
    time_result <- system.time(
      theta_smc <- SL_SMC(M, alpha, N, theta_d, s_obs, prior_sampler,
                          prior_func, sample_func, q_sigma, AM=TRUE,
                          theta_history=FALSE, gamma_history=FALSE)
    )
  } else if (func_type$type == "SMC(R)" & func_type$move == "std") {
    time_result <- system.time(
      theta_smc <- SL_SMC_resample(M, alpha, N, theta_d, s_obs, prior_sampler,
                                   prior_func, sample_func, q_sigma, AM=TRUE,
                                   theta_history=FALSE, gamma_history=FALSE)
    )
  } else if (func_type$type == "IBIS" & func_type$move == "std") {
    time_result <- system.time(
      theta_smc <- SL_IBIS(M, alpha, N, theta_d, s_obs, prior_sampler,
                           prior_func, sample_func, q_sigma, AM=TRUE,
                           theta_history=FALSE, gamma_history=FALSE)
    )
  } else if (func_type$type == "SMC(R)" &
             func_type$move == "waste-free" &
             func_type$run == "sequential") {
    time_result <- system.time(
      theta_smc <- SL_WF_SMC(M, alpha, N, N_sample, theta_d, s_obs, prior_sampler,
                             prior_func, sample_func, q_sigma, AM=TRUE,
                             theta_history=FALSE, gamma_history=FALSE)
    )
  } else if (func_type$type == "IBIS" &
             func_type$move == "waste-free" &
             func_type$run == "sequential") {
    time_result <- system.time(
      theta_smc <- SL_WF_IBIS(M, alpha, N, N_sample, theta_d, s_obs, prior_sampler,
                              prior_func, sample_func, q_sigma, AM=TRUE,
                              theta_history=FALSE, gamma_history=FALSE)
    )
  } else if (func_type$type == "SMC(R)" &
             func_type$move == "waste-free" &
             func_type$run == "parallel") {
    time_result <- system.time(
      theta_smc <- SL_WF_SMC.par(M, alpha, N, N_sample, theta_d, s_obs,
                                 prior_sampler, prior_func, sample_func,
                                 q_sigma, AM=TRUE, theta_history=FALSE,
                                 gamma_history=FALSE, future_seed=NULL)
    )
  } else if (func_type$type == "IBIS" &
             func_type$move == "waste-free" &
             func_type$run == "parallel") {
    time_result <- system.time(
      theta_smc <- SL_WF_IBIS.par(M, alpha, N, N_sample, theta_d, s_obs,
                                  prior_sampler, prior_func, sample_func,
                                  q_sigma, AM=TRUE, theta_history=FALSE,
                                  gamma_history=FALSE, future_seed=NULL)
    )
  }

  df_eff$time[i] <- as.numeric(time_result["elapsed"])
  df_eff$QU[i] <- quantile(theta_smc$theta[1, ], probs=0.95)
  df_eff$QL[i] <- quantile(theta_smc$theta[1, ], probs=0.05)
  df_eff$E[i] <- mean(theta_smc$theta[1, ])

  if (i %% 5 == 0) {print(paste0("Finish iteration ", i))}
}

saveRDS(df_eff, file="df_eff.rds")
