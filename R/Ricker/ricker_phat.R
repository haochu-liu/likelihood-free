library(mvtnorm)
library(coda)
library(matrixStats)
library(rmatio)
library(future.apply)
source("R/Ricker/simulate_ricker.R")
source("R/Ricker/ricker_summstats.R")


theta <- c(3.8, 10, 0.3)
N_0 <- 1

# Load obs data
github_url <- "https://raw.githubusercontent.com/cdrovandi/Bayesian-Synthetic-Likelihood/master/Ricker/data_ricker.mat"
data_list <- read.mat(github_url)

y_obs <- data_list$y
len_y <- length(y_obs)

# Setup
sample_func <- function(theta, n) {
  s_mat <- matrix(NA, nrow=13, ncol=n)
  for (i in 1:n) {
    x = simulate_ricker(theta, 1, 50)
    s_mat[, i] <- ricker_summstats(x, y_obs)
  }

  return(list(mean=rowMeans(s_mat),
              sigma=var(t(s_mat))))
}

# Change n between 20 and 250
set.seed(100)
n <- c(20, 30, 40, 50, 80, 100, 250)
s_obs <- ricker_summstats(y_obs, y_obs)
var_log_like <- rep(NA, length(n))
names(var_log_like) <- n
var_mean_square <- rep(NA, length(n))
names(var_mean_square) <- n
log_var_mean_square <- rep(NA, length(n))
names(log_var_mean_square) <- n

plan(multisession, workers = 10)

for (i in 1:length(n)) {
  n_val <- n[i]

  # --- PARALLELIZE SECOND K LOOP ---
  log_like_vec <- future_sapply(1:10000, function(k) {
    stats_n <- sample_func(theta, n_val)
    dmvnorm(x = s_obs,
            mean = stats_n$mean,
            sigma = stats_n$sigma,
            log = TRUE)
  }, future.seed = TRUE)

  var_log_like[i] <- var(log_like_vec)
  like_vec <- exp(log_like_vec)
  var_mean_square[i] <- var(like_vec) / (mean(like_vec)^2)
  log_var_mean_square[i] <- log(var(like_vec)) - 2*log(mean(like_vec))

  print(paste0("n = ", n_val, " finish."))
}

plan(sequential)


ricker_like <- list(var_log=var_log_like,
                    est=var_mean_square,
                    log_est=log_var_mean_square)

save(ricker_like, file="data/ricker_like.RData")
