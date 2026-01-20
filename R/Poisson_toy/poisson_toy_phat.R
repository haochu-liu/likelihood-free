library(mvtnorm)
library(matrixStats)
library(rmatio)
library(future.apply)


lambda <- 30
alpha <- 0.001
beta <- 0.001

# Load obs data
github_url <- "https://raw.githubusercontent.com/cdrovandi/Bayesian-Synthetic-Likelihood/master/Simple/data_poisson.mat"
data_list <- read.mat(github_url)

# Sampling and densities
y_obs <- data_list$y

# Setup
sample_func_fix_sigma <- function(theta, n) {
  x_mat <- matrix(rpois(n*N_val, theta), nrow=N_val, ncol=n)
  return(list(mean=rowMeans(x_mat),
              sigma=diag(30, nrow=N_val)))
}


# Change n and N
set.seed(100)
n <- c(1, 2, 5, 6, 7, 10, 15, 20)
N <- c(1, 2, 5, 6, 7, 10, 15, 20)

var_log.fix_var <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(var_log.fix_var) <- N
rownames(var_log.fix_var) <- n

est.fix_var <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(est.fix_var) <- N
rownames(est.fix_var) <- n

log_est.fix_var <- matrix(NA, nrow=length(n), ncol=length(N))
colnames(log_est.fix_var) <- N
rownames(log_est.fix_var) <- n

plan(multisession, workers = 10)

for (i in 1:length(N)) {
  N_val <- N[i]
  for (j in 1:length(n)) {
    n_val <- n[j]

    # --- PARALLELIZE SECOND K LOOP ---
    log_like_vec <- future_sapply(1:100, function(k) {
      stats_n <- sample_func_fix_sigma(lambda, n_val)
      dmvnorm(x = y_obs[1:N_val],
              mean = stats_n$mean,
              sigma = stats_n$sigma,
              log = TRUE)
    }, future.seed = TRUE)

    like_vec <- exp(log_like_vec)
    var_log.fix_var[j, i] <- var(log_like_vec)
    est.fix_var[j, i] <- var(like_vec) / (mean(like_vec)^2)
    log_est.fix_var[j, i] <- log(var(like_vec)) - 2*log(mean(like_vec))

    print(paste0("N = ", N_val, ", n = ", n_val, " finish."))
  }
}

plan(sequential)


toy_like_fix_var <- list(var_log=var_log.fix_var,
                         est=est.fix_var,
                         log_est=log_est.fix_var)

save(toy_like_fix_var, file="data/toy_like_fix_var.RData")
