library(simARG)
library(future.apply)
library(progressr)


set.seed(100)
tree <- clonal_genealogy(10L)

param_mat <- matrix(NA, nrow=3, ncol=1000)
summary_stats_mat <- matrix(NA, nrow=7, ncol=1000)

param_mat[1, ] <- runif(1000, min=0, max=0.2)  # rho_site
param_mat[2, ] <- runif(1000, min=1, max=500)  # delta
param_mat[3, ] <- runif(1000, min=0, max=0.2)  # theta_site

handlers(global = TRUE)
handlers("progress")

plan(multisession, workers = 10)

with_progress({
  p <- progressor(along = 1:1000)
  results_list <- future_lapply(1:1000, function(i) {
    result <- ClonalOrigin_pair_seq.simulator(
      tree,
      param_mat[1, i],
      param_mat[3, i],
      100000L,
      param_mat[2, i],
      100,
      c(50L, 200L, 2000L)
    )
    p(message = sprintf("Finish %g iteration", i))

    return(result)
  }, future.seed = TRUE)
})

plan(sequential)

summary_stats_mat <- do.call(cbind, results_list)


plot(param_mat[1, ], summary_stats_mat[1, ], pch=4, col="darkred")
points(param_mat[1, ], summary_stats_mat[2, ], pch=4, col="darkblue")
points(param_mat[1, ], summary_stats_mat[3, ], pch=4, col="darkgreen")

plot(param_mat[2, ], summary_stats_mat[1, ], pch=4, col="darkred")
points(param_mat[2, ], summary_stats_mat[2, ], pch=4, col="darkblue")
points(param_mat[2, ], summary_stats_mat[3, ], pch=4, col="darkgreen")

plot(param_mat[3, ], summary_stats_mat[1, ], pch=4, col="darkred")
points(param_mat[3, ], summary_stats_mat[2, ], pch=4, col="darkblue")
points(param_mat[3, ], summary_stats_mat[3, ], pch=4, col="darkgreen")


plot(param_mat[1, ], summary_stats_mat[4, ], pch=4, col="darkred")
points(param_mat[1, ], summary_stats_mat[5, ], pch=4, col="darkblue")
points(param_mat[1, ], summary_stats_mat[6, ], pch=4, col="darkgreen")

plot(param_mat[2, ], summary_stats_mat[4, ], pch=4, col="darkred")
points(param_mat[2, ], summary_stats_mat[5, ], pch=4, col="darkblue")
points(param_mat[2, ], summary_stats_mat[6, ], pch=4, col="darkgreen")

plot(param_mat[3, ], summary_stats_mat[4, ], pch=4, col="darkred")
points(param_mat[3, ], summary_stats_mat[5, ], pch=4, col="darkblue")
points(param_mat[3, ], summary_stats_mat[6, ], pch=4, col="darkgreen")
