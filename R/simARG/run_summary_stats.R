library(simARG)
library(future.apply)
library(progressr)


set.seed(100)
tree <- clonal_genealogy(15L)

x_mat <- matrix(NA, nrow=200, ncol=3)
x_mat[, 1] <- runif(200, min=0, max=0.2)   # rho_site
x_mat[, 2] <- runif(200, min=1, max=500)   # delta
x_mat[, 3] <- runif(200, min=0, max=0.2)   # theta_site

plan(multisession, workers = 10)

handlers(handler_txtprogressbar)
handlers(global = TRUE)

print("Starting simulation...")

summary_stats_list <- with_progress({
  p <- progressor(along = 1:nrow(x_mat))

  result <- future_lapply(1:nrow(x_mat), function(i) {

    rho_site   <- x_mat[i, 1]
    delta      <- x_mat[i, 2]
    theta_site <- x_mat[i, 3]

    s_vec <- ClonalOrigin_pair_seq.simulator_fast(
      tree,
      rho_site,
      theta_site,
      1e6L,
      delta,
      2000
    )

    p()

    return(s_vec)

  }, future.seed = TRUE, future.packages = "simARG")
})

summary_stats_mat <- do.call(rbind, summary_stats_list)

plan(sequential)

write.csv(x_mat, file = "x_mat.csv", row.names = FALSE, col.names = FALSE)
write.csv(summary_stats_mat, file = "summary_stats_mat.csv", row.names = FALSE, col.names = FALSE)

print("Finish!")
