library(simARG)


set.seed(100)
tree <- clonal_genealogy(10L)

param_mat <- matrix(NA, nrow=3, ncol=1000)
summary_stats_mat <- matrix(NA, nrow=7, ncol=1000)

param_mat[1, ] <- runif(1000, min=0, max=0.2)  # rho_site
param_mat[2, ] <- runif(1000, min=0, max=20)   # delta
param_mat[3, ] <- runif(1000, min=0, max=0.2)  # theta_site

for (i in 1:1000) {
  result <- ClonalOrigin_pair_seq.simulator(tree,
                                            param_mat[1, i],
                                            param_mat[3, i],
                                            100L,
                                            param_mat[2, i],
                                            100,
                                            c(2L, 10L, 75L))

  summary_stats_mat[, i] <- result
  print(paste("Finish", i, "iteration."))
}


