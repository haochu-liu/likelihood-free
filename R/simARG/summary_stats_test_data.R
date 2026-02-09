# devtools::install_github("haochu-liu/simARG")
library(simARG)
library(jsonlite)


set.seed(100)
r_vec <- rep(NA, 10)
g3_vec <- rep(NA, 10)
leaf_mat <- matrix(NA, nrow=10, ncol=20)
tree <- clonal_genealogy(10L)
rho_site <- 0.5
L <- 100L
delta <- 5
k <- 20L
tree_width <- tree$n
ARG <- ClonalOrigin_pair_seq(tree, rho_site, L, delta, k)
theta_site <- 0.5

for (i in 1:10) {
  ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
  mat <- ARG_mutated$node_site[1:tree_width, ]

  r_vec[i] <- LD_r(mat)
  g3_vec[i] <- G3_test(mat)
  leaf_mat[, (i*2-1):(i*2)] <- mat
}

summary_stats <- list(
  r = r_vec,
  g3 = g3_vec,
  leaf = leaf_mat
)

write_json(summary_stats, "pysimARG/test_data/summary_stats.json",
           pretty = TRUE, auto_unbox = TRUE)
