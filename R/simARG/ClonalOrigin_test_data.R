# devtools::install_github("haochu-liu/simARG")
library(simARG)


set.seed(100)
length_vec <- rep(NA, 10000)
height_vec <- rep(NA, 10000)
tree <- clonal_genealogy(10L)
rho_site <- 0.5
L <- 100L
delta <- 5
k <- 20L
for (i in 1:10000){
  ARG <- ClonalOrigin_pair_seq(tree, rho_site, L, delta, k)
  length_vec[i] <- sum(ARG$edge[, 3])
  height_vec[i] <- max(ARG$node_height)

  if (i %% 100 == 0) {
    print(paste("Complete", i, "iteration."))
  }
}

ARG_length_height <- data.frame(
  length = length_vec,
  height = height_vec
)

write.csv(ARG_length_height, file="pysimARG/test_data/ARG_length_height.csv",
          row.names=FALSE, col.names=FALSE)
