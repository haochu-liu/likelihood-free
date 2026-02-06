# devtools::install_github("haochu-liu/simARG")
library(simARG)


set.seed(100)
length_vec <- rep(NA, 10000)
height_vec <- rep(NA, 10000)
for (i in 1:10000){
  tree <- clonal_genealogy(10L)
  length_vec[i] <- sum(tree$edge[, 3])
  height_vec[i] <- max(tree$node)

  if (i %% 100 == 0) {
    print(paste("Complete", i, "iteration."))
  }
}

clonal_length_height <- data.frame(
  length = length_vec,
  height = height_vec
)

write.csv(clonal_length_height, file="pysimARG/test_data/clonal_length_height.csv",
          row.names=FALSE, col.names=FALSE)
