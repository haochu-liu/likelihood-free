# devtools::install_github("haochu-liu/simARG")
library(simARG)
library(jsonlite)


set.seed(100)
tree <- clonal_genealogy(15L)
ARG <- ClonalOrigin_pair_seq_fast(tree, 0.2, 100000L, 300, 50L)
num_mutations <- rep(NA, 10000)
for (i in 1:10000) {
  ARG.mutated <- FSM_mutation_fast(ARG, 0.2)
  leaf_sites <- ARG.mutated$node_site[1:15, ]
  num_mutations[i] <- sum(leaf_sites[, 1] | leaf_sites[, 2])

  if (i %% 100 == 0) {
    print(paste("Complete", i, "iteration."))
  }
}

ARG_data <- list(
  clonal_edge=tree$edge,
  clonal_node=tree$node,
  ARG_edge=ARG$edge,
  ARG_edge_mat=ARG$edge_mat,
  ARG_node_height=ARG$node_height,
  ARG_node_mat=ARG$node_mat,
  ARG_node_clonal=ARG$node_clonal,
  ARG_sum_time=ARG$sum_time
)

write.csv(num_mutations, file="pysimARG/test_data/num_mutations.csv",
          row.names=FALSE, col.names=FALSE)
write_json(ARG_data, "pysimARG/test_data/ARG_data.json",
           pretty = TRUE, auto_unbox = TRUE)
