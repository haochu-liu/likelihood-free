# devtools::install_github("haochu-liu/simARG")
library(simARG)


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

write.csv(num_mutations, file="pysimARG/test_data/num_mutations.csv",
          row.names=FALSE, col.names=FALSE)
