library(simARG)


set.seed(100)
tree <- clonal_genealogy(15L)
R_time <- rep(NA, 100)
for (i in 1:100) {
  set.seed(i)

  start_time <- Sys.time()
  result <- ClonalOrigin_pair_seq.simulator(
    tree,
    0.02,
    0.05,
    100000L,
    300,
    100,
    c(50L, 200L, 2000L)
  )
  end_time <- Sys.time()
  R_time[i] <- end_time - start_time

  print(paste("Complete", i, "iteration."))
}

write.csv(R_time, file="data/ClonalOrigin_time.csv",
          row.names=FALSE, col.names=FALSE)
