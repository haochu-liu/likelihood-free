library(simARG)


n <- 20L
rho_site <- 10 / 1e5
delta <- 30
time_vec <- rep(NA, 100)
for (i in 1:100) {
  set.seed(i)
  time_result <- system.time(
    simbac_ARG(n, rho_site, 1e5L, delta, node_max = 1000)
  )
  time_vec[i] <- time_result["elapsed"]

  print(paste("Complete", i, "iterations"))
}

print(time_vec)
save(time_vec, file="data/R_simbac_time.RData")
