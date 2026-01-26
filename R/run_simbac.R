library(simARG)


n <- 20L
rho_site <- 10 / 1e5
delta <- 30
R_time_vec <- rep(NA, 100)
Rcpp_time_vec <- rep(NA, 100)
for (i in 1:100) {
  set.seed(i)
  
  time_result <- system.time(
    simbac_ARG(n, rho_site, 1e5L, delta, node_max = 1000)
  )
  R_time_vec[i] <- time_result["elapsed"]

  time_result <- system.time(
    simbac_ARG.decimal(n, rho_site, 1e5L, delta, node_max = 1000)
  )
  Rcpp_time_vec[i] <- time_result["elapsed"]

  print(paste("Complete", i, "iterations"))
}

print(R_time_vec)
print(Rcpp_time_vec)
simbac_time <- data.frame(
  R=R_time_vec,
  R_Rcpp=Rcpp_time_vec,
  R_in_python=NA,
  R_Rcpp_in_python=NA
)
save(simbac_time, file="data/simbac_time.RData")
