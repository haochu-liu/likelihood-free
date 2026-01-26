library(readr)
library(simARG)
library(ggplot2)
library(tidyr)
library(dplyr)


n <- 20L
rho_site <- 10 / 1e5
delta <- 30
R_time_vec <- rep(NA, 100)
Rcpp_time_vec <- rep(NA, 100)
for (i in 1:100) {
  set.seed(i)

  start_time_R <- Sys.time()
  simbac_ARG(n, rho_site, 1e5L, delta, node_max = 1000)
  end_time_R <- Sys.time()
  R_time_vec[i] <- end_time_R - start_time_R

  start_time_Rcpp <- Sys.time()
  simbac_ARG.decimal(n, rho_site, 1e5L, delta, node_max = 1000)
  end_time_Rcpp <- Sys.time()
  Rcpp_time_vec[i] <- end_time_Rcpp - start_time_Rcpp

  print(paste("Complete", i, "iterations"))
}

print(R_time_vec)
print(Rcpp_time_vec)
python_simbac_time <- read_csv("data/python_simbac_time.csv")

simbac_time <- data.frame(
  R=R_time_vec,
  R_Rcpp=Rcpp_time_vec,
  R_in_python=python_simbac_time$R_in_python,
  R_Rcpp_in_python=python_simbac_time$R_Rcpp_in_python,
  clean_python=python_simbac_time$clean_time
)
save(simbac_time, file="data/simbac_time.RData")

load("data/simbac_time.RData")

simbac_time_long <- simbac_time %>%
  pivot_longer(cols = everything(), names_to = "Func", values_to = "Time") %>%
  mutate(Language = ifelse(Func %in% c("R", "R_Rcpp"), "R", "R_in_python"))

ggplot(simbac_time_long, aes(x = Func, y = Time, fill = Language)) +
  geom_boxplot() +
  labs(title = "Simbac function running time",
       x = "Functions",
       y = "Time(sec)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
