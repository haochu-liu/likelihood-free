# devtools::install_github("haochu-liu/simARG")
library(simARG)
library(ggplot2)
library(tidyr)
library(dplyr)


set.seed(100)
tree <- clonal_genealogy(15L)
rho_site <- 0.2
L <- 100000L
delta <- 300
k <- 2000L
theta_site <- 0.2


running_time <- matrix(NA, nrow=1000, ncol=4)
colnames(running_time) <- c("ClonalOrigin_pair", "ClonalOrigin_pair_fast",
                            "add_mutation", "add_mutation_fast")

for (i in 1:1000) {
  set.seed(i)
  start_time_1 <- Sys.time()
  ARG <- ClonalOrigin_pair_seq(tree, rho_site, L, delta, k)
  end_time_1 <- Sys.time()
  running_time[i, 1] <- end_time_1 - start_time_1

  set.seed(i)
  start_time_2 <- Sys.time()
  ARG <- ClonalOrigin_pair_seq_fast(tree, rho_site, L, delta, k)
  end_time_2 <- Sys.time()
  running_time[i, 2] <- end_time_2 - start_time_2

  set.seed(i)
  start_time_3 <- Sys.time()
  ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
  end_time_3 <- Sys.time()
  running_time[i, 3] <- end_time_3 - start_time_3

  set.seed(i)
  start_time_4 <- Sys.time()
  ARG_mutated <- FSM_mutation_fast(ARG, theta_site)
  end_time_4 <- Sys.time()
  running_time[i, 4] <- end_time_4 - start_time_4

  if (i %% 100 == 0) {print(paste("Complete", i, "iterations"))}
}


time_df <- as.data.frame(running_time) %>%
  pivot_longer(cols = everything(), names_to = "Func", values_to = "Time")

ggplot(time_df, aes(x = Func, y = Time)) +
  geom_boxplot() +
  labs(title = "Running time", x = "Functions", y = "Time (sec)") +
  theme_minimal()
