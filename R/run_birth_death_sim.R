library(readr)
library(simARG)
library(ggplot2)
library(tidyr)
library(dplyr)


time_vec <- rep(NA, 100)
for (i in 1:100) {
  set.seed(i)

  start_time <- Sys.time()
  result <- birth_death_sim(100L, 10)
  end_time <- Sys.time()
  time_vec[i] <- end_time - start_time

  print(paste("Complete", i, "iterations"))
}
python_time <- read_csv("data/python_bd_time.csv")

bd_time <- data.frame(
  R=time_vec,
  python=python_time$python
)
save(bd_time, file="data/bd_time.RData")

load("data/bd_time.RData")

bd_time_long <- bd_time %>%
  pivot_longer(cols = everything(), names_to = "Language", values_to = "Time")

ggplot(bd_time_long, aes(x = Language, y = Time, fill = Language)) +
  geom_boxplot() +
  labs(title = "Birth-death process running time",
       x = "Language",
       y = "Time(sec)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
