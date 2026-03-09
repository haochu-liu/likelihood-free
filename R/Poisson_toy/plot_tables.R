# library(openxlsx)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(scales)


load("data/poisson_toy_table.RData")

# out_tables <- poisson_toy_table
# out_tables$ess <- round(out_tables$ess, 0)
# out_tables$norm_ess <- round(out_tables$norm_ess, 2)
# out_tables$err_mean <- round(out_tables$err_mean, 3)
# out_tables$err_var <- round(out_tables$err_var, 3)
# out_tables$var_log <- round(out_tables$var_log, 2)
# out_tables$relative_var <- round(out_tables$relative_var, 2)
# write.xlsx(out_tables, "data/poisson_toy_origin.xlsx",
#            colNames = TRUE, rowNames = TRUE)


# Standardized norm ESS -> exact values but standardized color map
# var-log and relative var -> exact values and exact color map with red boxes

# Norm ESS
df_norm_ess <- as.data.frame(poisson_toy_table$norm_ess)
std_norm_ess <- poisson_toy_table$norm_ess
for (i in 1:ncol(std_norm_ess)) {
  std_norm_ess[, i] <- order(poisson_toy_table$norm_ess[, i])
}
df_std_norm_ess <- as.data.frame(std_norm_ess)

df_norm_ess_long <- df_norm_ess %>%
  rownames_to_column(var = "Row") %>%
  pivot_longer(cols = -Row,
               names_to = "Column",
               values_to = "norm_ess") %>%
  mutate(
    Row = as.numeric(Row),
    Column = as.numeric(Column)
  )

df_std_norm_ess_long <- df_std_norm_ess %>%
  rownames_to_column(var = "Row") %>%
  pivot_longer(cols = -Row,
               names_to = "Column",
               values_to = "std_colormap") %>%
  mutate(
    Row = as.numeric(Row),
    Column = as.numeric(Column)
  )
sum(df_std_norm_ess_long$Row == df_norm_ess_long$Row)
sum(df_std_norm_ess_long$Column == df_norm_ess_long$Column)

df_norm_ess_long$std_colormap <- df_std_norm_ess_long$std_colormap

ggplot(df_norm_ess_long, aes(x = Column, y = Row)) +
  geom_tile(aes(fill = std_colormap), color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 50) +
  # geom_text(aes(label = round(norm_ess, 2)), size = 4) +
  theme_minimal() +
  labs(title = "Heatmap of normalized ESS",
       fill = "Rank of norm ESS \nfor each dimension",
       x = "Dimension",
       y = "Number of samples") +
  theme(plot.title = element_text(face = "bold", size = 12))

# Var-log
df_var_log <- as.data.frame(poisson_toy_table$var_log)
df_var_log_long <- df_var_log %>%
  rownames_to_column(var = "Row") %>%
  pivot_longer(cols = -Row,
               names_to = "Column",
               values_to = "var_log") %>%
  mutate(
    Row = as.numeric(Row),
    Column = as.numeric(Column)
  )

df_var_log_highlight <- df_var_log_long %>%
  filter(var_log >= 1 & var_log <= 3)

ggplot(df_var_log_long, aes(x = Column, y = Row)) +
  geom_tile(aes(fill = log(var_log)), color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1,
                       limits = c(0, 10),
                       oob = scales::squish) +
  # geom_text(aes(label = round(Text_Value, 2)), size = 4) +
  geom_tile(data = df_var_log_highlight, fill = NA, color = "black", linewidth = 0.2) +
  theme_minimal() +
  labs(title = "Heatmap of var[log-likelihood]\n(values between 1 and 3 are highlighted by boxes)",
       fill = "Variance values \nthat squished to [0, 10]",
       x = "Dimension",
       y = "Number of samples") +
  theme(plot.title = element_text(face = "bold", size = 12))

# Relative-var
df_relative_var <- as.data.frame(poisson_toy_table$relative_var)
df_relative_var_long <- df_relative_var %>%
  rownames_to_column(var = "Row") %>%
  pivot_longer(cols = -Row,
               names_to = "Column",
               values_to = "relative_var") %>%
  mutate(
    Row = as.numeric(Row),
    Column = as.numeric(Column)
  )

df_relative_var_highlight <- df_relative_var_long %>%
  filter(relative_var >= 1 & relative_var <= 2)

ggplot(df_relative_var_long, aes(x = Column, y = Row)) +
  geom_tile(aes(fill = log(relative_var)), color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1,
                       limits = c(0, 10),
                       oob = scales::squish) +
  # geom_text(aes(label = round(Text_Value, 2)), size = 4) +
  geom_tile(data = df_relative_var_highlight, fill = NA, color = "black", linewidth = 0.2) +
  theme_minimal() +
  labs(title = "Heatmap of relative variance\n(values between 1 and 2 are highlighted by boxes)",
       fill = "Relative variance values \nthat squished to [0, 10]",
       x = "Dimension",
       y = "Number of samples") +
  theme(plot.title = element_text(face = "bold", size = 12))


