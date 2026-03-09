library(openxlsx)


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
