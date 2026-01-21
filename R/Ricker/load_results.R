load("data/ricker_50_short_1.RData")
df1 <- ricker_50

load("data/ricker_50_short_2.RData")
df2 <- ricker_50

load("data/ricker_50_short_3.RData")
df3 <- ricker_50

load("data/ricker_50_short_4.RData")
df4 <- ricker_50

load("data/ricker_50_short_5.RData")
df5 <- ricker_50

n <- c(df1$n[1], df2$n, df1$n[2], df3$n, df4$n, df5$n)
acc_rate <- c(df1$acc_rate[1], df2$acc_rate, df1$acc_rate[2],
              df3$acc_rate, df4$acc_rate, df5$acc_rate)

ess <- matrix(NA, nrow=3, ncol=7)
ess[, c(1, 4)] <- df1$ess
ess[, c(2, 3)] <- df2$ess
ess[, 5] <- df3$ess
ess[, 6] <- df4$ess
ess[, 7] <- df5$ess
colnames(ess) <- n

norm_ess <- matrix(NA, nrow=3, ncol=7)
norm_ess[, c(1, 4)] <- df1$norm_ess
norm_ess[, c(2, 3)] <- df2$norm_ess
norm_ess[, 5] <- df3$norm_ess
norm_ess[, 6] <- df4$norm_ess
norm_ess[, 7] <- df5$norm_ess
colnames(norm_ess) <- n

post_mean <- matrix(NA, nrow=3, ncol=7)
post_mean[, c(1, 4)] <- df1$post_mean
post_mean[, c(2, 3)] <- df2$post_mean
post_mean[, 5] <- df3$post_mean
post_mean[, 6] <- df4$post_mean
post_mean[, 7] <- df5$post_mean
colnames(post_mean) <- n
