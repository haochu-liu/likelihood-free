delta = 20
L = 100

P1 <- delta / (delta + L - 1)
P2 <- delta * (1-(1-1/delta)^L) / L

s <- 0
n <- 1000000
for (i in 1:n) {
  start_pos <- sample(1:L, 1, replace=TRUE)
  tract_length <- rgeom(1, prob = 1/delta) + 1

  start_to_first <- L - start_pos + 1
  if (start_to_first <= tract_length) {
    s <- s + 1
  }
}
est_P <- s/n
print(est_P)
print(abs(est_P - P1))
print(abs(est_P - P2))
