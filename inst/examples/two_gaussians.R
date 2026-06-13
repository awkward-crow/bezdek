# Two Gaussian clusters
#
# Rscript examples/two_gaussians.R

library(bezdek)

set.seed(10337142)

n = 150
X <- rbind(
    sweep(matrix(rnorm(2 * n), nrow = n), 2, c(-1.5, 0), "+"),
    sweep(matrix(rnorm(2 * n), nrow = n), 2, c(1.5, 0), "+")
)

r <- fcm(X, 2)

cat("True centers:  c(-1.5, 0.0)  and  c(1.5, 0.0)\n")
for (k in seq_len(nrow(r$centers)))
    cat(sprintf("FCM  center %d:  [%.3f, %.3f]\n", k, r$centers[k, 1], r$centers[k, 2]))
cat(sprintf("Converged in %d iterations\n", length(r$obj)))

entropy_row <- function(u) -sum(u * log(u + 1e-16))
ent <- apply(r$U, 1, entropy_row)
cat(sprintf("Max membership entropy: %.4f\n", max(ent)))
