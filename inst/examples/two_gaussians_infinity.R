# Gaussian clusters with a cluster at infinity
#
# Rscript examples/two_gaussians_infinity.R

library(bezdek)

set.seed(42)
X <- rbind(
    matrix(rnorm(300) + rep(c(-1.5, 0.0), each = 150), 150, 2),
    matrix(rnorm(300) + rep(c( 1.5, 0.0), each = 150), 150, 2),
    matrix(rnorm(10) * 5.0,                               5, 2)
)

r <- fcm(X, 2, rho = 15.0, seed = 42L)

cat("True centers:  [-1.5, 0.0]  and  [1.5, 0.0]\n")
for (k in seq_len(nrow(r$centers)))
    cat(sprintf("FCM  center %d:  [%.3f, %.3f]\n", k, r$centers[k, 1], r$centers[k, 2]))
cat(sprintf("Converged in %d iterations\n", length(r$obj)))

U_inf <- r$U[, ncol(r$U)]
blob  <- U_inf[1:300]
noise <- U_inf[301:305]
cat(sprintf("Inf membership — blobs (mean/max): %.4f / %.4f\n", mean(blob), max(blob)))
cat(sprintf("Inf membership — noise (mean/max): %.4f / %.4f\n", mean(noise), max(noise)))
