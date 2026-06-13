# Seeds dataset: geometric measurements of wheat kernels (UCI repository).
# Charytanowicz et al. (2010), "A Complete Gradient Clustering Algorithm for
# Features Analysis of X-ray Images", Information Technologies in Biomedicine,
# Springer, pp. 15-24.
#
# 210 samples, 7 features, 3 varieties: Kama (1), Rosa (2), Canadian (3).
# Kama and Rosa overlap substantially in feature space, making membership
# values near their boundary genuinely informative.
#
# Rscript examples/seeds.R [seed]

library(bezdek)

args  <- commandArgs(trailingOnly = TRUE)
seed  <- if (length(args) >= 1) as.integer(args[1]) else 42L

url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/00236/seeds_dataset.txt"
raw <- tryCatch(
    read.table(url, header = FALSE),
    error = function(e) stop("Could not download seeds dataset: ", conditionMessage(e))
)
X     <- as.matrix(raw[, 1:7])
truth <- raw[, 8]

varieties <- c("Kama", "Rosa", "Canadian")

r <- fcm(X, 3, seed = seed)

# Confusion matrix: rows = FCM cluster, cols = true variety
hard <- max.col(r$U)
C    <- table(FCM = hard, Variety = truth)
rownames(C) <- paste("FCM", 1:3)
colnames(C) <- varieties
cat("Confusion matrix (FCM cluster x true variety):\n")
print(C)

# Most ambiguous points by membership entropy
entropy_row <- function(u) {
    u <- u[u > 0]
    -sum(u * log(u))
}
ent  <- apply(r$U, 1, entropy_row)
top5 <- order(ent, decreasing = TRUE)[1:5]

cat("\nMost ambiguous points (highest membership entropy):\n")
cat(sprintf("  %3s  %-10s  %6s   %6s   %6s\n", "i", "variety", "U[1]", "U[2]", "U[3]"))
for (i in top5) {
    cat(sprintf("  %3d  %-10s  %6.3f   %6.3f   %6.3f\n",
                i, varieties[truth[i]], r$U[i, 1], r$U[i, 2], r$U[i, 3]))
}
