# bezdek

[![R-CMD-check](https://github.com/awkward-crow/bezdek/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/awkward-crow/bezdek/actions/workflows/R-CMD-check.yaml)

Fuzzy c-means clustering (Bezdek 1981) with a *cluster at
infinity* to identify outliers.

## Usage

### Standard FCM

FCM partitions *n* points into *k* clusters by minimising

$$J_m = \sum_i \sum_k u_{ik}^m \|x_i - v_k\|^2$$

over the membership matrix **U** ($u_{ik} \in [0,1]$, $\sum_k u_{ik} = 1$) and
the cluster centres **V**. The parameter $m > 1$ controls softness: $m \to 1$
recovers hard k-means; larger $m$ spreads memberships more evenly.

```r
library(bezdek)

set.seed(1)
X <- rbind(
    matrix(rnorm(300) + rep(c(-3, 0), each = 150), 150, 2),
    matrix(rnorm(300) + rep(c( 3, 0), each = 150), 150, 2)
)

r <- fcm(X, k = 2, seed = 1L)
r
#> fcm: 2 clusters, 2 features, 23 iterations

r$centers   # k × d matrix of cluster centres
r$U         # n × k membership matrix (rows sum to 1)
r$obj       # objective value at each iteration
```

### Cluster at infinity

Supply `rho` to add a *k*+1-th cluster whose centre can be thought of as lying at
infinity. Its distance to point $x_i$ is defined as

$$d_\infty(x_i) = \rho - \min_k \|x_i - v_k\|$$

where `rho` is chosen to exceed every point's distance to its nearest centre. If it does
not, a `rho_violation` condition is signalled. Points far from all centres have
a small $d_\infty$, which drives their membership of the cluster at infinity toward 1. The centre
update uses only the regular memberships (normalised to sum to one), so
outliers exert little influence on the centres.

The last column of `$U` holds the membership of the cluster at infinity.

```r
X_noisy <- rbind(X, matrix(c(0, 20), 1, 2))   # one extreme outlier

r <- fcm(X_noisy, k = 2, rho = 25, seed = 1L)
r
#> fcm: 2 clusters, 2 features, 23 iterations (with cluster at infinity)

r$U[nrow(r$U), ]   # outlier has near-unit membership in the last column
```

If `rho` is too small:

```r
tryCatch(
    fcm(X_noisy, k = 2, rho = 1),
    rho_violation = function(e) cat("Violated at point", e$point, "\n")
)
```

## Arguments

| Argument  | Default | Description |
|-----------|---------|-------------|
| `X`       | —       | Numeric matrix, *n* × *d* |
| `k`       | —       | Number of clusters (≥ 2) |
| `rho`     | `NULL`  | Scale for the infinity cluster; `NULL` for standard FCM |
| `m`       | `2.0`   | Fuzzifier exponent (> 1) |
| `tol`     | `1e-6`  | Convergence tolerance on the membership matrix |
| `maxiter` | `300`   | Maximum iterations |
| `seed`    | `NULL`  | RNG seed for reproducibility |

## References

Bezdek, J. C. (1981). *Pattern Recognition with Fuzzy Objective Function
Algorithms*. Plenum Press.

Davé, R. N. (1991). Characterization and detection of noise in clustering.
*Pattern Recognition Letters*, 12(11), 657–664. *(cluster at infinity)*
