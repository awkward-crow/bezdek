test_that("argument validation: m and k", {
    X <- matrix(runif(20), 10, 2)
    expect_error(fcm(X, 2, m = 1.0),               "m must be > 1")
    expect_error(fcm(X, 2, m = 0.5),               "m must be > 1")
    expect_error(fcm(X, 1),                         "k must be >= 2")
    expect_error(fcm(matrix(runif(6), 3, 2), 5),    "samples must be >= k")
})

test_that("output shapes: standard FCM", {
    r <- fcm(matrix(runif(200), 50, 4), 3, seed = 1L)
    expect_equal(dim(r$centers), c(3L, 4L))
    expect_equal(dim(r$U),       c(50L, 3L))
})

test_that("membership in [0,1] and rows sum to 1: standard FCM", {
    r <- fcm(matrix(runif(240), 80, 3), 4, seed = 0L)
    expect_true(all(r$U >= 0 & r$U <= 1))
    expect_true(all(abs(rowSums(r$U) - 1) < 1e-10))
})

test_that("objective is non-increasing", {
    r <- fcm(matrix(runif(200), 100, 2), 3, seed = 2L)
    expect_true(all(diff(r$obj) <= 1e-10))
})

test_that("well-separated clusters yield near-hard memberships", {
    X <- rbind(matrix(rnorm(200),      100, 2),
               matrix(rnorm(200) + 30, 100, 2))
    r <- fcm(X, 2, seed = 7L)
    expect_true(all(apply(r$U, 1, max) >= 0.99))
})

test_that("reproducibility with same seed", {
    X  <- matrix(runif(200), 50, 4)
    r1 <- fcm(X, 3, seed = 42L)
    r2 <- fcm(X, 3, seed = 42L)
    expect_identical(r1$U,       r2$U)
    expect_identical(r1$centers, r2$centers)
    expect_identical(r1$obj,     r2$obj)
})

test_that("rho: argument validation", {
    X <- matrix(runif(20), 10, 2)
    expect_error(fcm(X, 2, rho = 0.0),              "rho must be a positive scalar")
    expect_error(fcm(X, 2, rho = -1.0),             "rho must be a positive scalar")
    expect_error(fcm(X, 2, rho = 1e6, m = 1.0),    "m must be > 1")
    expect_error(fcm(X, 2, rho = 1e6, m = 0.5),    "m must be > 1")
    expect_error(fcm(X, 1, rho = 1e6),              "k must be >= 2")
    expect_error(fcm(matrix(runif(6), 3, 2), 5, rho = 1e6), "samples must be >= k")
})

test_that("rho: output shapes", {
    r <- fcm(matrix(runif(200), 50, 4), 3, rho = 1e6, seed = 1L)
    expect_equal(dim(r$centers), c(3L, 4L))
    expect_equal(dim(r$U),       c(50L, 4L))
})

test_that("rho: membership in [0,1] and rows sum to 1", {
    r <- fcm(matrix(runif(240), 80, 3), 4, rho = 1e6, seed = 0L)
    expect_true(all(r$U >= 0 & r$U <= 1))
    expect_true(all(abs(rowSums(r$U) - 1) < 1e-10))
})

test_that("rho: RhoViolation is raised with correct fields", {
    set.seed(42)
    X <- rbind(
        matrix(rnorm(60) + rep(c(-1.5, 0), each = 30), 30, 2),
        matrix(rnorm(60) + rep(c( 1.5, 0), each = 30), 30, 2),
        matrix(c(20.0, 0.0), 1, 2)
    )
    err <- tryCatch(
        fcm(X, 2, rho = 5.0, seed = 42L),
        rho_violation = function(e) e
    )
    expect_s3_class(err, "rho_violation")
    expect_equal(err$rho, 5.0)
    expect_gt(err$nearest_dist, 5.0)
    expect_true(err$point >= 1L && err$point <= nrow(X))
})

test_that("rho: large rho converges close to standard FCM", {
    X     <- matrix(runif(120), 60, 2)
    r_std <- fcm(X, 3,       seed = 3L)
    r_inf <- fcm(X, 3, rho = 1e6, seed = 3L)
    expect_equal(r_inf$U[, 1:3], r_std$U,      tolerance = 1e-4)
    expect_equal(r_inf$centers,  r_std$centers, tolerance = 1e-4)
    expect_true(all(r_inf$U[, 4] < 1e-4))
})

test_that("rho: centers unaffected by extreme outlier", {
    blob1 <- matrix(rnorm(200) + rep(c(-10, 0), each = 100), 100, 2)
    blob2 <- matrix(rnorm(200) + rep(c( 10, 0), each = 100), 100, 2)
    X     <- rbind(blob1, blob2, matrix(c(0.0, 20.0), 1, 2))
    r     <- fcm(X, 2, rho = 25.0, seed = 7L)
    cx    <- r$centers[order(r$centers[, 1]), ]
    expect_equal(cx[1, ], c(-10.0, 0.0), tolerance = 0.3)
    expect_equal(cx[2, ], c( 10.0, 0.0), tolerance = 0.3)
    expect_gt(r$U[201, ncol(r$U)], 0.9)
})
