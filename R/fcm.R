new_fcm <- function(centers, U, obj, call) {
    structure(list(centers = centers, U = U, obj = obj, call = call),
              class = "fcm")
}

print.fcm <- function(x, ...) {
    nc      <- nrow(x$centers)
    nd      <- ncol(x$centers)
    has_inf <- ncol(x$U) == nc + 1L
    cat(sprintf("fcm: %d cluster%s, %d feature%s, %d iteration%s%s\n",
                nc, if (nc == 1L) "" else "s",
                nd, if (nd == 1L) "" else "s",
                length(x$obj), if (length(x$obj) == 1L) "" else "s",
                if (has_inf) " (with cluster at infinity)" else ""))
    invisible(x)
}

summary.fcm <- function(object, ...) {
    stop("summary.fcm is not yet implemented")
}

plot.fcm <- function(x, ...) {
    stop("plot.fcm is not yet implemented")
}

update.fcm <- function(object, ...) {
    call   <- object$call
    extras <- match.call(expand.dots = FALSE)$`...`
    for (nm in names(extras)) call[[nm]] <- extras[[nm]]
    eval(call, parent.frame())
}

new_rho_violation <- function(rho, point, nearest_dist) {
    structure(
        class = c("rho_violation", "error", "condition"),
        list(
            message = sprintf(
                "point %d has nearest-center distance %.6g > rho=%.6g",
                point, nearest_dist, rho
            ),
            rho          = rho,
            point        = point,
            nearest_dist = nearest_dist
        )
    )
}

#' Fuzzy c-means clustering
#'
#' @param X numeric matrix, n x d (samples x features)
#' @param k integer >= 2, number of clusters
#' @param rho NULL for standard FCM, or a positive scalar to add a cluster at
#'   infinity for outlier handling
#' @param m numeric > 1, fuzzifier exponent (default 2.0)
#' @param tol convergence tolerance on the membership matrix (default 1e-6)
#' @param maxiter maximum number of iterations (default 300)
#' @param seed integer seed for reproducibility, or NULL
#' @return An object of class \code{fcm} with components \code{centers} (k x d),
#'   \code{U} (n x k, or n x (k+1) when \code{rho} is supplied), \code{obj}
#'   (objective value per iteration), and \code{call}.
#' @export
fcm <- function(X, k, rho = NULL, m = 2.0, tol = 1e-6, maxiter = 300L,
                seed = NULL) {
    cl <- match.call()

    if (!is.numeric(X) || !is.matrix(X))
        stop("X must be a numeric matrix")
    if (!is.numeric(m) || length(m) != 1L || m <= 1.0)
        stop(sprintf("m must be > 1, got %g", m))
    k <- as.integer(k)
    if (length(k) != 1L || is.na(k) || k < 2L)
        stop(sprintf("k must be >= 2, got %d", k))
    n <- nrow(X)
    d <- ncol(X)
    if (n < k)
        stop(sprintf("n=%d samples must be >= k=%d clusters", n, k))
    if (!is.null(rho)) {
        if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || rho <= 0)
            stop(sprintf("rho must be a positive scalar, got %g",
                         if (is.numeric(rho) && length(rho) == 1L) rho else NaN))
    }
    maxiter <- as.integer(maxiter)

    if (!is.null(seed)) set.seed(seed)
    if (!is.double(X)) X <- matrix(as.double(X), nrow = n, ncol = d)

    if (is.null(rho)) {
        raw <- .Call(fcm_standard_c, X, k, as.double(m),
                     as.double(tol), maxiter)
    } else {
        raw <- .Call(fcm_infinity_c, X, k, as.double(m),
                     as.double(tol), maxiter, as.double(rho))
        if (identical(raw[["type"]], "rho_violation"))
            stop(new_rho_violation(raw[["rho"]], raw[["point"]],
                                   raw[["nearest_dist"]]))
    }

    new_fcm(raw[["centers"]], raw[["U"]], raw[["obj"]], cl)
}
