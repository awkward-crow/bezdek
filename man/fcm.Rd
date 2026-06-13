\name{fcm}
\alias{fcm}
\alias{print.fcm}
\alias{update.fcm}
\title{Fuzzy C-Means Clustering}
\description{
  Fuzzy c-means clustering of a numeric data matrix into \code{c} clusters,
  with an optional cluster at infinity for outlier handling.
}
\usage{
fcm(X, k, rho = NULL, m = 2.0, tol = 1e-6, maxiter = 300L, seed = NULL)
}
\arguments{
  \item{X}{numeric matrix, \eqn{n \times d} (samples \eqn{\times} features).}
  \item{k}{integer \eqn{\ge 2}, number of clusters.}
  \item{rho}{\code{NULL} for standard FCM, or a positive scalar to add a
    cluster at infinity whose distance from point \eqn{x_i} is
    \eqn{\rho - \min_k \|x_i - v_k\|}. Must exceed every point-to-center
    distance; a \code{rho_violation} error is raised if it does not.}
  \item{m}{fuzzifier exponent, numeric \eqn{> 1} (default \code{2.0}).
    Larger values produce softer memberships.}
  \item{tol}{convergence tolerance on the membership matrix (default \code{1e-6}).}
  \item{maxiter}{maximum number of iterations (default \code{300}).}
  \item{seed}{integer seed passed to \code{\link{set.seed}} before the
    algorithm runs, or \code{NULL} to use the current RNG state.}
}
\value{
  An object of class \code{"fcm"} with components:
  \describe{
    \item{\code{centers}}{numeric matrix, \eqn{k \times d}, cluster centers.}
    \item{\code{U}}{numeric matrix, \eqn{n \times k} (standard) or
      \eqn{n \times (k+1)} (with \code{rho}); last column is the
      infinity-cluster membership when \code{rho} is supplied.}
    \item{\code{obj}}{numeric vector of objective values, one per iteration.}
    \item{\code{call}}{the matched call.}
  }
}
\section{Errors}{
  Argument errors (invalid \code{m}, \code{c}, \code{rho}, or incompatible
  dimensions) are raised as plain errors.

  A \code{rho_violation} condition (class \code{c("rho_violation", "error",
  "condition")}) is raised when the nearest-center distance of some point
  exceeds \code{rho}. The condition has fields \code{rho}, \code{point}
  (1-based index), and \code{nearest_dist}.
}
\seealso{
  \code{\link{update}} for re-running with modified arguments.
}
\examples{
set.seed(1)
X <- matrix(rnorm(200), 100, 2)

# standard FCM
r <- fcm(X, 3, seed = 1L)
print(r)

# cluster at infinity (handles outliers)
r2 <- fcm(X, 3, rho = 20.0, seed = 1L)

# re-run with a different fuzzifier
r3 <- update(r, m = 1.5)

# catch rho_violation
tryCatch(
  fcm(rbind(X, c(1000, 1000)), 3, rho = 1.0, seed = 1L),
  rho_violation = function(e) message("violation at point ", e$point)
)
}
\keyword{cluster}
