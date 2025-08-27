#' LP-based latent separation detection
#'
#' Check for perfect or quasi-complete separation in multivariate predictors.
#'
#' @param y Numeric vector of binary outcomes (0/1).
#' @param X Numeric matrix or data frame of predictors (columns = variables).
#' @param test_combinations Logical; if TRUE, tests every subset of variables of size ≥ \code{min_vars}.
#' @param min_vars Integer; minimum number of variables in subsets when \code{test_combinations=TRUE}.
#' @return For single X: a list with \code{type} ("perfect separation", "quasi-complete separation", or "no separation problem"), and optionally \code{removed} and \code{message}. For combinations: named list of results per subset.
#' @importFrom lpSolve lp
#' @export
latent_separation <- function(y, X, test_combinations = FALSE, min_vars = 2, only_perfect = FALSE) {
  var_names <- colnames(X)

  if (test_combinations) {
    p <- ncol(X); results <- list()

    for (k in seq(min_vars, p)) {
      for (cols in combn(p, k, simplify = FALSE)) {
        subset_name <- paste(var_names[cols], collapse = "_")
        results[[subset_name]] <- latent_separation(y, X[, cols, drop = FALSE], FALSE)
      }
    }
    if (only_perfect) {
      # keep only those with perfect separation
      results <- Filter(function(z) identical(z$type, "perfect separation"), results)
    }

    return(results)
  }
  full_res <- feasibility_check_LP(y, as.matrix(X))
  if (full_res$ge$status == 0 || full_res$le$status == 0) {
    return(list(
      type    = "perfect separation",
      message = "Perfect separation found on the full dataset."
    ))
  }
  quasi_inds <- which(vapply(seq_along(y), function(i) {
    sub <- feasibility_check_LP(y[-i], as.matrix(X[-i, , drop = FALSE]))
    sub$ge$status == 0 || sub$le$status == 0
  }, logical(1)))
  if (length(quasi_inds) > 0) {
    msg <- sprintf(
      "Quasi-complete separation: removing any single observation from {%s} individually achieves perfect separation.",
      paste(quasi_inds, collapse = ", ")
    )
    return(list(
      type    = "quasi-complete separation",
      removed = quasi_inds,
      message = msg
    ))
  }
  list(
    type    = "no separation problem",
    message = "No separation problem detected."
  )
}


# Internal LP feasibility check
#' @keywords internal
#' @noRd
feasibility_check_LP <- function(y, X) {
  p <- ncol(X); n <- length(y); eps <- 1e-5
  f.obj <- rep(0, 2*p + 1)
  A <- matrix(0, nrow = n+1, ncol = 2*p + 1)
  A[1:n, 1] <- 1
  for (i in seq_len(n)) {
    A[i, 2:(p+1)]     <- X[i, ];
    A[i, (p+2):(2*p+1)] <- -X[i, ]
  }
  A[n+1, 2:(p+1)] <- 1; A[n+1, (p+2):(2*p+1)] <- -1
  dir <- c(ifelse(y == 1, ">=", "<="), ">=")
  rhs <- c(rep(0, n), eps)
  lp1 <- lpSolve::lp("min", f.obj, A, dir, rhs)
  dir[n+1] <- "<="; rhs[n+1] <- -eps
  lp2 <- lpSolve::lp("min", f.obj, A, dir, rhs)
  list(ge = list(status = lp1$status, sol = lp1$solution),
       le = list(status = lp2$status, sol = lp2$solution))
}
