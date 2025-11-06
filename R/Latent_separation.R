#' LP-based latent separation detection (margin-max LP; quasi via non-strict feasibility)
#'
#' @param y Binary outcome (0/1 or 2-level factor/character/logical).
#' @param X Matrix or data frame of predictors (columns = variables).
#' @param test_combinations If TRUE, test every subset of size ≥ min_vars (legacy exhaustive mode).
#' @param min_vars Minimum subset size to consider.
#' @param epsilon Ignored in margin-max LP; kept for API compatibility.
#' @param only_perfect Legacy flag for exhaustive mode (kept for compatibility).
#' @param find_minimal If TRUE, return inclusion-minimal separating subsets using pruning.
#' @param mode "either", "perfect", or "quasi" – what to count as a hit.
#' @param max_vars Optional upper bound on subset size (default: all the way to p).
#' @param stop_at_first If TRUE and find_minimal=TRUE, stop after the first minimal hit.
#' @param missing "complete" or "impute" for predictors (never impute outcome).
#' @param impute_args Options for imputation when missing="impute".
#' @param scale_X If TRUE, scale encoded predictors with base scale() (recommended).
#' @return Single-result list or a list with $minimal_subsets (same API as before).
#' @export
latent_separation <- function(
    y, X,
    test_combinations = FALSE,
    min_vars = 2,
    epsilon = 1e-5,
    only_perfect = FALSE,
    find_minimal = FALSE,
    mode = c("either","perfect","quasi"),
    max_vars = NULL,
    stop_at_first = FALSE,
    missing = c("complete","impute"),
    impute_args = list(),
    scale_X = FALSE
) {
  mode <- match.arg(mode)
  missing <- match.arg(missing)

  # ---- coerce X to data.frame ----
  if (is.matrix(X)) {
    X_raw <- as.data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (!is.data.frame(X)) {
    X_raw <- data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    X_raw <- data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  }
  var_names <- colnames(X_raw)
  p <- ncol(X_raw); if (p == 0L) stop("X has no columns.")

  if (is.null(max_vars)) max_vars <- p
  min_vars <- max(1L, as.integer(min_vars))
  max_vars <- max(min_vars, as.integer(max_vars))

  y_raw <- y  # keep for per-subset missing handling

  # ---------- helpers ----------
  run_one <- function(cols, mode_local = mode) {
    mh <- .handle_missing_for_subset(y = y_raw, X = X_raw[, cols, drop = FALSE],
                                     method = missing, impute_args = impute_args)

    y1 <- encode_outcome(mh$y)           # expect 0/1 numeric (or logical)
    if (length(unique(y1)) < 2L) {
      return(list(
        hit = FALSE, type = "no separation problem", removed = NULL,
        missing_info = list(method = mh$params_used$method, params = mh$params_used,
                            rows_used = mh$rows_used, n_used = length(mh$rows_used)),
        message = "Fewer than two outcome classes after missing-data handling."
      ))
    }

    X1 <- mh$X
    X1_mm <- encode_predictors_lp(X1)    # numeric, no NAs
    if (isTRUE(scale_X)) X1_mm <- scale(X1_mm)

    # margin-max LP (Step 1) + non-strict feasibility (Step 2)
    stat <- separation_status_LP(y = y1, X = X1_mm)

    is_hit <- switch(mode_local,
                     "perfect" = identical(stat$type, "perfect separation"),
                     "quasi"   = identical(stat$type, "quasi-complete separation"),
                     "either"  = stat$type %in% c("perfect separation", "quasi-complete separation")
    )

    if (!is_hit) {
      return(list(
        hit = FALSE, type = "no separation problem", removed = NULL,
        diagnostics = stat,
        missing_info = list(method = mh$params_used$method, params = mh$params_used,
                            rows_used = mh$rows_used, n_used = length(mh$rows_used))
      ))
    }

    list(
      hit = TRUE,
      type = stat$type,
      removed = NULL,                          # no row-deletion needed with this test
      diagnostics = stat,
      missing_info = list(method = mh$params_used$method, params = mh$params_used,
                          rows_used = mh$rows_used, n_used = length(mh$rows_used))
    )
  }

  if (!test_combinations && !find_minimal) {
    cols_all <- seq_len(p)
    res <- run_one(cols_all, if (only_perfect) "perfect" else mode)

    if (!isTRUE(res$hit)) {
      available <- character(0)
      if (identical(res$diagnostics$type, "perfect separation")) available <- c(available, "perfect")
      if (identical(res$diagnostics$type, "quasi-complete separation")) available <- c(available, "quasi")

      if (length(available)) {
        return(list(
          type = "no separation (requested mode unsatisfied)",
          satisfied = FALSE,
          available_types = available,
          removed = NULL,
          message = sprintf(
            "Separation exists (%s), but it does not satisfy mode = \"%s\".",
            paste(available, collapse = " & "),
            if (isTRUE(only_perfect)) "perfect" else mode
          ),
          missing_info = res$missing_info
        ))
      }

      return(list(
        type = "no separation problem",
        satisfied = FALSE,
        available_types = character(0),
        message = "No separation problem detected.",
        missing_info = res$missing_info
      ))
    }

    return(list(
      type = res$type,
      satisfied = TRUE,
      available_types = if (identical(res$type, "perfect separation")) c("perfect") else c("quasi"),
      removed = res$removed,
      message = if (identical(res$type, "perfect separation"))
        "Perfect separation found on the given predictors."
      else
        "Quasi-complete separation detected (non-strict feasible; margin is zero).",
      missing_info = res$missing_info
    ))
  }

  # legacy exhaustive mode
  if (test_combinations && !find_minimal) {
    results <- list()
    for (k in seq.int(min_vars, p)) {
      for (cols in combn(p, k, simplify = FALSE)) {
        nm <- paste(var_names[cols], collapse = "_")
        out <- run_one(cols, if (only_perfect) "perfect" else mode)
        if (out$hit) {
          results[[nm]] <- list(
            type = out$type,
            vars = var_names[cols],
            removed = out$removed,
            missing_info = out$missing_info
          )
        }
      }
    }
    if (only_perfect) {
      results <- Filter(function(z) identical(z$type, "perfect separation"), results)
    }
    return(results)
  }

  # minimal-subset search with pruning
  minimal <- list()
  minimal_idx <- list()
  contains_any_minimal <- function(cols) {
    if (!length(minimal_idx)) return(FALSE)
    for (m in minimal_idx) if (all(m %in% cols)) return(TRUE)
    FALSE
  }

  for (k in seq.int(min_vars, max_vars)) {
    cand <- combn(p, k, simplify = FALSE)
    cand <- Filter(function(cols) !contains_any_minimal(cols), cand)
    if (!length(cand)) next

    for (cols in cand) {
      out <- run_one(cols, mode)
      if (out$hit) {
        nm <- paste(var_names[cols], collapse = "_")
        minimal[[nm]] <- list(
          type = out$type,
          vars = var_names[cols],
          idx  = cols,
          removed = out$removed,
          missing_info = out$missing_info
        )
        minimal_idx[[length(minimal_idx) + 1L]] <- cols
        if (isTRUE(stop_at_first)) return(list(minimal_subsets = minimal))
      }
    }
  }

  list(minimal_subsets = minimal)
}

# -------------------------------------------------------------------
# INTERNALS: margin-max LP + quasi feasibility
# -------------------------------------------------------------------

# Map y in {0,1}/logical/factor(2) -> {-1,+1}
.y_to_pm1 <- function(y) {
  if (is.logical(y)) return(ifelse(y, 1L, -1L))
  if (is.numeric(y)) {
    y01 <- as.integer(round(y))
    if (!all(y01 %in% c(0L,1L))) stop("Outcome must be binary 0/1 if numeric.")
    return(ifelse(y01 == 1L, 1L, -1L))
  }
  if (is.factor(y) || is.character(y)) {
    f <- factor(y)
    if (nlevels(f) != 2L) stop("Outcome must have exactly two levels.")
    return(ifelse(as.integer(f) == 2L, 1L, -1L))  # second level -> +1
  }
  stop("Unsupported outcome type.")
}

# Solve: maximize epsilon s.t. y*(alpha + X beta) >= epsilon, ||beta||_1 <= 1
# Variables (all >=0 by lpSolve):
#   [Apos, Aneg, Bpos_1..p, Bneg_1..p, EPS]
# alpha = Apos - Aneg ; beta_j = Bpos_j - Bneg_j
.margin_max_lp <- function(ysign, X) {
  n <- nrow(X); p <- ncol(X)
  nvar <- 2 + 2*p + 1
  idx_eps <- nvar

  f.obj <- rep(0, nvar)
  f.obj[idx_eps] <- 1  # maximize EPS

  # Row-wise constraints: y_i*(alpha + x_i^T beta) - EPS >= 0
  A <- matrix(0, nrow = n + 1, ncol = nvar)

  # Apos/Aneg
  A[1:n, 1] <-  ysign
  A[1:n, 2] <- -ysign

  # Bpos/Bneg
  A[1:n, 3:(2+p)]         <-  sweep(X, 1, ysign, `*`)
  A[1:n, (3+p):(2+2*p)]   <- -sweep(X, 1, ysign, `*`)

  # -EPS
  A[1:n, idx_eps] <- -1

  # L1 norm constraint on beta: sum(Bpos + Bneg) <= 1
  A[n+1, 3:(2+p)]       <- 1
  A[n+1, (3+p):(2+2*p)] <- 1

  dir <- c(rep(">=", n), "<=")
  rhs <- c(rep(0, n), 1)

  out <- lpSolve::lp("max", f.obj, A, dir, rhs)
  list(status = out$status, eps = if (out$status == 0L) out$objval else NA_real_)
}

# Feasibility for quasi: y*(alpha + X beta) >= 0 with ||beta||_1 = 1
._quasi_feasible_lp <- function(ysign, X) {
  n <- nrow(X); p <- ncol(X)
  nvar <- 2 + 2*p  # no EPS here

  f.obj <- rep(0, nvar)  # feasibility only

  A <- matrix(0, nrow = n + 2, ncol = nvar)
  # rows: margins >= 0
  A[1:n, 1] <-  ysign
  A[1:n, 2] <- -ysign
  A[1:n, 3:(2+p)]       <-  sweep(X, 1, ysign, `*`)
  A[1:n, (3+p):(2+2*p)] <- -sweep(X, 1, ysign, `*`)

  # ||beta||_1 = 1  (two inequalities to enforce equality)
  A[n+1, 3:(2+p)]       <- 1
  A[n+1, (3+p):(2+2*p)] <- 1
  A[n+2, 3:(2+p)]       <- 1
  A[n+2, (3+p):(2+2*p)] <- 1

  dir <- c(rep(">=", n), ">=", "<=")
  rhs <- c(rep(0, n), 1, 1)

  out <- lpSolve::lp("min", f.obj, A, dir, rhs)
  out$status == 0L
}

# Main status function used by run_one()
separation_status_LP <- function(y, X) {
  # y expected 0/1 or logical/factor 2-level; convert to {-1,+1}
  ysign <- .y_to_pm1(y)
  X <- as.matrix(X); storage.mode(X) <- "double"
  if (any(!is.finite(X))) stop("X must be finite.")

  # STEP 1: margin-max
  mm <- .margin_max_lp(ysign, X)
  if (mm$status == 0L && is.finite(mm$eps) && mm$eps > 0) {
    return(list(type = "perfect separation", eps = mm$eps))
  }

  # STEP 2: non-strict feasibility (quasi)
  if (._quasi_feasible_lp(ysign, X)) {
    return(list(type = "quasi-complete separation", eps = 0))
  }

  list(type = "no separation problem", eps = 0)
}
