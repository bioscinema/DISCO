#' LP-based latent separation detection (with optional minimal-subset search)
#'
#' @param y Binary outcome (0/1 or 2-level factor/character/logical).
#' @param X Matrix or data frame of predictors (columns = variables).
#' @param test_combinations If TRUE, test every subset of size ≥ min_vars (legacy exhaustive mode).
#' @param min_vars Minimum subset size to consider.
#' @param epsilon Numeric scalar > 0. Margin used in LP feasibility check (default 1e-5). Larger values make the test stricter; smaller values approach the limiting definition.
#' @param only_perfect Legacy flag for exhaustive mode (kept for compatibility).
#' @param find_minimal If TRUE, return inclusion-minimal separating subsets using pruning.
#' @param mode "either", "perfect", or "quasi" – what to count as a hit.
#' @param max_vars Optional upper bound on subset size (default: all the way to p).
#' @param stop_at_first If TRUE and find_minimal=TRUE, stop after the first minimal hit.
#' @param missing How to treat missing data: "complete" (drop rows w/ NA in the subset) or "impute" (impute predictors; never impute outcome).
#' @param impute_args Optional list of imputation settings when
#'   `missing = "impute"`. Recognized keys:
#'   \itemize{
#'     \item \code{numeric_method}: "median" (default) or "mean"
#'     \item \code{categorical_method}: "mode" (default) or "missing"
#'     \item \code{logical_method}: "mode" (default) or "missing"
#'     \item \code{custom_fn}: function(data.frame) -> imputed data.frame
#'   }
#' @param scale_X Logical; if TRUE, standardize encoded predictors with base \code{scale()} (default FALSE).
#' @return If not searching: single result list {type, removed?, message, missing_info}.
#'         If searching: list with $minimal_subsets (named by vars joined with "_"), each containing
#'         {type, vars, removed?, missing_info}.
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

  # coerce X to data.frame early (we will subset columns and then encode)
  if (is.matrix(X)) {
    X_raw <- as.data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (!is.data.frame(X)) {
    X_raw <- data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    X_raw <- data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  }
  var_names <- colnames(X_raw)
  p <- ncol(X_raw)
  if (p == 0L) stop("X has no columns.")

  if (is.null(max_vars)) max_vars <- p
  min_vars <- max(1L, as.integer(min_vars))
  max_vars <- max(min_vars, as.integer(max_vars))

  y_raw <- y  # keep original for per-subset missing handling

  # helper: run a single check on a given column index set, with missing-handling + encoding
  run_one <- function(cols, mode_local = mode) {
    # missing handling on this subset only
    mh <- .handle_missing_for_subset(y = y_raw, X = X_raw[, cols, drop = FALSE],
                                     method = missing, impute_args = impute_args)

    # after missing handling, encode outcome and predictors for LP
    y1 <- encode_outcome(mh$y)
    if (length(unique(y1)) < 2L) {
      return(list(
        hit = FALSE, type = "no separation problem",
        removed = NULL,
        missing_info = list(method = mh$params_used$method, params = mh$params_used,
                            rows_used = mh$rows_used, n_used = length(mh$rows_used)),
        message = "Fewer than two outcome classes after missing-data handling."
      ))
    }

    X1 <- mh$X
    # encode_predictors_lp expects no NAs; ensured by missing handling
    X1_mm <- encode_predictors_lp(X1)

    # optional scaling with base scale()
    if (isTRUE(scale_X)) {
      X1_mm <- scale(X1_mm)
    }

    # Evaluate both, then filter by mode (strict quasi semantics)
    chk <- .check_sep_both(y1, X1_mm, epsilon = epsilon)

    is_hit <- switch(mode_local,
                     "perfect" = chk$perfect,
                     "quasi"   = (chk$quasi && !chk$perfect),  # STRICT quasi-only
                     "either"  = (chk$perfect || chk$quasi)
    )

    if (!is_hit) {
      return(list(
        hit = FALSE,
        type = "no separation problem",
        removed = NULL,
        diagnostics = chk,
        missing_info = list(
          method = mh$params_used$method,
          params = mh$params_used,
          rows_used = mh$rows_used,
          n_used = length(mh$rows_used)
        )
      ))
    }

    out_type <- if (chk$perfect) "perfect separation" else "quasi-complete separation"
    list(
      hit = TRUE,
      type = out_type,
      removed = if (chk$perfect) NULL else chk$removed,
      diagnostics = chk,
      missing_info = list(
        method = mh$params_used$method,
        params = mh$params_used,
        rows_used = mh$rows_used,
        n_used = length(mh$rows_used)
      )
    )
  }

  #   out <- .check_sep(y1, X1_mm, mode = mode_local)
  #   out$missing_info <- list(
  #     method = mh$params_used$method,
  #     params = mh$params_used,
  #     rows_used = mh$rows_used,
  #     n_used = length(mh$rows_used)
  #   )
  #   out
  # }

  # === No search: single set check (back-compat path) ===========================
  if (!test_combinations && !find_minimal) {
    cols_all <- seq_len(p)
    res <- run_one(cols_all, if (only_perfect) "perfect" else mode)

    # If requested mode is unsatisfied but some other type exists, say so
    if (!isTRUE(res$hit)) {
      available <- character(0)
      if (isTRUE(res$diagnostics$perfect)) available <- c(available, "perfect")
      if (isTRUE(res$diagnostics$quasi))   available <- c(available, "quasi")

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

    # Hit: requested mode satisfied
    return(list(
      type = res$type,
      satisfied = TRUE,
      available_types = if (identical(res$type, "perfect separation"))
        c("perfect") else c("quasi"),
      removed = res$removed,
      message = if (identical(res$type, "perfect separation"))
        "Perfect separation found on the given predictors."
      else
        sprintf("Quasi-complete separation: removing any of {%s} yields perfect.",
                paste(res$removed, collapse = ", ")),
      missing_info = res$missing_info
    ))
  }

  # === Legacy exhaustive enumeration (kept) =====================================
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

  # === Minimal-subset search with pruning (recommended) =========================
  minimal <- list()      # named by "X1_X2_..."
  minimal_idx <- list()  # store integer indices for fast subset checks

  contains_any_minimal <- function(cols) {
    if (!length(minimal_idx)) return(FALSE)
    for (m in minimal_idx) {
      if (all(m %in% cols)) return(TRUE)
    }
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
        if (isTRUE(stop_at_first)) {
          return(list(minimal_subsets = minimal))
        }
      }
    }
  }

  list(minimal_subsets = minimal)
}

# -------------------------------------------------------------------
# Internal: evaluate both perfect and quasi (row-deletion) conditions
# -------------------------------------------------------------------
#' @keywords internal
.check_sep_both <- function(y, X, epsilon = 1e-5) {
  # helper using the LP feasibility check
  .is_perfect <- function(y, X) {
    res <- feasibility_check_LP(y, X, epsilon = epsilon)
    # lpSolve status 0 = feasible
    (res$ge$status == 0L) || (res$le$status == 0L)
  }

  perfect <- .is_perfect(y, X)

  removed <- integer(0)
  if (!perfect) {
    n <- length(y)
    for (i in seq_len(n)) {
      yi <- y[-i]
      Xi <- X[-i, , drop = FALSE]
      if (.is_perfect(yi, Xi)) removed <- c(removed, i)
    }
  }

  list(
    perfect = perfect,
    quasi   = length(removed) > 0L,
    removed = if (perfect) NULL else removed
  )
}


# -------------------------------------------------------------------
# LP‐based perfect‐separation feasibility check
# -------------------------------------------------------------------
feasibility_check_LP <- function(y, X, epsilon = 1e-5) {
  p <- ncol(X)
  n <- nrow(X)

  # Objective: zeros for alpha + 2*p slack vars
  f.obj <- rep(0, 2 * p + 1)

  # Build constraint matrix A
  A <- matrix(0, nrow = n + 1, ncol = 2 * p + 1)
  A[1:n, 1] <- 1
  for (i in seq_len(n)) {
    A[i, 2:(p+1)]     <- X[i, ]
    A[i, (p+2):(2*p+1)] <- -X[i, ]
  }
  # last constraint: sum(S1_i - S2_i)
  A[n+1, 2:(p+1)]     <- 1
  A[n+1, (p+2):(2*p+1)] <- -1

  # directions and RHS
  dir <- ifelse(y == 1, ">=", "<=")
  dir <- c(dir, ">=")
  rhs <- c(rep(0, n), epsilon)

  # LP1: sum(S1 - S2) ≥ ε
  lp1 <- lpSolve::lp("min", f.obj, A, dir, rhs)
  res_ge <- list(status = lp1$status, sol = lp1$solution)

  # LP2: sum(S1 - S2) ≤ -ε
  dir[n+1] <- "<="
  rhs[n+1] <- -epsilon
  lp2 <- lpSolve::lp("min", f.obj, A, dir, rhs)
  res_le <- list(status = lp2$status, sol = lp2$solution)

  list(ge = res_ge, le = res_le)
}
