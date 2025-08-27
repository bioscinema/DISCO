#' LP-based latent separation detection (with optional minimal-subset search)
#'
#' @param y Binary outcome (0/1 or 2-level factor/character/logical).
#' @param X Matrix or data frame of predictors (columns = variables).
#' @param test_combinations If TRUE, test every subset of size ≥ min_vars (legacy exhaustive mode).
#' @param min_vars Minimum subset size to consider.
#' @param only_perfect Legacy flag for exhaustive mode (kept for compatibility).
#' @param find_minimal If TRUE, return inclusion-minimal separating subsets using pruning.
#' @param mode "either", "perfect", or "quasi" – what to count as a hit.
#' @param max_vars Optional upper bound on subset size (default: all the way to p).
#' @param stop_at_first If TRUE and find_minimal=TRUE, stop after the first minimal hit.
#' @param missing How to treat missing data: "complete" (drop rows w/ NA in the subset) or "impute" (impute predictors; never impute outcome).
#' @param impute_args Optional list (see \code{uni_separation}).
#' @return If not searching: single result list {type, removed?, message, missing_info}.
#'         If searching: list with $minimal_subsets (named by vars joined with "_"), each containing
#'         {type, vars, removed?, missing_info}.
#' @export
latent_separation <- function(
    y, X,
    test_combinations = FALSE,
    min_vars = 2,
    only_perfect = FALSE,
    find_minimal = FALSE,
    mode = c("either","perfect","quasi"),
    max_vars = NULL,
    stop_at_first = FALSE,
    missing = c("complete","impute"),
    impute_args = list()
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

    out <- .check_sep(y1, X1_mm, mode = mode_local)
    out$missing_info <- list(
      method = mh$params_used$method,
      params = mh$params_used,
      rows_used = mh$rows_used,
      n_used = length(mh$rows_used)
    )
    out
  }

  # === No search: single set check (back-compat path) ===========================
  if (!test_combinations && !find_minimal) {
    cols_all <- seq_len(p)
    res <- run_one(cols_all, if (only_perfect) "perfect" else mode)
    if (res$hit) {
      return(list(
        type = res$type, removed = res$removed,
        message = if (identical(res$type, "perfect separation"))
          "Perfect separation found on the given predictors."
        else
          sprintf("Quasi-complete separation: removing any of {%s} yields perfect.",
                  paste(res$removed, collapse = ", ")),
        missing_info = res$missing_info
      ))
    }
    return(list(
      type = "no separation problem",
      message = "No separation problem detected.",
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
