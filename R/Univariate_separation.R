# R/uni_separation.R

#' Univariate Separation Detection
#'
#' Detects perfect, quasi-complete, or no separation for a binary outcome
#' against one predictor. Uses Hungarian matching for initial cluster
#' alignment, a vectorized Rand index, a non-negative boundary threshold,
#' and a continuous severity score.
#'
#' @param data A data frame or tibble.
#' @param predictor Name of the predictor column (string).
#' @param outcome Name of the binary outcome column (string, default "Y").
#' @param missing How to handle missing data for the tested variables.
#'   `"complete"` (recommended) drops rows with NA in `predictor` or `outcome`.
#'   `"impute"` imputes the predictor only; outcome NA is always dropped.
#'   Imputation is intended as a sensitivity analysis because it can create
#'   or mask separation artifacts.
#' @param impute_args Optional list of imputation settings when
#'   `missing = "impute"`. Recognized keys:
#'   \itemize{
#'     \item \code{numeric_method}: "median" (default) or "mean"
#'     \item \code{categorical_method}: "mode" (default) or "missing"
#'     \item \code{logical_method}: "mode" (default) or "missing"
#'     \item \code{custom_fn}: function(data.frame) -> imputed data.frame
#'   }
#'
#' @details
#' Procedure:
#' \enumerate{
#'   \item Validate inputs and slice only the needed columns.
#'   \item Handle missing data per user choice (never impute the outcome).
#'   \item Early exit if the outcome or predictor is constant after handling.
#'   \item Build an initial 2-cluster assignment by sorting the predictor,
#'         align with the outcome via Hungarian matching (\code{clue::solve_LSAP}).
#'   \item Compute the Rand index and detect single-tie boundaries.
#'   \item Apply a non-negative boundary threshold and compute a severity score.
#' }
#'
#' @return A list with fields:
#' \itemize{
#'   \item \code{predictor}, \code{outcome}
#'   \item \code{separation_type} ("Perfect separation",
#'         "Quasi-complete separation", "No separation problem",
#'         "Constant outcome", "Constant predictor")
#'   \item \code{separation_index} (Rand index)
#'   \item \code{severity_score} in \code{[0,1]}
#'   \item \code{rand_details} (list from \code{rand_index})
#'   \item \code{single_tie_boundary} (logical), \code{tie_rows_boundary} (integer)
#'   \item \code{boundary_threshold} (numeric)
#'   \item \code{missing_info}: list with \code{method}, \code{params},
#'         \code{rows_used} (indices in the original data), and \code{n_used}
#' }
#'
#' @importFrom clue solve_LSAP
#' @export
uni_separation <- function(
    data,
    predictor,
    outcome = "Y",
    missing = c("complete", "impute"),
    impute_args = list()
) {
  missing <- match.arg(missing)

  if (!predictor %in% names(data)) {
    stop(sprintf("Predictor column '%s' not found.", predictor), call. = FALSE)
  }
  if (!outcome %in% names(data)) {
    stop(sprintf("Outcome column '%s' not found.", outcome), call. = FALSE)
  }

  if (missing == "impute") {
    warning(
      "uni_separation(): missing = \"impute\" imputes the predictor only. ",
      "For separation diagnosis, missing = \"complete\" is recommended; ",
      "treat imputation results as sensitivity analysis.",
      call. = FALSE
    )
  }

  # Slice only what's needed (preserve original row indexing)
  df <- data[, c(predictor, outcome), drop = FALSE]
  y0 <- df[[outcome]]
  X0 <- df[predictor]  # keep as 1-col data.frame for missing handler consistency

  # Handle missing data per user choice (never impute outcome)
  mh <- .handle_missing_for_subset(
    y = y0,
    X = X0,
    method = missing,
    impute_args = impute_args
  )

  # Normalize outcome to 0/1 (integer)
  y <- encode_outcome(mh$y)

  # Prepare x as numeric vector if possible
  x_raw <- mh$X[[1L]]
  if (is.numeric(x_raw)) {
    x <- x_raw
  } else {
    ux <- unique(x_raw)
    k <- length(ux)
    if (is.logical(x_raw) || k == 2L) {
      x <- as.integer(factor(x_raw, levels = sort(ux))) - 1L
    } else {
      # Multi-level categorical predictor: map via empirical P(Y=1|level)
      x <- encode_predictors(x_raw, y)
    }
  }

  # Early exits after missing handling
  if (length(unique(y)) == 1L) {
    return(list(
      predictor = predictor,
      outcome = outcome,
      separation_type = "Constant outcome",
      separation_index = NA_real_,
      severity_score = NA_real_,
      rand_details = NULL,
      single_tie_boundary = NA,
      tie_rows_boundary = NA_integer_,
      boundary_threshold = NA_real_,
      message = sprintf(
        "All outcomes %s = %s after missing-data handling.",
        outcome, paste(unique(y), collapse = ", ")
      ),
      missing_info = list(
        method = mh$params_used$method,
        params = mh$params_used,
        rows_used = mh$rows_used,
        n_used = length(mh$rows_used)
      )
    ))
  }
  if (length(unique(x)) == 1L) {
    return(list(
      predictor = predictor,
      outcome = outcome,
      separation_type = "Constant predictor",
      separation_index = NA_real_,
      severity_score = NA_real_,
      rand_details = NULL,
      single_tie_boundary = NA,
      tie_rows_boundary = NA_integer_,
      boundary_threshold = NA_real_,
      message = sprintf(
        "All %s = %s after missing-data handling.",
        predictor, paste(unique(x), collapse = ", ")
      ),
      missing_info = list(
        method = mh$params_used$method,
        params = mh$params_used,
        rows_used = mh$rows_used,
        n_used = length(mh$rows_used)
      )
    ))
  }

  # --- Use encoded y to split groups ---
  idx0 <- which(y == 0L)
  idx1 <- which(y == 1L)

  # Sort by predictor and build initial clusters
  ord_idx <- order(x)
  y_ord <- y[ord_idx]
  n <- length(y_ord)
  n0 <- sum(y_ord == 0L)
  initial_clusters <- c(rep(1L, n0), rep(2L, n - n0))

  cm <- table(factor(y_ord, levels = 0:1), factor(initial_clusters, levels = 1:2))
  cost_mat <- max(cm) - cm
  opt_map <- clue::solve_LSAP(cost_mat)
  preds <- as.vector(opt_map[initial_clusters] - 1L)

  ri <- rand_index(y_ord, preds)
  sep_idx <- ri$rand_index

  X0u <- sort(unique(x[idx0]))
  X1u <- sort(unique(x[idx1]))
  shared <- intersect(X0u, X1u)

  eps_tie <- 1e-8
  tol_perfect <- 1e-12

  single_tie <- length(shared) == 1L && (
    (all(X0u <= shared + eps_tie) && all(X1u >= shared - eps_tie)) ||
      (all(X0u >= shared - eps_tie) && all(X1u <= shared + eps_tie))
  )
  tie_count <- if (isTRUE(single_tie)) sum(abs(x - shared) < eps_tie) else 0L

  raw_thresh <- 1 - (2 * tie_count / n)
  boundary_thresh <- max(raw_thresh, 0)

  is_perfect <- !is.na(sep_idx) && (sep_idx >= 1 - tol_perfect) && !single_tie
  sep_type <- if (isTRUE(is_perfect)) {
    "Perfect separation"
  } else if (!is.na(sep_idx) && sep_idx > boundary_thresh) {
    "Quasi-complete separation"
  } else {
    "No separation problem"
  }

  sev_score <- severity_scale(sep_idx, boundary_thresh, tie_count, n, is_perfect)

  out <- list(
    predictor = predictor,
    outcome = outcome,
    separation_type = sep_type,
    separation_index = sep_idx,
    severity_score = sev_score,
    rand_details = ri,
    single_tie_boundary = single_tie,
    tie_rows_boundary = tie_count,
    boundary_threshold = boundary_thresh,
    missing_info = list(
      method = mh$params_used$method,
      params = mh$params_used,
      rows_used = mh$rows_used,
      n_used = length(mh$rows_used)
    )
  )

  class(out) <- "uni_separation"
  out
}

#' @export
print.uni_separation <- function(x, digits = 3, ...) {
  fmt <- function(z) {
    if (is.na(z)) return("NA")
    sprintf("%.*f", digits, z)
  }

  cat("Univariate separation for", x$predictor, "vs", x$outcome, "\n")

  if (!is.null(x$message)) {
    cat("  ", x$message, "\n", sep = "")
    return(invisible(x))
  }

  cat("  separation type:       ", x$separation_type, "\n", sep = "")
  cat("  separation index (SI): ", fmt(x$separation_index), "\n", sep = "")
  cat("  severity score:        ", fmt(x$severity_score), "\n", sep = "")
  cat("  boundary threshold τ+: ", fmt(x$boundary_threshold), "\n", sep = "")
  cat("  boundary tie rows nb:  ", x$tie_rows_boundary, "\n", sep = "")
  cat("  single-tie boundary:   ", x$single_tie_boundary, "\n", sep = "")
  invisible(x)
}
