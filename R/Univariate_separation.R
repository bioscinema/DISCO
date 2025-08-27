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
#'   One of `"complete"` (drop rows with NA in `predictor` or `outcome`) or
#'   `"impute"` (impute the predictor only; outcome NA is always dropped).
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
#'   \item \code{severity_score} in [0,1]
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

  # Slice only what's needed (preserve original row indexing)
  df   <- data[, c(predictor, outcome), drop = FALSE]
  y0   <- df[[outcome]]
  X0   <- df[predictor]

  # Handle missing data per user choice (never impute outcome)
  mh <- .handle_missing_for_subset(
    y = y0, X = X0,
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
    ux <- unique(x_raw); k <- length(ux)
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
      predictor = predictor, outcome = outcome,
      separation_type = "Constant outcome",
      message = sprintf("All outcomes %s = %s after missing-data handling.", outcome, unique(y)),
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
      predictor = predictor, outcome = outcome,
      separation_type = "Constant predictor",
      message = sprintf("All %s = %s after missing-data handling.", predictor, unique(x)),
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
  y_ord   <- y[ord_idx]
  n       <- length(y_ord)
  n0      <- sum(y_ord == 0L)
  initial_clusters <- c(rep(1L, n0), rep(2L, n - n0))

  cm <- table(factor(y_ord, levels = 0:1), factor(initial_clusters, levels = 1:2))
  cost_mat <- max(cm) - cm
  opt_map  <- clue::solve_LSAP(cost_mat)
  preds    <- as.vector(opt_map[initial_clusters] - 1L)

  ri      <- rand_index(y_ord, preds)
  sep_idx <- ri$rand_index

  X0u <- sort(unique(x[idx0]))
  X1u <- sort(unique(x[idx1]))
  shared <- intersect(X0u, X1u)

  single_tie <- length(shared) == 1L && (
    (all(X0u <= shared) && all(X1u >= shared)) ||
      (all(X0u >= shared) && all(X1u <= shared))
  )
  tie_count <- if (isTRUE(single_tie)) sum(x == shared) else 0L

  raw_thresh      <- 1 - (2 * tie_count / n)
  boundary_thresh <- max(raw_thresh, 0)

  is_perfect <- isTRUE(!is.na(sep_idx)) && identical(sep_idx, 1) && !single_tie
  sep_type <- if (isTRUE(is_perfect)) {
    "Perfect separation"
  } else if (!is.na(sep_idx) && sep_idx > boundary_thresh) {
    "Quasi-complete separation"
  } else {
    "No separation problem"
  }

  sev_score <- severity_scale(sep_idx, boundary_thresh, tie_count, n, is_perfect)

  list(
    predictor           = predictor,
    outcome             = outcome,
    separation_type     = sep_type,
    separation_index    = sep_idx,
    severity_score      = sev_score,
    rand_details        = ri,
    single_tie_boundary = single_tie,
    tie_rows_boundary   = tie_count,
    boundary_threshold  = boundary_thresh,
    missing_info        = list(
      method = mh$params_used$method,
      params = mh$params_used,
      rows_used = mh$rows_used,
      n_used = length(mh$rows_used)
    )
  )
}
