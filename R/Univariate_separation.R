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
#'
#' @details
#' The procedure:
#' \enumerate{
#'   \item Validates inputs; drops rows with missing \code{predictor} or \code{outcome}.
#'   \item Early exits if the outcome or predictor is constant.
#'   \item Optionally computes group overlap via kernel density (\pkg{overlapping}).
#'   \item Builds an initial 2-cluster assignment by sorting \code{predictor},
#'         aligns with \code{outcome} using Hungarian matching
#'         (\code{clue::solve_LSAP}).
#'   \item Computes Rand index and detects boundary single-ties.
#'   \item Applies a non-negative threshold and computes a severity score.
#' }
#'
#' @return A list with fields:
#' \itemize{
#'   \item \code{predictor}, \code{outcome}
#'   \item \code{separation_type} ("Perfect separation", "Quasi-complete separation",
#'         "No separation problem", "Constant outcome", "Constant predictor")
#'   \item \code{separation_index} (Rand index)
#'   \item \code{severity_score} in \eqn{[0,1]}
#'   \item \code{rand_details} (list from \code{rand_index})
#'   \item \code{single_tie_boundary} (logical), \code{tie_rows_boundary} (integer)
#'   \item \code{boundary_threshold} (numeric)
#'   \item \code{overlap_prop} (numeric in [0,1], optional)
#' }
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom clue solve_LSAP
#' @export
uni_separation <- function(data, predictor, outcome = "Y") {
  if (!predictor %in% names(data)) stop(sprintf("Predictor column '%s' not found.", predictor), call. = FALSE)
  if (!outcome %in% names(data))   stop(sprintf("Outcome column '%s' not found.", outcome),   call. = FALSE)

  # drop missing
  data <- dplyr::filter(data, !is.na(.data[[predictor]]), !is.na(.data[[outcome]]))

  # normalize outcome to 0/1 (integer)
  y <- encode_outcome(data[[outcome]])

  # Prepare x as numeric vector if possible
  x_raw <- data[[predictor]]
  if (is.numeric(x_raw)) {
    x <- x_raw
  } else {
    ux <- unique(x_raw); k <- length(ux)
    if (is.logical(x_raw) || k == 2L) {
      x <- as.integer(factor(x_raw, levels = sort(ux))) - 1L
    } else {
      # Multi-level categorical predictor: run LP on its dummies and summarize
      Xdum <- encode_predictors(data[predictor, drop = FALSE])
      lp_res <- latent_separation(y = y, X = Xdum, test_combinations = FALSE)
      return(list(
        predictor           = predictor,
        outcome             = outcome,
        separation_type     = switch(lp_res$type,
                                     "perfect separation" = "Perfect separation",
                                     "quasi-complete separation" = "Quasi-complete separation",
                                     "no separation problem" = "No separation problem",
                                     lp_res$type),
        separation_index    = NA_real_,
        severity_score      = NA_real_,
        rand_details        = NULL,
        single_tie_boundary = FALSE,
        tie_rows_boundary   = 0L,
        boundary_threshold  = NA_real_,
        overlap_prop        = NA_real_,
        note                = "Categorical predictor auto-encoded to dummies for LP-based check."
      ))
    }
  }


  # Early exits
  if (length(unique(y)) == 1L) {
    return(list(predictor=predictor, outcome=outcome, separation_type="Constant outcome",
                message=sprintf("All outcomes %s = %s", outcome, unique(y))))
  }
  if (length(unique(x)) == 1L) {
    return(list(predictor=predictor, outcome=outcome, separation_type="Constant predictor",
                message=sprintf("All %s = %s", predictor, unique(x))))
  }

  # --- Use encoded y to split groups (fixes your bug) ---
  idx0 <- which(y == 0L)
  idx1 <- which(y == 1L)

  # Optional: compute overlap safely
  overlap_prop <- NA_real_
  if (length(idx0) >= 2L && length(idx1) >= 2L) {
    overlap_prop <- as.numeric(round(
      overlapping::overlap(list(x1 = x[idx0], x2 = x[idx1]), type = "2", kernel = "gaussian")$OV, 5
    ))
  }

  # Sort by predictor and build initial clusters
  ord_idx <- order(x)
  y_ord   <- y[ord_idx]
  x_ord   <- x[ord_idx]
  n       <- length(y_ord)
  n0      <- sum(y_ord == 0L)
  initial_clusters <- c(rep(1L, n0), rep(2L, n - n0))

  cm <- table(factor(y_ord, levels = 0:1), factor(initial_clusters, levels = 1:2))
  cost_mat <- max(cm) - cm
  opt_map  <- clue::solve_LSAP(cost_mat)
  preds    <- as.vector(opt_map[initial_clusters] - 1L)

  ri      <- rand_index(y_ord, preds)
  sep_idx <- ri$rand_index

  X0 <- sort(unique(x[idx0]))
  X1 <- sort(unique(x[idx1]))
  shared <- intersect(X0, X1)

  single_tie <- length(shared) == 1L && (
    (all(X0 <= shared) && all(X1 >= shared)) ||
      (all(X0 >= shared) && all(X1 <= shared))
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
    overlap_prop        = overlap_prop
  )
}
