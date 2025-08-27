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
#'   \item \code{separation_type} ("Perfect separation", "Quasi-complete separation", "No separation problem",
#'         "Constant outcome", "Constant predictor")
#'   \item \code{separation_index} (Rand index)
#'   \item \code{severity_score} in \eqn{[0,1]}
#'   \item \code{rand_details} (list from \code{rand_index})
#'   \item \code{single_tie_boundary} (logical), \code{tie_rows_boundary} (integer)
#'   \item \code{boundary_threshold} (numeric)
#'   \item \code{overlap_prop} (numeric in [0,1], optional)
#' }
#'
#' @section Edge Cases:
#' - Constant \code{outcome}: returns early with \code{separation_type = "Constant outcome"}.
#' - Constant \code{predictor}: returns early with \code{separation_type = "Constant predictor"}.
#'
#' @examples
#' toy1 <- data.frame(Y = c(0,0,0,0,1,1,1), X = c(1,2,3,4,4,5,6))
#' toy2 <- data.frame(Y = c(1,1,1,1,0,0,0,0), X = c(1,2,3,4,4,5,6,7))
#' separation(toy1, "X")
#' separation(toy2, "X")
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom clue solve_LSAP
#' @export
uni_separation <- function(data, predictor, outcome = "Y") {
  if (!predictor %in% names(data)) {
    stop(sprintf("Predictor column '%s' not found.", predictor), call. = FALSE)
  }
  if (!outcome %in% names(data)) {
    stop(sprintf("Outcome column '%s' not found.", outcome), call. = FALSE)
  }

  data <- dplyr::filter(
    data,
    !is.na(.data[[predictor]]),
    !is.na(.data[[outcome]])
  )

  y <- data[[outcome]]
  x <- data[[predictor]]

  y <- encode_outcome(y)
  x <- encode_predictors(x)

  # Early exits
  if (length(unique(y)) == 1L) {
    return(list(
      predictor       = predictor,
      outcome         = outcome,
      separation_type = "Constant outcome",
      message         = sprintf("All outcomes %s = %s", outcome, unique(y))
    ))
  }
  if (length(unique(x)) == 1L) {
    return(list(
      predictor       = predictor,
      outcome         = outcome,
      separation_type = "Constant predictor",
      message         = sprintf("All %s = %s", predictor, unique(x))
    ))
  }

  # Overlap proportion (optional diagnostic)
  group0 <- dplyr::filter(data, .data[[outcome]] == 0)
  group1 <- dplyr::filter(data, .data[[outcome]] == 1)
  overlap_prop <- as.numeric(round(
    overlapping::overlap(
      list(x1 = group0[[predictor]], x2 = group1[[predictor]]),
      type   = "2",
      kernel = "gaussian"
    )$OV, 5
  ))

  # Sort by predictor and make initial 2 clusters by outcome balance
  ord_idx <- order(x)
  y_ord <- y[ord_idx]
  x_ord <- x[ord_idx]
  n <- length(y_ord)
  n0 <- sum(y_ord == 0)
  initial_clusters <- c(rep(1L, n0), rep(2L, n - n0))

  cm <- table(
    factor(y_ord, levels = 0:1),
    factor(initial_clusters, levels = 1:2)
  )
  cost_mat <- max(cm) - cm
  opt_map <- clue::solve_LSAP(cost_mat)
  preds <- as.vector(opt_map[initial_clusters] - 1L)

  ri <- rand_index(y_ord, preds)
  sep_idx <- ri$rand_index

  X0 <- sort(unique(group0[[predictor]]))
  X1 <- sort(unique(group1[[predictor]]))
  shared <- intersect(X0, X1)

  single_tie <- length(shared) == 1L && (
    (all(X0 <= shared) && all(X1 >= shared)) ||
      (all(X0 >= shared) && all(X1 <= shared))
  )
  tie_count <- if (isTRUE(single_tie)) sum(x == shared) else 0L

  raw_thresh <- 1 - (2 * tie_count / n)
  boundary_thresh <- max(raw_thresh, 0)

  is_perfect <- isTRUE(!is.na(sep_idx)) && identical(sep_idx, 1) && !single_tie
  sep_type <- if (isTRUE(is_perfect)) {
    "Perfect separation"
  } else if (!is.na(sep_idx) && sep_idx > boundary_thresh) {
    "Quasi-complete separation"
  } else {
    "No separation problem"
  }

  sev_score <- severity_scale(
    sep_idx = sep_idx,
    threshold = boundary_thresh,
    tie_count = tie_count,
    n = n,
    is_perfect = is_perfect
  )

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
