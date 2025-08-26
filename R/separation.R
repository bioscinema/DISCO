#' Detect separation for a single predictor
#'
#' Evaluate whether a numeric predictor exhibits complete, quasi-complete, or no separation with respect to a binary outcome.
#'
#' @param data A data frame containing columns \code{Y} (0/1) and the predictor.
#' @param order_by Character. Name of the predictor column to analyze.
#' @return A list with:
#' \itemize{
#'   \item separation_type: "Complete (Perfect) Separation", "Quasi Complete Separation", or "No Separation Problem".
#'   \item separation_index: Rand index value.
#'  \item rand_details: output of \code{rand_index}.
#'   \item proportion_overlap: numeric overlap proportion from the \code{overlapping} package.
#'   \item shared_values: predictor values present in both outcome groups.
#'   \item single_tie_boundary: logical indicating a single boundary tie.
#'   \item tie_rows_boundary: count of rows at the boundary value.
#'   \item boundary_threshold: threshold used to categorize quasi separation.
#' }
#' @import dplyr overlapping
#' @export
separation <- function(data, order_by) {
  if (!order_by %in% names(data)) stop("Column 'order_by' not found in data.")
  vals <- unique(data[[order_by]])
  if (length(vals) == 1L) {
    return(list(
      separation_type = "Same value for all observations",
      message = paste("All observations have", order_by, "=", vals)
    ))
  }
  # Prepare data and prediction vector
  if (length(vals) == 2L) {
    ordered_df <- data[order(data$Y, data[[order_by]]), ]
    true_vec <- ordered_df$Y
    pred_vec <- ordered_df[[order_by]]
  } else {
    ordered_df <- data[order(data$Y), ]
    ordered_df <- ordered_df[order(ordered_df[[order_by]]), ]
    true_vec <- ordered_df$Y
    n0 <- sum(true_vec == 0)
    n1 <- length(true_vec) - n0
    pred_vec <- c(rep(0, n0), rep(1, n1))
  }
  # Compute Rand index and details
  rand_details <- rand_index(true_vec, pred_vec)
  ri <- rand_details$rand_index
  # Boundary logic
  group0_vals <- sort(unique(data[[order_by]][data$Y == 0]))
  group1_vals <- sort(unique(data[[order_by]][data$Y == 1]))
  shared <- intersect(group0_vals, group1_vals)
  single_tie <- length(shared) == 1L && (
    all(group0_vals <= shared) && all(group1_vals >= shared) ||
      all(group0_vals >= shared) && all(group1_vals <= shared)
  )
  tie_count <- if (single_tie) sum(data[[order_by]] == shared) else 0L
  threshold <- 1 - (tie_count / nrow(data))
  # Determine separation type
  separation_type <- if (ri == 1 && !single_tie) {
    "Complete (Perfect) Separation"
  } else if (ri > threshold) {
    "Quasi Complete Separation"
  } else {
    "No Separation Problem"
  }
  list(
    separation_type     = separation_type,
    separation_index    = ri,
#    rand_details        = rand_details,
#    shared_values       = shared,
    single_tie_boundary = single_tie,
    tie_rows_boundary   = tie_count,
    boundary_threshold  = threshold
  )
}
