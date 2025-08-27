# --- Helper: encode non-numeric predictors to numeric matrix of dummies ----
#' @keywords internal
#' @noRd
encode_predictors <- function(X) {
  # Accept matrix or data.frame and return a numeric matrix with no intercept
  if (is.matrix(X)) {
    X <- as.data.frame(X, check.names = FALSE, stringsAsFactors = FALSE)
  }
  stopifnot(is.data.frame(X))

  # model.matrix handles factors/characters/logicals nicely
  # Use contrasts = FALSE to avoid dropping a baseline — for separation
  # diagnostics it's often preferable to keep the full indicator set.
  mm <- stats::model.matrix(~ . - 1, data = X,
                            contrasts.arg = lapply(X, function(z) FALSE))
  # Ensure it's strictly numeric matrix
  storage.mode(mm) <- "double"
  mm
}
