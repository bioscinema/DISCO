# --- Helpers ---------------------------------------------------------------

#' Encode predictors to a numeric matrix of dummies (no intercept)
#' Keeps full indicator set (no baseline dropped), useful for separation checks.
#' @keywords internal
#' @noRd
encode_predictors <- function(X) {
  if (is.matrix(X)) X <- as.data.frame(X, check.names = FALSE)
  stopifnot(is.data.frame(X))
  mm <- stats::model.matrix(
    ~ . - 1,
    data = X,
    contrasts.arg = lapply(X, function(z) FALSE)
  )
  storage.mode(mm) <- "double"
  mm
}
