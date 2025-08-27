# --- Helpers ---------------------------------------------------------------

#' Normalize binary outcome y to 0/1 integer
#' @keywords internal
#' @noRd
encode_outcome <- function(y) {
  if (is.logical(y)) return(as.integer(y))
  if (is.factor(y) || is.character(y)) {
    f <- factor(y)
    if (nlevels(f) != 2L) {
      stop("Outcome must be binary (2 levels). Got levels: ",
           paste(levels(f), collapse = ", "))
    }
    return(as.integer(f) - 1L)  # first level -> 0, second -> 1
  }
  if (is.numeric(y)) {
    u <- sort(unique(y))
    if (!all(u %in% c(0, 1))) {
      stop("Numeric outcome must be coded 0/1. Got: ", paste(u, collapse = ", "))
    }
    return(as.integer(y))
  }
  stop("Unsupported outcome type: ", class(y)[1])
}
