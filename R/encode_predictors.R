# --- Helpers ---------------------------------------------------------------

#' @keywords internal
#' @noRd
# Map any predictor to a numeric axis:
# - numeric: as is
# - logical / 2-level: map to {0,1} in order of increasing P(Y=1|level)
# - >2 levels: score by p_hat = P(Y=1|level) and use that numeric score
encode_predictors <- function(x, y01) {
  if (is.numeric(x)) return(as.numeric(x))

  f <- if (is.factor(x)) x else factor(x)
  levs <- levels(f)
  k <- length(levs)

  # empirical P(Y=1 | level)
  p_hat <- tapply(y01, f, mean)

  # use the probabilities themselves as the numeric axis
  x_num <- as.numeric(p_hat[as.character(f)])

  # ensure deterministic behavior if two levels have identical p_hat:
  # add a tiny, level-specific jitter that preserves equality within a level
  # (doesn't change boundary logic)
  if (anyDuplicated(p_hat)) {
    tiny <- 1e-12 * match(as.character(f), levs)
    x_num <- x_num + tiny
  }

  x_num
}
