#' Severity Score Rescaling
#'
#' Rescales a separation index above a threshold into \eqn{[0,1]}, with
#' optional attenuation when a single boundary tie inflates the index.
#'
#' @param sep_idx Numeric separation index (e.g., Rand index).
#' @param threshold Numeric threshold in \eqn{[0,1]} (already clipped at 0).
#' @param tie_count Integer number of observations at the single boundary tie.
#' @param n Integer sample size used to compute \code{sep_idx}.
#' @param is_perfect Logical; set to \code{TRUE} when perfect separation holds.
#'
#' @return A numeric severity score in \eqn{[0,1]}.
#' @keywords internal
#' @noRd
severity_scale <- function(sep_idx, threshold, tie_count, n, is_perfect) {
  if (isTRUE(is_perfect)) return(1)
  if (is.na(sep_idx) || sep_idx <= threshold) return(0)

  raw <- (sep_idx - threshold) / (1 - threshold)
  raw * (1 - tie_count / n)
}

