#' Compute the Rand index for two clusterings
#'
#' @param y_true A vector of true class labels (0/1 or factors).
#' @param y_pred A vector of predicted class labels or assignments.
#' @return A list with elements \code{rand_index}, \code{a}, \code{b}, \code{c}, \code{d}.
#' @keywords internal
#' @noRd
rand_index <- function(y_true, y_pred) {
  n <- length(y_true)
  a <- b <- c <- d <- 0
  for (i in seq_len(n-1)) {
    for (j in (i+1):n) {
      same_true <- y_true[i] == y_true[j]
      same_pred <- y_pred[i] == y_pred[j]
      if (same_true && same_pred) {
        a <- a + 1
      } else if (!same_true && !same_pred) {
        b <- b + 1
      } else if (same_true && !same_pred) {
        c <- c + 1
      } else if (!same_true && same_pred) {
        d <- d + 1
      }
    }
  }
  rand_val <- (a + b) / (a + b + c + d)
  list(rand_index = rand_val, a = a, b = b, c = c, d = d)
}
