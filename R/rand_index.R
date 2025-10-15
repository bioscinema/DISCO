#' Vectorized Rand Index (symmetric, NA-safe)
#'
#' Computes the Rand index and pair counts (a, b, c, d) using a contingency
#' table, skipping any pairs with missing labels.
#'
#' @param y_true A vector of true cluster labels (integer, logical, or factor).
#' @param y_pred A vector of predicted cluster labels (integer, logical, or factor).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{rand_index} Rand index in \code{[0,1]} or \code{NA} if \eqn{n < 2}.
#'   \item \code{a}, \code{b}, \code{c}, \code{d} Pair counts.
#' }
#' @examples
#' ri <- rand_index(c(0,0,1,1), c(0,1,1,1))
#' ri$rand_index
#' @export
rand_index <- function(y_true, y_pred) {
  keep <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[keep]
  y_pred <- y_pred[keep]

  n <- length(y_true)
  if (n < 2) {
    return(list(rand_index = NA_real_, a = 0, b = 0, c = 0, d = 0))
  }

  tab <- table(
    factor(y_true, levels = sort(unique(y_true))),
    factor(y_pred, levels = sort(unique(y_pred)))
  )

  a <- sum(choose(tab, 2))
  total_pairs <- choose(n, 2)
  same_true <- sum(choose(rowSums(tab), 2))
  same_pred <- sum(choose(colSums(tab), 2))
  c <- same_true - a
  d <- same_pred - a
  b <- total_pairs - (a + c + d)
  rand_val <- (a + b) / total_pairs

  list(rand_index = as.numeric(rand_val), a = a, b = b, c = c, d = d)
}
