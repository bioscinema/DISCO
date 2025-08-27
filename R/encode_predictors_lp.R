# outcome-free encoder for LP (drop-in)
encode_predictors_lp <- function(X) {
  # 1) If numeric matrix, pass through
  if (is.matrix(X) && is.numeric(X)) {
    storage.mode(X) <- "double"
    if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
    return(X)
  }

  # 2) Otherwise coerce to data.frame and recover numerics
  if (is.matrix(X)) {
    X <- as.data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (!is.data.frame(X)) {
    X <- data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    X <- data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  }

  # Turn numeric-looking strings back into numeric; keep non-numeric as char
  X <- type.convert(X, as.is = TRUE)

  # Char -> factor; numeric stays numeric
  chr <- vapply(X, is.character, logical(1))
  if (any(chr)) X[chr] <- lapply(X[chr], factor)

  # 3) Reference coding (k-1 dummies) WITH intercept, then drop intercept
  #    because your LP already has its own intercept column.
  mm <- stats::model.matrix(~ ., data = X)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  }

  storage.mode(mm) <- "double"
  mm
}
