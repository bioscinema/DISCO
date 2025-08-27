# --- Missing-data helpers ------------------------------------------------------

#' @keywords internal
#' @noRd
.mode_value <- function(x) {
  ux <- unique(x[!is.na(x)])
  if (length(ux) == 0L) return(NA)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Simple, dependency-free imputation for a data.frame of predictors
#' @param X data.frame or matrix (will be returned with original classes where possible)
#' @param numeric_method "median" or "mean"
#' @param categorical_method "mode" or "missing" (treat NA as its own level)
#' @param logical_method "mode" or "missing"
#' @param custom_fn optional user-supplied function(X) -> X_imputed
#' @keywords internal
#' @noRd
.impute_predictors_df <- function(
    X,
    numeric_method = c("median","mean"),
    categorical_method = c("mode","missing"),
    logical_method = c("mode","missing"),
    custom_fn = NULL
) {
  if (!is.null(custom_fn)) return(custom_fn(X))

  numeric_method      <- match.arg(numeric_method)
  categorical_method  <- match.arg(categorical_method)
  logical_method      <- match.arg(logical_method)

  X <- as.data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)

  for (j in names(X)) {
    col <- X[[j]]

    # logical
    if (is.logical(col)) {
      if (anyNA(col)) {
        if (identical(logical_method, "missing")) {
          # represent missing as FALSE/TRUE/NA; for modeling we prefer no NA, so use mode
          m <- .mode_value(col)
          col[is.na(col)] <- if (is.na(m)) FALSE else m
        } else {
          m <- .mode_value(col)
          col[is.na(col)] <- if (is.na(m)) FALSE else m
        }
      }
      X[[j]] <- col
      next
    }

    # numeric / integer
    if (is.numeric(col)) {
      if (anyNA(col)) {
        fill <- if (identical(numeric_method, "median")) stats::median(col, na.rm = TRUE) else mean(col, na.rm = TRUE)
        if (is.na(fill)) fill <- 0
        col[is.na(col)] <- fill
      }
      X[[j]] <- as.numeric(col)
      next
    }

    # everything else -> treat as categorical
    if (!is.factor(col)) {
      col <- factor(col, exclude = NULL)
    }

    if (anyNA(col)) {
      if (identical(categorical_method, "missing")) {
        # explicitly add a "Missing" level
        if (!("Missing" %in% levels(col))) levels(col) <- c(levels(col), "Missing")
        col[is.na(col)] <- "Missing"
      } else {
        m <- .mode_value(col)
        if (is.na(m)) {
          # if everything NA, create a fill level
          levels(col) <- c(levels(col), "Filled")
          col[is.na(col)] <- "Filled"
        } else {
          col[is.na(col)] <- m
        }
      }
    }
    X[[j]] <- droplevels(col)
  }

  X
}

#' Handle missing data for a given outcome + predictors
#' @param y outcome vector (may contain NA; NA rows are always dropped)
#' @param X predictor frame for *current* test (columns already chosen)
#' @param method "complete" or "impute"
#' @param impute_args list: numeric_method, categorical_method, logical_method, custom_fn
#' @return list(y=, X=, rows_used=, params_used=)
#' @keywords internal
#' @noRd
.handle_missing_for_subset <- function(
    y, X,
    method = c("complete","impute"),
    impute_args = list()
) {
  method <- match.arg(method)

  # always drop rows with missing outcome
  keep_y <- !is.na(y)
  y0 <- y[keep_y]
  X0 <- X[keep_y, , drop = FALSE]
  rows_after_y <- which(keep_y)

  if (identical(method, "complete")) {
    keep <- stats::complete.cases(X0)
    list(
      y = y0[keep],
      X = X0[keep, , drop = FALSE],
      rows_used = rows_after_y[keep],
      params_used = list(method = "complete")
    )
  } else {
    # impute predictors only (never impute outcome)
    numeric_method     <- impute_args$numeric_method     %||% "median"
    categorical_method <- impute_args$categorical_method %||% "mode"
    logical_method     <- impute_args$logical_method     %||% "mode"
    custom_fn          <- impute_args$custom_fn          %||% NULL

    X_imp <- .impute_predictors_df(
      X0,
      numeric_method = numeric_method,
      categorical_method = categorical_method,
      logical_method = logical_method,
      custom_fn = custom_fn
    )

    list(
      y = y0,
      X = X_imp,
      rows_used = rows_after_y,
      params_used = list(
        method = "impute",
        numeric_method = if (is.null(custom_fn)) numeric_method else NA,
        categorical_method = if (is.null(custom_fn)) categorical_method else NA,
        logical_method = if (is.null(custom_fn)) logical_method else NA,
        custom_fn = if (!is.null(custom_fn)) "user_supplied_function" else NULL
      )
    )
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b


# ---- Internal: simple rule-based imputer (fallback) ---------------------------
.impute_rule_based <- function(X, numeric_method = "median",
                               categorical_method = "mode",
                               logical_method = "mode") {
  X <- as.data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)

  # Ensure characters are factors for "mode"/"missing" handling
  chr <- vapply(X, is.character, logical(1))
  if (any(chr)) X[chr] <- lapply(X[chr], factor, exclude = NULL)

  for (j in names(X)) {
    xj <- X[[j]]
    if (is.numeric(xj)) {
      if (numeric_method == "mean") fill <- mean(xj, na.rm = TRUE) else fill <- stats::median(xj, na.rm = TRUE)
      if (is.na(fill)) fill <- 0
      xj[is.na(xj)] <- fill
    } else if (is.logical(xj)) {
      if (logical_method == "missing") {
        xj <- factor(xj, exclude = NULL)
        if (!("Missing" %in% levels(xj))) levels(xj) <- c(levels(xj), "Missing")
        xj[is.na(xj)] <- "Missing"
        xj <- droplevels(xj)
      } else {
        # mode (most common non-NA), default FALSE if all NA
        tab <- sort(table(xj, useNA = "no"), decreasing = TRUE)
        fill <- if (length(tab)) as.logical(names(tab)[1]) else FALSE
        xj[is.na(xj)] <- fill
      }
    } else { # treat as categorical
      f <- if (is.factor(xj)) xj else factor(xj, exclude = NULL)
      if (categorical_method == "missing") {
        if (!("Missing" %in% levels(f))) levels(f) <- c(levels(f), "Missing")
        f[is.na(f)] <- "Missing"
      } else {
        tab <- sort(table(f, useNA = "no"), decreasing = TRUE)
        if (length(tab)) f[is.na(f)] <- names(tab)[1]
      }
      xj <- droplevels(f)
    }
    X[[j]] <- xj
  }
  X
}

.impute_predictors_df <- function(X, impute_args = list()) {
  # 1) custom function wins
  if (!is.null(impute_args$custom_fn)) {
    fn <- impute_args$custom_fn
    stopifnot(is.function(fn))
    Xi <- fn(X)
    return(as.data.frame(Xi, check.names = FALSE))
  }

  # 2) backend if provided
  if (!is.null(impute_args$backend)) {
    Xi <- .impute_with_backend(X, backend = impute_args$backend,
                               backend_args = impute_args$backend_args %||% list())
    return(as.data.frame(Xi, check.names = FALSE))
  }

  # 3) default rule-based
  .impute_rule_based(
    X,
    numeric_method     = tolower(impute_args$numeric_method     %||% "median"),
    categorical_method = tolower(impute_args$categorical_method %||% "mode"),
    logical_method     = tolower(impute_args$logical_method     %||% "mode")
  )
}

# ---- Internal: per-subset missing handler (y never imputed) -------------------
.handle_missing_for_subset <- function(y, X, method = c("complete","impute"),
                                       impute_args = list()) {
  method <- match.arg(method)
  X <- as.data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)

  # Always drop rows with NA in outcome
  keep_y <- !is.na(y)
  y2 <- y[keep_y]
  X2 <- X[keep_y, , drop = FALSE]
  idx_y <- which(keep_y)

  if (method == "complete") {
    keep <- stats::complete.cases(X2)
    list(
      y = y2[keep],
      X = X2[keep, , drop = FALSE],
      rows_used = idx_y[keep],
      params_used = list(method = "complete")
    )
  } else {
    Xi <- .impute_predictors_df(X2, impute_args = impute_args)
    params <- list(
      method = "impute",
      numeric_method     = impute_args$numeric_method     %||% NULL,
      categorical_method = impute_args$categorical_method %||% NULL,
      logical_method     = impute_args$logical_method     %||% NULL,
      backend            = impute_args$backend            %||% NULL,
      backend_args       = impute_args$backend_args       %||% NULL,
      custom_fn          = if (!is.null(impute_args$custom_fn)) "custom_fn" else NULL
    )
    list(
      y = y2,
      X = Xi,
      rows_used  = idx_y,
      params_used = params
    )
  }
}
