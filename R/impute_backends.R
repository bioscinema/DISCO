.impute_with_backend <- function(X, backend, backend_args = list()) {
  backend <- tolower(backend)
  X <- as.data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)

  # convert characters to factors (backends expect factors, not chars)
  chr <- vapply(X, is.character, logical(1))
  if (any(chr)) X[chr] <- lapply(X[chr], factor, exclude = NULL)

  p <- ncol(X)

  # --- missForest ------------------------------------------------------------
  if (backend == "missforest") {
    if (p < 2) {
      warning("missForest requires ≥2 columns; falling back to rule-based imputer.")
      return(.impute_rule_based(X))
    }
    if (!requireNamespace("missForest", quietly = TRUE))
      stop("Please install 'missForest' to use backend = 'missForest'.", call. = FALSE)

    res <- try(
      do.call(missForest::missForest, c(list(xmis = X, verbose = FALSE), backend_args)),
      silent = TRUE
    )
    if (inherits(res, "try-error")) {
      warning("missForest failed on this data; falling back to rule-based imputer.")
      return(.impute_rule_based(X))
    }
    return(res$ximp)
  }

  # --- VIM::kNN --------------------------------------------------------------
  if (backend %in% c("vim_knn","vim-knn","knn")) {
    if (p < 2) {
      warning("VIM::kNN requires ≥2 columns; falling back to rule-based imputer.")
      return(.impute_rule_based(X))
    }
    if (!requireNamespace("VIM", quietly = TRUE))
      stop("Please install 'VIM' to use backend = 'vim_knn'.", call. = FALSE)
    res <- try(
      do.call(VIM::kNN, c(list(data = X, variable = colnames(X), imp_var = FALSE), backend_args)),
      silent = TRUE
    )
    if (inherits(res, "try-error")) {
      warning("VIM::kNN failed on this data; falling back to rule-based imputer.")
      return(.impute_rule_based(X))
    }
    return(res)
  }

  # --- mice ------------------------------------------------------------------
  if (backend == "mice") {
    if (p < 2) {
      warning("mice requires ≥2 columns; falling back to rule-based imputer.")
      return(.impute_rule_based(X))
    }
    if (!requireNamespace("mice", quietly = TRUE))
      stop("Please install 'mice' to use backend = 'mice'.", call. = FALSE)
    res <- try(
      do.call(mice::mice, c(list(data = X, m = 1, maxit = 5, printFlag = FALSE), backend_args)),
      silent = TRUE
    )
    if (inherits(res, "try-error")) {
      warning("mice failed on this data; falling back to rule-based imputer.")
      return(.impute_rule_based(X))
    }
    return(mice::complete(res, 1))
  }

  # --- recipes ---------------------------------------------------------------
  if (backend == "recipes") {
    if (!requireNamespace("recipes", quietly = TRUE))
      stop("Please install 'recipes' to use backend = 'recipes'.", call. = FALSE)

    method <- tolower(backend_args$method %||% "median_mode")
    # If KNN is requested but only 1 column, downgrade to median/mode
    if (p < 2 && method == "knn") {
      warning("recipes::step_impute_knn requires ≥2 columns; using median/mode instead.")
      method <- "median_mode"
    }

    rec <- recipes::recipe(~ ., data = X)
    if (method == "median_mode") {
      rec <- rec |>
        recipes::step_impute_median(recipes::all_numeric_predictors()) |>
        recipes::step_impute_mode(recipes::all_nominal_predictors(), recipes::all_logical_predictors())
    } else if (method == "knn") {
      rec <- rec |>
        recipes::step_impute_knn(recipes::all_predictors(),
                                 neighbors = backend_args$neighbors %||% 5)
    } else if (method %in% c("bag","bagimpute","bagged")) {
      rec <- rec |>
        recipes::step_impute_bag(recipes::all_predictors())
    } else {
      stop("recipes backend: unknown method '", method, "'.", call. = FALSE)
    }

    res <- try(recipes::prep(rec, training = X, verbose = FALSE), silent = TRUE)
    if (inherits(res, "try-error")) {
      warning("recipes backend failed on this data; falling back to rule-based imputer.")
      return(.impute_rule_based(X))
    }
    return(recipes::bake(res, new_data = NULL))
  }

  stop("Unknown imputation backend: '", backend,
       "'. Supported: 'mice', 'vim_knn', 'missForest', 'recipes'.", call. = FALSE)
}
