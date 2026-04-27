# ============================================================
# LP-based latent separation detection (full updated version)
#   Updates:
#   1) missing_scope = "global" option to fix samples across subsets
#   2) minimal_strategy = auto/forward/backward
#   3) backward search supports beam search and layered exhaustive search
#   4) beam_width contro_
#   5) verbose option for progress and stop-reason messages
# ============================================================

#' LP-based latent separation detection with optional subset search
#'
#' Detects latent complete or quasi-complete separation for a binary outcome
#' using linear programming in the multivariate predictor space.
#'
#' The method uses a two-stage LP pipeline.
#'
#' Stage A solves a max-margin linear program with L1 equality normalization
#' to identify complete separation.
#'
#' Stage B solves a severity linear program to compute the multivariate
#' severity lower bound, `K_relax`, for quasi-complete separation.
#'
#' Quasi-complete separation is treated as non-substantive and labeled as
#' no separation when `K_relax / n` is large. By default, this threshold is
#' controlled by `quasi_to_none_if = 0.5`. This rule is intended to reduce
#' false positives from weak or non-substantive quasi-complete separation.
#'
#' A legacy row-deletion diagnostic for quasi-complete separation is retained.
#' The returned `removed` indices indicate observations whose deletion would
#' yield complete separation, but these indices do not determine the quasi
#' decision.
#'
#' Missing data can be handled either globally or separately within each
#' subset evaluation. With `missing_scope = "global"`, missing-data handling is
#' performed once on the full predictor set, and the resulting fixed sample is
#' reused for every subset. This is the recommended setting when comparing
#' different predictor subsets. With `missing_scope = "subset"`, missing-data
#' handling is performed separately for each subset, which preserves the legacy
#' behavior but may evaluate different subsets on different samples.
#'
#' Minimal subset search is controlled by `minimal_strategy`. The `"forward"`
#' strategy enumerates subsets in increasing size and returns minimal separating
#' subsets, but it can be computationally expensive. The `"backward"` strategy
#' starts from the full predictor set and removes predictors layer by layer.
#' By default, backward search uses a beam-style search controlled by
#' `beam_width`, retaining the best separating subsets at each layer. Ties at
#' the beam cutoff are retained. If `backward_exhaustive = TRUE`, backward
#' search instead performs layered exhaustive search, keeping all separating
#' subsets at each layer, subject to the evaluation limit. The `"auto"` strategy
#' uses forward search when `p <= small_p_threshold` and backward search
#' otherwise.
#'
#' Runtime can be controlled with:
#'
#' `options(latent_separation.eval_limit = N)`
#'
#' `options(latent_separation.progress_every = M)`
#'
#' Progress reporting can be enabled with `verbose = TRUE`. This prints progress
#' updates and key stopping reasons during minimal-subset search. Advanced users
#' may also use `options(latent_separation.show_progress = TRUE/FALSE)`, but
#' `verbose` is the recommended user-facing switch.
#'
#' @param y Binary outcome, supplied as 0/1, a two-level factor, a two-level
#'   character vector, or a logical vector.
#' @param X Matrix or data frame of predictors. Columns correspond to variables.
#' @param test_combinations Logical. If `TRUE`, test all predictor subsets of
#'   size at least `min_vars`.
#' @param min_vars Minimum subset size to consider. Default is 2.
#' @param epsilon Numeric scalar greater than 0. Legacy margin parameter; used
#'   as the default value for `eps_boundary` when `eps_boundary = NULL`.
#' @param only_perfect Logical. Legacy flag for exhaustive mode. If `TRUE`, only
#'   perfect separation is returned.
#' @param find_minimal Logical. If `TRUE`, search for separating subsets using
#'   the strategy specified by `minimal_strategy`.
#' @param mode Character string specifying the target separation type. One of
#'   `"either"`, `"perfect"`, or `"quasi"`.
#' @param max_vars Optional upper bound on subset size. Defaults to the total
#'   number of predictors.
#' @param stop_at_first Logical. If `TRUE`, stop after the first minimal hit in
#'   forward search or backward layered exhaustive search.
#' @param missing Character string specifying how missing data are handled. One
#'   of `"complete"` or `"impute"`.
#' @param impute_args Optional list of imputation settings used when
#'   `missing = "impute"`.
#' @param scale_X Logical. If `TRUE`, encoded predictors are standardized using
#'   `scale()`.
#' @param tau_complete Numeric threshold on `delta_hat` used to declare complete
#'   separation. Default is `1e-6`.
#' @param eps_boundary Target margin delta for the severity LP. If `NULL`, it
#'   defaults to `epsilon`.
#' @param quasi_to_none_if Numeric value in `(0, 1]`. If
#'   `K_relax >= quasi_to_none_if * n`, quasi-complete separation is treated as
#'   no separation.
#' @param missing_scope Character string specifying whether missing-data handling
#'   is applied globally or separately by subset. One of `"global"` or
#'   `"subset"`. Default is `"global"`.
#' @param minimal_strategy Character string specifying the minimal-subset search
#'   strategy. One of `"auto"`, `"forward"`, or `"backward"`. Default is
#'   `"auto"`.
#' @param small_p_threshold Integer threshold used by
#'   `minimal_strategy = "auto"`. Forward search is used when
#'   `p <= small_p_threshold`; backward search is used otherwise. Default is 15.
#' @param backward_exhaustive Logical. If `TRUE` and
#'   `minimal_strategy = "backward"`, perform layered exhaustive backward search
#'   instead of beam-style backward search. Default is `FALSE`.
#' @param beam_width Positive integer controlling the number of best separating
#'   subsets retained at each layer in beam-style backward search. Ties at the
#'   cutoff are retained. Default is 10.
#' @param verbose Logical. If `TRUE`, print progress messages and stop reasons
#'   during minimal-subset search. Defaults to
#'   `getOption("latent_separation.verbose", FALSE)`.
#'
#' @return
#' In single-set mode, a list with fields:
#' \itemize{
#'   \item `type`: detected separation type or no-separation label.
#'   \item `satisfied`: logical indicator for whether the requested mode is satisfied.
#'   \item `available_types`: separation types available under the fitted diagnostics.
#'   \item `removed`: diagnostic row indices for quasi-complete separation.
#'   \item `message`: human-readable summary.
#'   \item `missing_info`: information about missing-data handling.
#'   \item `diagnostics`: LP diagnostics including `delta_hat`, `K_relax`, and `n`.
#' }
#'
#' In exhaustive mode, a named list of separating subsets. Each entry contains:
#' \itemize{
#'   \item `type`: separation type.
#'   \item `vars`: predictor names in the separating subset.
#'   \item `removed`: diagnostic row indices for quasi-complete separation.
#'   \item `missing_info`: information about missing-data handling.
#'   \item `diagnostics`: LP diagnostics.
#' }
#'
#' In minimal mode, a list with `minimal_subsets`. Each entry contains:
#' \itemize{
#'   \item `type`: separation type.
#'   \item `vars`: predictor names in the subset.
#'   \item `idx`: column indices of the subset.
#'   \item `removed`: diagnostic row indices for quasi-complete separation.
#'   \item `missing_info`: information about missing-data handling.
#'   \item `diagnostics`: LP diagnostics.
#'   \item strategy metadata such as `strategy`, `k`, `beam_width`,
#'     `tie_policy`, and `eval_count`, depending on the search strategy.
#' }
#'
#' For backward beam search, the returned object may also include `best_k`,
#' `beam_width`, `tie_policy`, and `eval_count`.
#'
#' For backward layered exhaustive search, the returned object may also include
#' `best_k`, `eval_count`, and `exhaustive_type`.
#'
#' @export
latent_separation <- function(
    y, X,
    test_combinations = FALSE,
    min_vars = 2,
    epsilon = 1e-5,
    only_perfect = FALSE,
    find_minimal = FALSE,
    mode = c("either","perfect","quasi"),
    max_vars = NULL,
    stop_at_first = FALSE,
    missing = c("complete","impute"),
    impute_args = list(),
    scale_X = TRUE,
    tau_complete = 1e-6,
    eps_boundary = NULL,
    quasi_to_none_if = 0.5,
    missing_scope = c("global","subset"),
    minimal_strategy = c("auto","forward","backward"),
    small_p_threshold = 15L,
    backward_exhaustive = FALSE,
    beam_width = 10L,
    verbose = getOption("latent_separation.verbose", FALSE)
) {
  mode <- match.arg(mode)
  missing <- match.arg(missing)
  missing_scope <- match.arg(missing_scope)
  minimal_strategy <- match.arg(minimal_strategy)
  verbose <- isTRUE(verbose)
  if (is.null(eps_boundary)) eps_boundary <- epsilon

  beam_width <- as.integer(beam_width)
  if (is.na(beam_width) || beam_width <= 0L) stop("beam_width must be a positive integer.")

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  vcat <- function(...) {
    if (verbose) message(sprintf(...))
    invisible(NULL)
  }

  # coerce X to data.frame
  if (is.matrix(X)) {
    X_raw <- as.data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (!is.data.frame(X)) {
    X_raw <- data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    X_raw <- data.frame(X, stringsAsFactors = FALSE, check.names = FALSE)
  }

  var_names <- colnames(X_raw)
  p <- ncol(X_raw)
  if (p == 0L) stop("X has no columns.")
  if (is.null(max_vars)) max_vars <- p
  min_vars <- max(1L, as.integer(min_vars))
  max_vars <- max(min_vars, as.integer(max_vars))
  if (min_vars > p) stop("min_vars cannot exceed number of columns in X.")
  if (max_vars > p) max_vars <- p

  y_raw <- y

  .as01 <- function(y){
    if (is.logical(y)) return(as.integer(y))
    if (is.factor(y))  return(as.integer(y == levels(y)[2L]))
    if (is.character(y)) {
      lev <- factor(y)
      if (nlevels(lev) != 2L) stop("Character outcome must have exactly 2 levels.")
      return(as.integer(lev == levels(lev)[2L]))
    }
    y <- as.integer(y)
    if (!all(y %in% c(0L,1L))) stop("Outcome y must be binary (0/1, factor/logical/2-level char).")
    y
  }

  # ============================================================
  # Stage A: max delta (tp/tm + L1 equality)
  # ============================================================
  .lp_max_delta <- function(y01, X_int, tau_complete = 1e-6) {
    ypm <- ifelse(y01 == 1L, 1, -1)
    n <- nrow(X_int)
    p_loc <- ncol(X_int) - 1L
    nvar <- 2 * p_loc + 2
    id_tp <- 1:p_loc
    id_tm <- (p_loc + 1):(2 * p_loc)
    id_b0 <- 2 * p_loc + 1
    id_de <- 2 * p_loc + 2

    A <- list()
    b <- c()
    d <- c()

    for (i in seq_len(n)) {
      row <- numeric(nvar)
      row[id_tp] <- -ypm[i] * X_int[i, -1, drop = FALSE]
      row[id_tm] <-  ypm[i] * X_int[i, -1, drop = FALSE]
      row[id_b0] <-  ypm[i]
      row[id_de] <-  1
      A[[length(A) + 1L]] <- row
      b <- c(b, 0)
      d <- c(d, "<=")
    }

    row_le <- numeric(nvar)
    row_ge <- numeric(nvar)
    row_le[c(id_tp, id_tm)] <-  1
    row_ge[c(id_tp, id_tm)] <- -1
    A[[length(A) + 1L]] <- row_le; b <- c(b,  1); d <- c(d, "<=")
    A[[length(A) + 1L]] <- row_ge; b <- c(b, -1); d <- c(d, "<=")

    rowd <- numeric(nvar)
    rowd[id_de] <- -1
    A[[length(A) + 1L]] <- rowd; b <- c(b, 0); d <- c(d, "<=")

    A <- do.call(rbind, A)
    obj <- numeric(nvar)
    obj[id_de] <- 1

    res <- lpSolve::lp("max", obj, A, d, b)
    if (res$status != 0) {
      return(list(status = "infeasible", delta_hat = NA_real_))
    }

    list(
      status = if (res$objval > tau_complete) "complete" else "borderline",
      delta_hat = res$objval,
      tp = res$solution[id_tp],
      tm = res$solution[id_tm],
      b0 = res$solution[id_b0]
    )
  }

  # ============================================================
  # Stage B: feasibility at delta = 0
  # ============================================================
  .lp_feasible_delta0 <- function(y01, X_int) {
    ypm <- ifelse(y01 == 1L, 1, -1)
    n <- nrow(X_int)
    p_loc <- ncol(X_int) - 1L
    nvar <- 2 * p_loc + 2
    id_tp <- 1:p_loc
    id_tm <- (p_loc + 1):(2 * p_loc)
    id_b0 <- 2 * p_loc + 1
    id_de <- 2 * p_loc + 2

    A <- list()
    b <- c()
    d <- c()

    for (i in seq_len(n)) {
      row <- numeric(nvar)
      row[id_tp] <- -ypm[i] * X_int[i, -1, drop = FALSE]
      row[id_tm] <-  ypm[i] * X_int[i, -1, drop = FALSE]
      row[id_b0] <-  ypm[i]
      row[id_de] <-  1
      A[[length(A) + 1L]] <- row
      b <- c(b, 0)
      d <- c(d, "<=")
    }

    row_le <- numeric(nvar)
    row_ge <- numeric(nvar)
    row_le[c(id_tp, id_tm)] <-  1
    row_ge[c(id_tp, id_tm)] <- -1
    A[[length(A) + 1L]] <- row_le; b <- c(b,  1); d <- c(d, "<=")
    A[[length(A) + 1L]] <- row_ge; b <- c(b, -1); d <- c(d, "<=")

    rowd <- numeric(nvar)
    rowd[id_de] <- -1
    A[[length(A) + 1L]] <- rowd; b <- c(b, 0); d <- c(d, "<=")

    row1 <- numeric(nvar)
    row2 <- numeric(nvar)
    row1[id_de] <-  1
    row2[id_de] <- -1
    A[[length(A) + 1L]] <- row1; b <- c(b, 0); d <- c(d, "<=")
    A[[length(A) + 1L]] <- row2; b <- c(b, 0); d <- c(d, "<=")

    A <- do.call(rbind, A)
    obj0 <- numeric(nvar)
    res <- lpSolve::lp("max", obj0, A, d, b)
    list(feasible = (res$status == 0))
  }

  # ============================================================
  # Severity LP: K_relax
  # ============================================================
  .lp_multivar_severity <- function(y01, X_int, eps = 1e-3) {
    ypm <- ifelse(y01 == 1L, 1, -1)
    n <- nrow(X_int)
    p_loc <- ncol(X_int) - 1L
    nvar <- 2 * p_loc + 1 + n
    id_tp <- 1:p_loc
    id_tm <- (p_loc + 1):(2 * p_loc)
    id_b0 <- 2 * p_loc + 1
    id_t <- (2 * p_loc + 2):(2 * p_loc + 1 + n)

    A <- list()
    b <- c()
    d <- c()

    for (i in seq_len(n)) {
      row <- numeric(nvar)
      row[id_tp] <- -ypm[i] * X_int[i, -1, drop = FALSE]
      row[id_tm] <-  ypm[i] * X_int[i, -1, drop = FALSE]
      row[id_b0] <-  ypm[i]
      row[id_t[i]] <- -1
      A[[length(A) + 1L]] <- row
      b <- c(b, -eps)
      d <- c(d, "<=")
    }

    row_le <- numeric(nvar)
    row_ge <- numeric(nvar)
    row_le[c(id_tp, id_tm)] <-  1
    row_ge[c(id_tp, id_tm)] <- -1
    A[[length(A) + 1L]] <- row_le; b <- c(b,  1); d <- c(d, "<=")
    A[[length(A) + 1L]] <- row_ge; b <- c(b, -1); d <- c(d, "<=")

    for (i in seq_len(n)) {
      row <- numeric(nvar)
      row[id_t[i]] <- -1
      A[[length(A) + 1L]] <- row
      b <- c(b, 0)
      d <- c(d, "<=")
    }

    A <- do.call(rbind, A)
    obj <- numeric(nvar)
    obj[id_t] <- -1
    res <- lpSolve::lp("max", obj, A, d, b)

    if (res$status != 0) return(list(K_relax = Inf, sum_t = NA_real_))

    sum_t <- sum(res$solution[id_t])
    list(K_relax = ceiling(sum_t / eps), sum_t = sum_t)
  }

  # ============================================================
  # Missing handling
  # ============================================================
  global_mh <- NULL
  if (missing_scope == "global") {
    global_mh <- DISCO:::.handle_missing_for_subset(
      y = y_raw,
      X = X_raw,
      method = missing,
      impute_args = impute_args
    )
  }

  # ============================================================
  # Main per-subset runner
  # ============================================================
  run_one <- function(cols, mode_local = mode, compute_removed = TRUE) {
    mh <- if (!is.null(global_mh)) {
      list(
        y = global_mh$y,
        X = global_mh$X[, cols, drop = FALSE],
        rows_used = global_mh$rows_used,
        params_used = global_mh$params_used
      )
    } else {
      DISCO:::.handle_missing_for_subset(
        y = y_raw,
        X = X_raw[, cols, drop = FALSE],
        method = missing,
        impute_args = impute_args
      )
    }

    y1 <- DISCO:::encode_outcome(mh$y)
    if (length(unique(y1)) < 2L) {
      return(list(
        hit = FALSE,
        type = "no separation problem",
        removed = NULL,
        missing_info = list(
          method = mh$params_used$method,
          params = mh$params_used,
          rows_used = mh$rows_used,
          n_used = length(mh$rows_used),
          scope = missing_scope
        ),
        message = "Fewer than two outcome classes after missing-data handling."
      ))
    }

    X1_mm <- DISCO:::encode_predictors_lp(mh$X)
    if (isTRUE(scale_X)) X1_mm <- scale(X1_mm)

    Xi <- cbind("(Intercept)" = 1, as.matrix(X1_mm))
    y01 <- .as01(y1)
    n <- nrow(Xi)

    # Stage A
    maxd <- .lp_max_delta(y01, Xi, tau_complete = tau_complete)
    if (identical(maxd$status, "complete")) {
      return(list(
        hit = (mode_local %in% c("either", "perfect")),
        type = "perfect separation",
        removed = NULL,
        diagnostics = list(
          perfect = TRUE,
          quasi = FALSE,
          delta_hat = maxd$delta_hat,
          K_relax = 0L,
          n = n
        ),
        missing_info = list(
          method = mh$params_used$method,
          params = mh$params_used,
          rows_used = mh$rows_used,
          n_used = length(mh$rows_used),
          scope = missing_scope
        ),
        message = "Perfect separation found on the given predictors."
      ))
    }

    # Stage B
    feas <- .lp_feasible_delta0(y01, Xi)
    sev <- .lp_multivar_severity(y01, Xi, eps = eps_boundary)
    Krel <- sev$K_relax

    removed <- integer(0)
    if (isTRUE(compute_removed) && feas$feasible) {
      for (i in seq_len(n)) {
        yi <- y01[-i]
        Xi_i <- Xi[-i, , drop = FALSE]
        maxd_i <- .lp_max_delta(yi, Xi_i, tau_complete = tau_complete)
        if (identical(maxd_i$status, "complete")) removed <- c(removed, i)
      }
    }

    quasi_flag <- (feas$feasible && is.finite(Krel) && (Krel < quasi_to_none_if * n))
    perfect_flag <- FALSE

    is_hit <- switch(
      mode_local,
      "perfect" = perfect_flag,
      "quasi"   = (quasi_flag && !perfect_flag),
      "either"  = (perfect_flag || quasi_flag)
    )

    if (!is_hit) {
      available <- character(0)
      if (perfect_flag) available <- c(available, "perfect")
      if (feas$feasible && (Krel < quasi_to_none_if * n)) available <- c(available, "quasi")

      return(list(
        hit = FALSE,
        type = "no separation problem",
        removed = NULL,
        diagnostics = list(
          perfect = perfect_flag,
          quasi = (feas$feasible && (Krel < quasi_to_none_if * n)),
          delta_hat = maxd$delta_hat,
          K_relax = Krel,
          n = n
        ),
        missing_info = list(
          method = mh$params_used$method,
          params = mh$params_used,
          rows_used = mh$rows_used,
          n_used = length(mh$rows_used),
          scope = missing_scope
        ),
        message = if (length(available)) {
          sprintf(
            "Separation exists (%s), but it does not satisfy mode = \"%s\".",
            paste(available, collapse = " & "),
            mode_local
          )
        } else {
          "No separation problem detected."
        }
      ))
    }

    list(
      hit = TRUE,
      type = if (perfect_flag) "perfect separation" else "quasi-complete separation",
      removed = if (perfect_flag) NULL else removed,
      diagnostics = list(
        perfect = perfect_flag,
        quasi = !perfect_flag,
        delta_hat = maxd$delta_hat,
        K_relax = Krel,
        n = n
      ),
      missing_info = list(
        method = mh$params_used$method,
        params = mh$params_used,
        rows_used = mh$rows_used,
        n_used = length(mh$rows_used),
        scope = missing_scope
      ),
      message = if (perfect_flag) {
        "Perfect separation found on the given predictors."
      } else {
        if (length(removed)) {
          sprintf(
            "Quasi-complete separation (latent). K_relax = %d; removing any of {%s} yields perfect (diagnostic).",
            Krel, paste(removed, collapse = ", ")
          )
        } else {
          sprintf("Quasi-complete separation (latent). K_relax = %d.", Krel)
        }
      }
    )
  }

  # ============================================================
  # Cache
  # ============================================================
  subset_cache <- new.env(parent = emptyenv())

  subset_key <- function(cols, mode_local, compute_removed) {
    paste0(
      paste(sort(cols), collapse = ","),
      "|mode=", mode_local,
      "|removed=", as.integer(isTRUE(compute_removed))
    )
  }

  run_one_cached <- function(cols, mode_local = mode, compute_removed = FALSE) {
    key <- subset_key(cols, mode_local, compute_removed)
    if (exists(key, envir = subset_cache, inherits = FALSE)) {
      return(get(key, envir = subset_cache, inherits = FALSE))
    }
    ans <- run_one(cols, mode_local = mode_local, compute_removed = compute_removed)
    assign(key, ans, envir = subset_cache)
    ans
  }

  candidate_score <- function(out) {
    dg <- out$diagnostics %||% list()
    Krel <- dg$K_relax %||% Inf
    dhat <- dg$delta_hat %||% -Inf

    if (identical(mode, "perfect")) {
      return(c(-dhat))
    }

    c(Krel, -dhat)
  }

  keep_all_ties_at_cutoff <- function(hit_records, beam_width) {
    if (!length(hit_records)) {
      return(integer(0))
    }

    score_mat <- do.call(
      rbind,
      lapply(hit_records, function(z) candidate_score(z$out))
    )

    if (is.vector(score_mat)) {
      score_mat <- matrix(score_mat, ncol = length(candidate_score(hit_records[[1]]$out)))
    }

    ord <- do.call(order, as.data.frame(score_mat))

    if (length(ord) <= beam_width) {
      return(ord)
    }

    cutoff_pos <- beam_width
    cutoff_score <- score_mat[ord[cutoff_pos], , drop = FALSE]

    same_as_cutoff <- apply(score_mat, 1, function(x) {
      isTRUE(all.equal(as.numeric(x), as.numeric(cutoff_score), tolerance = 0))
    })

    keep_pool <- ord[seq_len(cutoff_pos)]
    tied_extra <- setdiff(which(same_as_cutoff), keep_pool)

    keep_idx <- c(keep_pool, tied_extra)
    keep_idx <- ord[ord %in% keep_idx]
    keep_idx
  }

  # ============================================================
  # Single-set path
  # ============================================================
  if (!test_combinations && !find_minimal) {
    cols_all <- seq_len(p)
    res <- run_one(cols_all, if (only_perfect) "perfect" else mode)

    if (!isTRUE(res$hit)) {
      available <- character(0)
      if (isTRUE(res$diagnostics$perfect)) available <- c(available, "perfect")
      if (isTRUE(res$diagnostics$quasi))   available <- c(available, "quasi")

      if (length(available)) {
        return(list(
          type = "no separation (requested mode unsatisfied)",
          satisfied = FALSE,
          available_types = available,
          removed = NULL,
          message = if (!is.null(res$message)) {
            res$message
          } else {
            sprintf(
              "Separation exists (%s), but it does not satisfy mode = \"%s\".",
              paste(available, collapse = " & "),
              if (isTRUE(only_perfect)) "perfect" else mode
            )
          },
          missing_info = res$missing_info,
          diagnostics = res$diagnostics %||% list(
            K_relax = NA_real_,
            n = length(res$missing_info$rows_used %||% integer())
          ),
          vars_all = var_names,
          k_all = p
        ))
      }

      return(list(
        type = "no separation problem",
        satisfied = FALSE,
        available_types = character(0),
        message = if (!is.null(res$message)) res$message else "No separation problem detected.",
        missing_info = res$missing_info,
        diagnostics = res$diagnostics %||% list(
          K_relax = NA_real_,
          n = length(res$missing_info$rows_used %||% integer())
        ),
        vars_all = var_names,
        k_all = p
      ))
    }

    return(list(
      type = res$type,
      satisfied = TRUE,
      available_types = if (identical(res$type, "perfect separation")) c("perfect") else c("quasi"),
      removed = res$removed,
      message = res$message,
      missing_info = res$missing_info,
      diagnostics = res$diagnostics
    ))
  }

  # ============================================================
  # Full exhaustive enumeration by size
  # ============================================================
  if (test_combinations && !find_minimal) {
    results <- list()
    for (k in seq.int(min_vars, p)) {
      for (cols in combn(p, k, simplify = FALSE)) {
        nm <- paste(var_names[cols], collapse = "_")
        out <- run_one(cols, if (only_perfect) "perfect" else mode)
        if (out$hit) {
          results[[nm]] <- list(
            type = out$type,
            vars = var_names[cols],
            removed = out$removed,
            missing_info = out$missing_info,
            diagnostics = out$diagnostics
          )
        }
      }
    }
    if (only_perfect) {
      results <- Filter(function(z) identical(z$type, "perfect separation"), results)
    }
    return(results)
  }

  # ============================================================
  # Minimal-subset search
  # ============================================================
  if (find_minimal) {
    full_res <- run_one_cached(seq_len(p), mode, compute_removed = FALSE)
    if (!isTRUE(full_res$hit)) {
      vcat("minimal search: early stop because full predictor set is not a hit under mode = \"%s\".", mode)
      return(list(minimal_subsets = list()))
    }

    if (minimal_strategy == "auto") {
      small_p_threshold <- as.integer(small_p_threshold)
      if (is.na(small_p_threshold) || small_p_threshold <= 0L) small_p_threshold <- 15L
      minimal_strategy <- if (p <= small_p_threshold) "forward" else "backward"
    }

    eval_count <- 0L
    eval_limit <- getOption("latent_separation.eval_limit", Inf)
    show_progress <- verbose || isTRUE(getOption("latent_separation.show_progress", FALSE))
    progress_every <- as.integer(getOption("latent_separation.progress_every", 200L))
    if (is.na(progress_every) || progress_every <= 0L) progress_every <- 200L

    vcat(
      "minimal search: strategy=%s, p=%d, min_vars=%d, eval_limit=%s",
      minimal_strategy, p, min_vars,
      if (is.finite(eval_limit)) eval_limit else "Inf"
    )

    tick_progress <- function(hits_count, best_k) {
      if (!show_progress) return(invisible(NULL))
      if (eval_count %% progress_every == 0L) {
        message(sprintf(
          "minimal search: strategy=%s, evaluated=%d, hits=%d, best_k=%s",
          minimal_strategy, eval_count, hits_count,
          if (is.finite(best_k)) best_k else "Inf"
        ))
      }
      invisible(NULL)
    }

    # ------------------------------------------------------------
    # Forward minimal
    # ------------------------------------------------------------
    if (minimal_strategy == "forward") {
      minimal <- list()
      minimal_idx <- list()
      best_k <- Inf
      hits_count <- 0L

      contains_any_minimal <- function(cols) {
        if (!length(minimal_idx)) return(FALSE)
        for (m in minimal_idx) if (all(m %in% cols)) return(TRUE)
        FALSE
      }

      for (k in seq.int(min_vars, max_vars)) {
        if (k > best_k) break

        for (cols in combn(p, k, simplify = FALSE)) {
          if (contains_any_minimal(cols)) next

          eval_count <- eval_count + 1L
          if (eval_count > eval_limit) {
            warning("Reached evaluation limit in minimal-subset search; results may be incomplete.")
            vcat("minimal search: stopped at eval_limit=%s.", as.character(eval_limit))
            break
          }
          tick_progress(hits_count, best_k)

          out_fast <- run_one_cached(cols, mode, compute_removed = FALSE)
          if (!isTRUE(out_fast$hit)) next

          best_k <- min(best_k, k)
          hits_count <- hits_count + 1L

          out <- run_one_cached(cols, mode, compute_removed = TRUE)
          nm <- paste(var_names[cols], collapse = "_")
          minimal[[nm]] <- list(
            type = out$type,
            vars = var_names[cols],
            idx = cols,
            removed = out$removed,
            missing_info = out$missing_info,
            diagnostics = out$diagnostics,
            strategy = "forward",
            eval_count = eval_count
          )

          minimal_idx[[length(minimal_idx) + 1L]] <- cols

          if (isTRUE(stop_at_first)) return(list(minimal_subsets = minimal))
        }

        if (eval_count > eval_limit) break
      }

      return(list(minimal_subsets = minimal))
    }

    # ------------------------------------------------------------
    # Backward search
    # ------------------------------------------------------------
    if (minimal_strategy == "backward") {
      hits_count <- 0L

      frontier <- list(list(
        cols = seq_len(p),
        out = full_res
      ))

      best_layer <- frontier
      best_k <- p

      if (!isTRUE(backward_exhaustive)) {
        vcat(
          "backward beam: starting from full set, k=%d, beam_width=%d (soft cutoff with ties), mode=%s",
          p, beam_width, mode
        )

        repeat {
          current_k <- length(frontier[[1L]]$cols)
          if (current_k <= min_vars) {
            vcat("backward beam: stop because current_k=%d reached min_vars=%d.", current_k, min_vars)
            break
          }

          candidate_map <- new.env(parent = emptyenv())

          # layer idea: only generate from previous-layer hits
          for (node in frontier) {
            parent_cols <- node$cols
            for (j in parent_cols) {
              child <- setdiff(parent_cols, j)
              if (length(child) < min_vars) next
              if (length(child) > max_vars) next
              key <- paste(child, collapse = ",")
              if (!exists(key, envir = candidate_map, inherits = FALSE)) {
                assign(key, child, envir = candidate_map)
              }
            }
          }

          cand_keys <- ls(candidate_map, all.names = TRUE)
          if (!length(cand_keys)) {
            vcat("backward beam: no candidates generated from current frontier.")
            break
          }

          vcat(
            "backward beam: layer k=%d -> k=%d generated %d unique candidates before evaluation.",
            current_k, current_k - 1L, length(cand_keys)
          )

          hit_records <- list()

          # deduplicate first, then evaluate
          for (key in cand_keys) {
            eval_count <- eval_count + 1L
            if (eval_count > eval_limit) {
              warning("Reached evaluation limit in backward beam search; results may be incomplete.")
              vcat("backward beam: stopped at eval_limit=%s.", as.character(eval_limit))
              break
            }

            tick_progress(hits_count, best_k)

            cols <- get(key, envir = candidate_map, inherits = FALSE)
            out_fast <- run_one_cached(cols, mode, compute_removed = FALSE)

            if (!isTRUE(out_fast$hit)) next

            hits_count <- hits_count + 1L
            hit_records[[length(hit_records) + 1L]] <- list(
              key = key,
              cols = cols,
              out = out_fast
            )
          }

          if (eval_count > eval_limit) break

          if (!length(hit_records)) {
            vcat(
              "backward beam: no hits at next layer k=%d, stop. Final best_k=%d.",
              current_k - 1L, best_k
            )
            break
          }

          keep_idx <- keep_all_ties_at_cutoff(hit_records, beam_width)
          frontier <- hit_records[keep_idx]
          best_layer <- frontier
          best_k <- length(frontier[[1L]]$cols)

          vcat(
            "backward beam: retained %d hit(s) at k=%d for next frontier.",
            length(frontier), best_k
          )
        }

        minimal <- list()

        for (node in best_layer) {
          cols <- node$cols
          out <- run_one_cached(cols, mode, compute_removed = TRUE)
          nm <- paste(var_names[cols], collapse = "_")
          minimal[[nm]] <- list(
            type = out$type,
            vars = var_names[cols],
            idx = cols,
            removed = out$removed,
            missing_info = out$missing_info,
            diagnostics = out$diagnostics,
            strategy = "backward_beam",
            beam_width = beam_width,
            tie_policy = "keep_all_at_cutoff",
            k = length(cols),
            eval_count = eval_count
          )
        }

        vcat(
          "backward beam: returning %d hit(s) at best_k=%d.",
          length(minimal),
          if (length(minimal)) unique(vapply(minimal, function(z) z$k, integer(1)))[1L] else NA_integer_
        )

        return(list(
          minimal_subsets = minimal,
          best_k = if (length(minimal)) unique(vapply(minimal, function(z) z$k, integer(1)))[1L] else NA_integer_,
          beam_width = beam_width,
          tie_policy = "keep_all_at_cutoff",
          eval_count = eval_count
        ))
      }

      # ----------------------------------------------------------
      # Backward layered exhaustive
      # ----------------------------------------------------------
      vcat(
        "backward layered exhaustive: starting from full set, k=%d, mode=%s",
        p, mode
      )

      repeat {
        current_k <- length(frontier[[1L]]$cols)
        if (current_k <= min_vars) {
          vcat("backward layered exhaustive: stop because current_k=%d reached min_vars=%d.", current_k, min_vars)
          break
        }

        candidate_map <- new.env(parent = emptyenv())

        # layer idea: only generate from previous-layer hits
        # and keep ALL hits at next layer
        for (node in frontier) {
          parent_cols <- node$cols
          for (j in parent_cols) {
            child <- setdiff(parent_cols, j)
            if (length(child) < min_vars) next
            if (length(child) > max_vars) next
            key <- paste(child, collapse = ",")
            if (!exists(key, envir = candidate_map, inherits = FALSE)) {
              assign(key, child, envir = candidate_map)
            }
          }
        }

        cand_keys <- ls(candidate_map, all.names = TRUE)
        if (!length(cand_keys)) {
          vcat("backward layered exhaustive: no candidates generated from current frontier.")
          break
        }

        vcat(
          "backward layered exhaustive: layer k=%d -> k=%d generated %d unique candidates before evaluation.",
          current_k, current_k - 1L, length(cand_keys)
        )

        next_frontier <- list()

        for (key in cand_keys) {
          eval_count <- eval_count + 1L
          if (eval_count > eval_limit) {
            warning("Reached evaluation limit in backward layered exhaustive search; results may be incomplete.")
            vcat("backward layered exhaustive: stopped at eval_limit=%s.", as.character(eval_limit))
            break
          }

          tick_progress(hits_count, best_k)

          cols <- get(key, envir = candidate_map, inherits = FALSE)
          out_fast <- run_one_cached(cols, mode, compute_removed = FALSE)

          if (!isTRUE(out_fast$hit)) next

          hits_count <- hits_count + 1L
          next_frontier[[length(next_frontier) + 1L]] <- list(
            key = key,
            cols = cols,
            out = out_fast
          )

          if (isTRUE(stop_at_first)) {
            vcat("backward layered exhaustive: stop_at_first=TRUE, stopping at first hit in current layer.")
            break
          }
        }

        if (eval_count > eval_limit) break

        if (isTRUE(stop_at_first) && length(next_frontier)) {
          frontier <- next_frontier
          best_layer <- frontier
          best_k <- length(frontier[[1L]]$cols)
          break
        }

        if (!length(next_frontier)) {
          vcat(
            "backward layered exhaustive: no hits at next layer k=%d, stop. Final best_k=%d.",
            current_k - 1L, best_k
          )
          break
        }

        frontier <- next_frontier
        best_layer <- frontier
        best_k <- length(frontier[[1L]]$cols)

        vcat(
          "backward layered exhaustive: retained ALL %d hit(s) at k=%d for next frontier.",
          length(frontier), best_k
        )
      }

      minimal <- list()

      for (node in best_layer) {
        cols <- node$cols
        out <- run_one_cached(cols, mode, compute_removed = TRUE)
        nm <- paste(var_names[cols], collapse = "_")
        minimal[[nm]] <- list(
          type = out$type,
          vars = var_names[cols],
          idx = cols,
          removed = out$removed,
          missing_info = out$missing_info,
          diagnostics = out$diagnostics,
          strategy = "backward_layered_exhaustive",
          k = length(cols),
          eval_count = eval_count
        )
      }

      vcat(
        "backward layered exhaustive: returning %d hit(s) at best_k=%d.",
        length(minimal),
        if (length(minimal)) unique(vapply(minimal, function(z) z$k, integer(1)))[1L] else NA_integer_
      )

      return(list(
        minimal_subsets = minimal,
        best_k = if (length(minimal)) unique(vapply(minimal, function(z) z$k, integer(1)))[1L] else NA_integer_,
        eval_count = eval_count,
        exhaustive_type = "layered"
      ))
    }
  }

  stop("Unreachable state.")
}

