# ============================================================
# LP-based latent separation detection (full updated version)
#   Updates:
#   1) missing_scope = "global" option to fix samples across subsets
#   2) minimal_strategy = auto/forward/backward with backward options
# ============================================================

#' LP-based latent separation detection (with optional subset search)
#'
#' Detects latent (complete or quasi-complete) separation for a binary outcome
#' using linear programming in the multivariate space.
#'
#' Two-stage pipeline:
#'   Stage A: max-margin LP with L1 equality normalization to identify complete separation
#'   Stage B: severity LP to compute multivariate severity lower bound K_relax
#'
#' Notes:
#' - Quasi is treated as non-substantive and labeled no separation when
#'   K_relax / n is large (default quasi_to_none_if = 0.5) to reduce false positives.
#' - Legacy row-deletion diagnostic for quasi is kept (removed indices), but does not
#'   drive the quasi decision.
#'
#' Missingness handling across subsets:
#' - missing_scope = "global" (recommended): missing handling is done once on full X,
#'   then reused for every subset so rows are fixed across subset evaluations.
#' - missing_scope = "subset": missing handling is done per subset (legacy behavior).
#'
#' Minimal subset search strategies:
#' - minimal_strategy = "forward": enumerate increasing subset sizes (can be expensive)
#' - minimal_strategy = "backward": scalable greedy backward elimination by default
#' - minimal_strategy = "auto": forward if p <= small_p_threshold else backward
#' - backward_exhaustive = TRUE: enumerate subsets by decreasing size p-1, p-2, ...
#'   respecting eval_limit, can be expensive for moderate/large p.
#'
#' Runtime controls for minimal search:
#' - options(latent_separation.eval_limit = N)
#' - options(latent_separation.show_progress = TRUE/FALSE)
#' - options(latent_separation.progress_every = M)
#'
#' @param y Binary outcome (0/1 or 2-level factor/character/logical).
#' @param X Matrix or data frame of predictors (columns = variables).
#' @param test_combinations If TRUE, test every subset of size >= min_vars.
#' @param min_vars Minimum subset size to consider (default 2).
#' @param epsilon Numeric scalar > 0. Legacy margin parameter; used as default for eps_boundary.
#' @param only_perfect Legacy flag for exhaustive mode (kept for compatibility).
#' @param find_minimal If TRUE, return separating subsets based on strategy.
#' @param mode One of "either", "perfect", or "quasi".
#' @param max_vars Optional upper bound on subset size (default all the way to p).
#' @param stop_at_first If TRUE, stop after the first minimal hit (forward or backward_exhaustive).
#' @param missing How to treat missing data: "complete" or "impute".
#' @param impute_args Optional list of imputation settings when missing="impute".
#' @param scale_X Logical; if TRUE, standardize encoded predictors with scale().
#' @param tau_complete Threshold on delta_hat to declare complete separation. Default 1e-6.
#' @param eps_boundary Target margin delta for severity LP. If NULL, defaults to epsilon.
#' @param quasi_to_none_if Numeric in (0,1]. If K_relax >= rho * n, quasi treated as none.
#' @param missing_scope "global" or "subset". Default "global".
#' @param minimal_strategy "auto","forward","backward". Default "auto".
#' @param small_p_threshold Threshold used by minimal_strategy="auto". Default 15.
#' @param backward_exhaustive If TRUE and strategy is backward, enumerate by size. Default FALSE.
#'
#' @return
#' - Single-set mode: list with fields type, satisfied, available_types, removed, message,
#'   missing_info, diagnostics.
#' - Exhaustive mode: named list of hits with type, vars, removed, missing_info, diagnostics.
#' - Minimal mode: list with $minimal_subsets, each entry containing type, vars, idx, removed,
#'   missing_info, diagnostics plus strategy metadata.
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
    backward_exhaustive = FALSE
) {
  mode <- match.arg(mode)
  missing <- match.arg(missing)
  missing_scope <- match.arg(missing_scope)
  minimal_strategy <- match.arg(minimal_strategy)
  if (is.null(eps_boundary)) eps_boundary <- epsilon

  `%||%` <- function(a, b) if (!is.null(a)) a else b

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

  y_raw <- y

  .as01 <- function(y){
    if (is.logical(y)) return(as.integer(y))
    if (is.factor(y))  return(as.integer(y == levels(y)[2L]))
    if (is.character(y)) {
      lev <- factor(y)
      if (length(lev) != 2L) stop("Character outcome must have exactly 2 levels.")
      return(as.integer(y == lev[2L]))
    }
    y <- as.integer(y)
    if (!all(y %in% c(0L,1L))) stop("Outcome y must be binary (0/1, factor/logical/2-level char).")
    y
  }

  # ---- user must supply these utilities (already in your project) ----
  # .handle_missing_for_subset(y, X, method, impute_args)
  # encode_outcome(y)
  # encode_predictors_lp(X)
  # -------------------------------------------------------------------

  # Stage A: max delta (tp/tm + L1 equality)
  .lp_max_delta <- function(y01, X_int, tau_complete = 1e-6) {
    ypm <- ifelse(y01==1L, 1, -1)
    n <- nrow(X_int); p <- ncol(X_int) - 1L
    nvar <- 2*p + 2
    id_tp <- 1:p; id_tm <- (p+1):(2*p); id_b0 <- 2*p+1; id_de <- 2*p+2

    A <- list(); b <- c(); d <- c()
    for (i in seq_len(n)){
      row <- numeric(nvar)
      row[id_tp] <- - ypm[i] * X_int[i, -1, drop=FALSE]
      row[id_tm] <- + ypm[i] * X_int[i, -1, drop=FALSE]
      row[id_b0] <- + ypm[i]
      row[id_de] <- + 1
      A[[length(A)+1]] <- row; b <- c(b, 0); d <- c(d, "<=")
    }
    row_le <- numeric(nvar); row_le[c(id_tp,id_tm)] <-  1
    row_ge <- numeric(nvar); row_ge[c(id_tp,id_tm)] <- -1
    A[[length(A)+1]] <- row_le; b <- c(b,  1); d <- c(d, "<=")
    A[[length(A)+1]] <- row_ge; b <- c(b, -1); d <- c(d, "<=")
    rowd <- numeric(nvar); rowd[id_de] <- -1
    A[[length(A)+1]] <- rowd; b <- c(b, 0); d <- c(d, "<=")

    A <- do.call(rbind, A)
    obj <- numeric(nvar); obj[id_de] <- 1
    res <- lpSolve::lp("max", obj, A, d, b)
    if (res$status != 0) return(list(status="infeasible", delta_hat=NA_real_))
    list(
      status = if (res$objval > tau_complete) "complete" else "borderline",
      delta_hat = res$objval,
      tp = res$solution[id_tp],
      tm = res$solution[id_tm],
      b0 = res$solution[id_b0]
    )
  }

  # Stage B: feasibility at delta=0
  .lp_feasible_delta0 <- function(y01, X_int) {
    ypm <- ifelse(y01==1L, 1, -1)
    n <- nrow(X_int); p <- ncol(X_int) - 1L
    nvar <- 2*p + 2
    id_tp <- 1:p; id_tm <- (p+1):(2*p); id_b0 <- 2*p+1; id_de <- 2*p+2

    A <- list(); b <- c(); d <- c()
    for (i in seq_len(n)){
      row <- numeric(nvar)
      row[id_tp] <- - ypm[i] * X_int[i, -1, drop=FALSE]
      row[id_tm] <- + ypm[i] * X_int[i, -1, drop=FALSE]
      row[id_b0] <- + ypm[i]
      row[id_de] <- + 1
      A[[length(A)+1]] <- row; b <- c(b, 0); d <- c(d, "<=")
    }
    row_le <- numeric(nvar); row_le[c(id_tp,id_tm)] <-  1
    row_ge <- numeric(nvar); row_ge[c(id_tp,id_tm)] <- -1
    A[[length(A)+1]] <- row_le; b <- c(b,  1); d <- c(d, "<=")
    A[[length(A)+1]] <- row_ge; b <- c(b, -1); d <- c(d, "<=")
    rowd <- numeric(nvar); rowd[id_de] <- -1
    A[[length(A)+1]] <- rowd; b <- c(b, 0); d <- c(d, "<=")
    row1 <- numeric(nvar); row1[id_de] <- +1
    row2 <- numeric(nvar); row2[id_de] <- -1
    A[[length(A)+1]] <- row1; b <- c(b, 0); d <- c(d, "<=")
    A[[length(A)+1]] <- row2; b <- c(b, 0); d <- c(d, "<=")

    A <- do.call(rbind, A)
    obj0 <- numeric(nvar)
    res <- lpSolve::lp("max", obj0, A, d, b)
    list(feasible = (res$status == 0))
  }

  # Severity LP: K_relax
  .lp_multivar_severity <- function(y01, X_int, eps = 1e-3) {
    ypm <- ifelse(y01==1L, 1, -1)
    n <- nrow(X_int); p <- ncol(X_int) - 1L
    nvar <- 2*p + 1 + n
    id_tp <- 1:p; id_tm <- (p+1):(2*p); id_b0 <- 2*p+1; id_t <- (2*p+2):(2*p+1+n)

    A <- list(); b <- c(); d <- c()
    for (i in seq_len(n)){
      row <- numeric(nvar)
      row[id_tp] <- - ypm[i] * X_int[i, -1, drop=FALSE]
      row[id_tm] <- + ypm[i] * X_int[i, -1, drop=FALSE]
      row[id_b0] <- + ypm[i]
      row[id_t[i]] <- -1
      A[[length(A)+1]] <- row; b <- c(b, -eps); d <- c(d, "<=")
    }
    row_le <- numeric(nvar); row_le[c(id_tp,id_tm)] <-  1
    row_ge <- numeric(nvar); row_ge[c(id_tp,id_tm)] <- -1
    A[[length(A)+1]] <- row_le; b <- c(b,  1); d <- c(d, "<=")
    A[[length(A)+1]] <- row_ge; b <- c(b, -1); d <- c(d, "<=")
    for (i in seq_len(n)){
      row <- numeric(nvar); row[id_t[i]] <- -1
      A[[length(A)+1]] <- row; b <- c(b, 0); d <- c(d, "<=")
    }
    A <- do.call(rbind, A)
    obj <- numeric(nvar); obj[id_t] <- -1
    res <- lpSolve::lp("max", obj, A, d, b)
    if (res$status != 0) return(list(K_relax=Inf, sum_t=NA_real_))
    sum_t <- sum(res$solution[id_t])
    list(K_relax = ceiling(sum_t / eps), sum_t = sum_t)
  }

  # ---- missing handling: global preprocessing if requested ----
  global_mh <- NULL
  if (missing_scope == "global") {
    global_mh <- .handle_missing_for_subset(
      y = y_raw, X = X_raw,
      method = missing, impute_args = impute_args
    )
  }

  # main per-subset runner
  run_one <- function(cols, mode_local = mode, compute_removed = TRUE) {
    mh <- if (!is.null(global_mh)) {
      list(
        y = global_mh$y,
        X = global_mh$X[, cols, drop = FALSE],
        rows_used = global_mh$rows_used,
        params_used = global_mh$params_used
      )
    } else {
      .handle_missing_for_subset(
        y = y_raw,
        X = X_raw[, cols, drop = FALSE],
        method = missing, impute_args = impute_args
      )
    }

    y1 <- encode_outcome(mh$y)
    if (length(unique(y1)) < 2L) {
      return(list(
        hit = FALSE, type = "no separation problem",
        removed = NULL,
        missing_info = list(method = mh$params_used$method, params = mh$params_used,
                            rows_used = mh$rows_used, n_used = length(mh$rows_used),
                            scope = missing_scope),
        message = "Fewer than two outcome classes after missing-data handling."
      ))
    }

    X1_mm <- encode_predictors_lp(mh$X)
    if (isTRUE(scale_X)) X1_mm <- scale(X1_mm)

    Xi <- cbind("(Intercept)" = 1, as.matrix(X1_mm))
    y01 <- .as01(y1)
    n <- nrow(Xi)

    # Stage A
    maxd <- .lp_max_delta(y01, Xi, tau_complete = tau_complete)
    if (identical(maxd$status, "complete")) {
      return(list(
        hit = (mode_local %in% c("either","perfect")),
        type = "perfect separation",
        removed = NULL,
        diagnostics = list(perfect=TRUE, quasi=FALSE, delta_hat=maxd$delta_hat,
                           K_relax=0L, n=n),
        missing_info = list(method = mh$params_used$method, params = mh$params_used,
                            rows_used = mh$rows_used, n_used = length(mh$rows_used),
                            scope = missing_scope),
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

    is_hit <- switch(mode_local,
                     "perfect" = perfect_flag,
                     "quasi"   = (quasi_flag && !perfect_flag),
                     "either"  = (perfect_flag || quasi_flag))

    if (!is_hit) {
      available <- character(0)
      if (perfect_flag) available <- c(available, "perfect")
      if (feas$feasible && (Krel < quasi_to_none_if * n)) available <- c(available, "quasi")

      return(list(
        hit = FALSE,
        type = "no separation problem",
        removed = NULL,
        diagnostics = list(perfect=perfect_flag,
                           quasi=(feas$feasible && (Krel < quasi_to_none_if * n)),
                           delta_hat=maxd$delta_hat, K_relax=Krel, n=n),
        missing_info = list(method = mh$params_used$method, params = mh$params_used,
                            rows_used = mh$rows_used, n_used = length(mh$rows_used),
                            scope = missing_scope),
        message = if (length(available)) {
          sprintf("Separation exists (%s), but it does not satisfy mode = \"%s\".",
                  paste(available, collapse = " & "), mode_local)
        } else "No separation problem detected."
      ))
    }

    list(
      hit = TRUE,
      type = if (perfect_flag) "perfect separation" else "quasi-complete separation",
      removed = if (perfect_flag) NULL else removed,
      diagnostics = list(perfect=perfect_flag, quasi=!perfect_flag,
                         delta_hat=maxd$delta_hat, K_relax=Krel, n=n),
      missing_info = list(method = mh$params_used$method, params = mh$params_used,
                          rows_used = mh$rows_used, n_used = length(mh$rows_used),
                          scope = missing_scope),
      message = if (perfect_flag) {
        "Perfect separation found on the given predictors."
      } else {
        if (length(removed)) {
          sprintf("Quasi-complete separation (latent). K_relax = %d; removing any of {%s} yields perfect (diagnostic).",
                  Krel, paste(removed, collapse = ", "))
        } else {
          sprintf("Quasi-complete separation (latent). K_relax = %d.", Krel)
        }
      }
    )
  }

  # === Single-set path ===
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
          message = if (!is.null(res$message)) res$message else
            sprintf("Separation exists (%s), but it does not satisfy mode = \"%s\".",
                    paste(available, collapse = " & "),
                    if (isTRUE(only_perfect)) "perfect" else mode),
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

  # === Exhaustive enumeration ===
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
    if (only_perfect) results <- Filter(function(z) identical(z$type, "perfect separation"), results)
    return(results)
  }

  # === Minimal-subset search ===
  if (find_minimal) {
    # Early global check: if full set is not a hit under mode, no subset can be a hit under that mode
    full_res <- run_one(seq_len(p), mode, compute_removed = FALSE)
    if (!isTRUE(full_res$hit)) {
      if (isTRUE(getOption("latent_separation.show_progress", FALSE))) {
        message("minimal search: stopped early because full predictor set is not a hit under the requested mode; returning empty.")
      }
      return(list(minimal_subsets = list()))
    }

    if (minimal_strategy == "auto") {
      small_p_threshold <- as.integer(small_p_threshold)
      if (is.na(small_p_threshold) || small_p_threshold <= 0L) small_p_threshold <- 15L
      minimal_strategy <- if (p <= small_p_threshold) "forward" else "backward"
    }

    eval_count <- 0L
    eval_limit <- getOption("latent_separation.eval_limit", Inf)
    show_progress <- isTRUE(getOption("latent_separation.show_progress", FALSE))
    progress_every <- as.integer(getOption("latent_separation.progress_every", 200L))
    if (is.na(progress_every) || progress_every <= 0L) progress_every <- 200L

    tick_progress <- function(hits_count, best_k) {
      if (!show_progress) return(invisible(NULL))
      if (eval_count %% progress_every == 0L) {
        message(sprintf("minimal search: strategy=%s, evaluated=%d, hits=%d, best_k=%s",
                        minimal_strategy, eval_count, hits_count,
                        if (is.finite(best_k)) best_k else "Inf"))
      }
      invisible(NULL)
    }

    # ---- Forward minimal (your current behavior) ----
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
            break
          }
          tick_progress(hits_count, best_k)

          out_fast <- run_one(cols, mode, compute_removed = FALSE)
          if (!isTRUE(out_fast$hit)) next

          best_k <- min(best_k, k)
          hits_count <- hits_count + 1L

          out <- run_one(cols, mode, compute_removed = TRUE)
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

    # ---- Backward minimal ----
    # Default: greedy backward elimination that returns one locally minimal separating set
    # Optional: backward_exhaustive enumerates subsets by decreasing size p-1, p-2, ...
    if (!isTRUE(backward_exhaustive)) {
      current <- seq_len(p)
      improved <- TRUE
      hits_count <- 0L

      while (improved) {
        improved <- FALSE
        for (j in current) {
          candidate <- setdiff(current, j)
          if (length(candidate) < min_vars) next

          eval_count <- eval_count + 1L
          if (eval_count > eval_limit) {
            warning("Reached evaluation limit in backward greedy search; result may be incomplete.")
            improved <- FALSE
            break
          }
          tick_progress(hits_count, length(current))

          out <- run_one(candidate, mode, compute_removed = FALSE)
          if (isTRUE(out$hit)) {
            current <- candidate
            improved <- TRUE
            hits_count <- hits_count + 1L
            break
          }
        }
      }

      # local minimality check: removing any single var should break hit
      minimal_by_1 <- TRUE
      for (j in current) {
        cand2 <- setdiff(current, j)
        if (length(cand2) < min_vars) next

        eval_count <- eval_count + 1L
        if (eval_count > eval_limit) {
          minimal_by_1 <- NA
          warning("Reached evaluation limit during minimality verification; minimal_by_1 is NA.")
          break
        }
        out2 <- run_one(cand2, mode, compute_removed = FALSE)
        if (isTRUE(out2$hit)) {
          minimal_by_1 <- FALSE
          break
        }
      }

      final_out <- run_one(current, mode, compute_removed = TRUE)
      nm <- paste(var_names[current], collapse = "_")
      minimal <- list()
      minimal[[nm]] <- list(
        type = final_out$type,
        vars = var_names[current],
        idx = current,
        removed = final_out$removed,
        missing_info = final_out$missing_info,
        diagnostics = final_out$diagnostics,
        strategy = "backward_greedy",
        minimal_by_1 = minimal_by_1,
        eval_count = eval_count
      )
      return(list(minimal_subsets = minimal))
    }

    # ---- Backward exhaustive (return all hits at minimal size) ----
    minimal <- list()
    best_hits <- list()
    best_k <- NA_integer_
    hits_count <- 0L

    for (k in seq.int(p - 1L, min_vars, by = -1L)) {
      hits_at_k <- list()

      for (cols in combn(p, k, simplify = FALSE)) {
        eval_count <- eval_count + 1L
        if (eval_count > eval_limit) {
          warning("Reached evaluation limit in backward exhaustive search; results may be incomplete.")
          break
        }
        tick_progress(hits_count, k)

        out_fast <- run_one(cols, mode, compute_removed = FALSE)
        if (isTRUE(out_fast$hit)) {
          hits_count <- hits_count + 1L
          nm <- paste(var_names[cols], collapse = "_")
          hits_at_k[[nm]] <- cols
          if (isTRUE(stop_at_first)) break
        }
      }

      if (eval_count > eval_limit) break

      if (length(hits_at_k)) {
        best_hits <- hits_at_k
        best_k <- k
        if (isTRUE(stop_at_first)) break
        # keep going smaller to see if we can still hit with fewer vars
      } else {
        # once we found hits at some larger k, and now no hits at smaller k,
        # the minimal size is best_k and we can stop
        if (!is.na(best_k)) break
      }
    }

    if (!length(best_hits)) return(list(minimal_subsets = list(), best_k = NA_integer_))

    # Return ALL problematic combinations at size best_k
    for (nm in names(best_hits)) {
      cols <- best_hits[[nm]]
      out <- run_one(cols, mode, compute_removed = TRUE)
      minimal[[nm]] <- list(
        type = out$type,
        vars = var_names[cols],
        idx = cols,
        removed = out$removed,
        missing_info = out$missing_info,
        diagnostics = out$diagnostics,
        strategy = "backward_exhaustive_all_at_min_k",
        k = best_k,
        eval_count = eval_count
      )
    }

    return(list(minimal_subsets = minimal, best_k = best_k))

  }

  stop("Unreachable state.")
}
