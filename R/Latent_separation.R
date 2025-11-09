# ============================================================
# LP-based latent separation detection (updated: L1-eq + K_relax)
# ============================================================
#' LP-based latent separation detection (with optional minimal-subset search)
#'
#' This function detects latent (complete or quasi-complete) separation for a
#' binary outcome using **linear programming** in the *multivariate* space.
#' It upgrades the classic feasibility-only LP check to a **two-stage** pipeline:
#' (i) a **max-margin LP** with positive/negative split and an L1 *equality*
#' normalization to identify **complete separation**; and (ii) a **severity LP**
#' that computes a multivariate severity lower bound \eqn{K_relax}, i.e.,
#' the minimal “boundary tie budget” (LP relaxation) needed to lift the **smallest
#' margin** to a target \eqn{\delta} on the standardized scale. When
#' \eqn{K_relax/n} is large (by default \eqn{\ge \rho=0.9}), quasi is treated as
#' **non-substantive** and labeled **no separation** to avoid false positives.
#'
#' The function keeps your previous interface and search utilities:
#' - A single-set check (default).
#' - Exhaustive subset testing (`test_combinations = TRUE`).
#' - Minimal-subset search with pruning (`find_minimal = TRUE`).
#'
#' It also preserves the legacy **row-deletion** diagnostic for quasi
#' (indices of observations whose removal yields perfect separation),
#' reported in `removed`, but **does not** drive the quasi decision alone.
#'
#' @param y Binary outcome (0/1 or 2-level factor/character/logical).
#' @param X Matrix or data frame of predictors (columns = variables).
#' @param test_combinations If `TRUE`, test every subset of size ≥ `min_vars`
#'   (legacy exhaustive mode).
#' @param min_vars Minimum subset size to consider (default 2).
#' @param epsilon Numeric scalar > 0. Legacy margin parameter for the old LP
#'   feasibility; here it is used as the **default** for `eps_boundary`
#'   (the target margin in the severity LP). Default `1e-5`.
#' @param only_perfect Legacy flag for exhaustive mode (kept for compatibility).
#' @param find_minimal If `TRUE`, return inclusion-minimal separating subsets
#'   using pruning.
#' @param mode One of `"either"`, `"perfect"`, or `"quasi"` – which kind of
#'   separation to regard as a **hit**.
#' @param max_vars Optional upper bound on subset size (default: all the way to p).
#' @param stop_at_first If `TRUE` and `find_minimal=TRUE`, stop after the first
#'   minimal hit.
#' @param missing How to treat missing data: `"complete"` (drop rows with NA in
#'   the subset) or `"impute"` (impute predictors; never impute outcome).
#' @param impute_args Optional list of imputation settings when
#'   `missing = "impute"`. Recognized keys:
#'   \itemize{
#'     \item \code{numeric_method}: `"median"` (default) or `"mean"`
#'     \item \code{categorical_method}: `"mode"` (default) or `"missing"`
#'     \item \code{logical_method}: `"mode"` (default) or `"missing"`
#'     \item \code{custom_fn}: function(data.frame) -> imputed data.frame
#'   }
#' @param scale_X Logical; if `TRUE`, standardize encoded predictors with base
#'   \code{scale()} (default `FALSE`). Recommended `TRUE` when comparing
#'   `K_relax` across datasets.
#' @param tau_complete Numeric scalar > 0. Threshold on the **max-margin** value
#'   to declare **complete separation** (Stage A). Default `1e-6` (after scaling).
#' @param eps_boundary Numeric scalar > 0. Target margin \eqn{\delta} for the
#'   **severity LP** used to compute \eqn{K_relax}. If `NULL`, defaults to
#'   `epsilon`. Typical values on standardized features: `1e-3`–`1e-4`.
#' @param quasi_to_none_if Numeric in (0,1]. If \eqn{K_relax \ge \rho n}, then
#'   quasi is treated as **none** (non-substantive). Default `0.9`.
#'
#' @return
#' - **Single-set mode** (default): a list with fields
#'   \itemize{
#'     \item `type`: `"perfect separation"`, `"quasi-complete separation"`, or
#'       `"no separation problem"`.
#'     \item `satisfied`: logical, whether the requested `mode` is satisfied.
#'     \item `available_types`: which types exist regardless of `mode`.
#'     \item `removed`: (quasi only) integer indices of observations whose removal
#'           yields perfect (diagnostic; not used for the decision).
#'     \item `message`: human-readable message.
#'     \item `missing_info`: list with missing-data handling metadata.
#'     \item `diagnostics`: list with `delta_hat` (max-margin), `K_relax`, and `n`.
#'   }
#'
#' - **Exhaustive mode** (`test_combinations=TRUE`): a named list of hits (by
#'   concatenated variable names), each element containing
#'   `{type, vars, removed, missing_info, diagnostics}`.
#'
#' - **Minimal-subset mode** (`find_minimal=TRUE`): a list with
#'   `$minimal_subsets`, each entry `{type, vars, idx, removed, missing_info, diagnostics}`.
#'
#' @section Method (outline):
#' Stage A (**max-margin LP**): maximize \eqn{\delta} subject to
#' \eqn{y_i(\beta_0 + x_i^\top\beta) - \delta \ge 0} with
#' \eqn{\beta = s^+ - s^-}, \eqn{s^\pm \ge 0}, and \eqn{\| \beta \|_1 = \mathbf{1}^\top(s^+ + s^-) = 1}.
#' If \eqn{\hat\delta > \tau_\text{complete}}, declare **complete separation**.
#'
#' Stage B (**delta=0 feasibility**): check feasibility with \eqn{\delta=0} under
#' the same normalization. If feasible, compute **multivariate severity**
#' \eqn{K_\text{relax} = \lceil \sum_i t_i^\* / \delta \rceil} from the LP
#' \eqn{\min \sum_i t_i} s.t. \eqn{y_i(\beta_0 + x_i^\top\beta) + t_i \ge \delta},
#' \eqn{t_i \ge 0}, \eqn{\| \beta \|_1 = 1}. If \eqn{K_\text{relax} \ge \rho n}
#' (with \eqn{\rho=\texttt{quasi\_to\_none\_if}}), treat quasi as **none**.
#'
#' @section Practical tips:
#' - Use `scale_X = TRUE` (recommended) to stabilize thresholds and comparisons.
#' - Sensitivity-check `eps_boundary` (e.g., `1e-4`–`1e-2`); `K_relax` increases
#'   with larger target margins.
#' - `quasi_to_none_if` controls how strict you are in down-weighting weak quasi.
#'
#' @examples
#' # Simple no-separation example
#' X1 <- c(0, -10, 2, -9, 5)
#' Y1 <- c(0,   0, 0,  1, 1)
#' res <- latent_separation(Y1, cbind(X1), scale_X = TRUE)
#' str(res)
#' res$type; res$diagnostics
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
    scale_X = FALSE,
    # --- new knobs (backward compatible defaults) ---
    tau_complete = 1e-6,
    eps_boundary = NULL,
    quasi_to_none_if = 0.9
) {
  mode <- match.arg(mode)
  missing <- match.arg(missing)
  if (is.null(eps_boundary)) eps_boundary <- epsilon

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

  y_raw <- y  # keep original

  .as01 <- function(y){
    if (is.logical(y)) return(as.integer(y))
    if (is.factor(y))  return(as.integer(y == levels(y)[2L]))
    if (is.character(y)) {
      # 两水平字符，取第二个水平为 1
      lev <- unique(y)
      if (length(lev) != 2L) stop("Character outcome must have exactly 2 levels.")
      return(as.integer(y == lev[2L]))
    }
    y <- as.integer(y)
    if (!all(y %in% c(0L,1L))) stop("Outcome y must be binary (0/1, factor/logical/2-level char).")
    y
  }


  # --- internal helpers (rely on your own existing implementations) -----------
  # You already have:
  # - .handle_missing_for_subset(y, X, method, impute_args)
  # - encode_outcome(y)
  # - encode_predictors_lp(X)
  # The function below assumes they are available in your package/environment.

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

  # main per-subset runner (uses your missing+encoding utilities)
  run_one <- function(cols, mode_local = mode) {
    mh <- .handle_missing_for_subset(y = y_raw, X = X_raw[, cols, drop = FALSE],
                                     method = missing, impute_args = impute_args)

    y1 <- encode_outcome(mh$y)
    if (length(unique(y1)) < 2L) {
      return(list(
        hit = FALSE, type = "no separation problem",
        removed = NULL,
        missing_info = list(method = mh$params_used$method, params = mh$params_used,
                            rows_used = mh$rows_used, n_used = length(mh$rows_used)),
        message = "Fewer than two outcome classes after missing-data handling."
      ))
    }

    X1_mm <- encode_predictors_lp(mh$X)
    if (isTRUE(scale_X)) X1_mm <- scale(X1_mm)

    # add intercept
    Xi <- cbind("(Intercept)" = 1, as.matrix(X1_mm))
    y01 <- .as01(y1)
    n <- nrow(Xi)

    # Stage A: max delta
    maxd <- .lp_max_delta(y01, Xi, tau_complete = tau_complete)
    if (identical(maxd$status, "complete")) {
      return(list(
        hit = (mode_local %in% c("either","perfect")),
        type = "perfect separation",
        removed = NULL,
        diagnostics = list(perfect=TRUE, quasi=FALSE, delta_hat=maxd$delta_hat,
                           K_relax=0L, n=n),
        missing_info = list(method = mh$params_used$method, params = mh$params_used,
                            rows_used = mh$rows_used, n_used = length(mh$rows_used)),
        message = "Perfect separation found on the given predictors."
      ))
    }

    # Stage B: feasibility at delta=0
    feas <- .lp_feasible_delta0(y01, Xi)
    # severity
    sev <- .lp_multivar_severity(y01, Xi, eps = eps_boundary)
    Krel <- sev$K_relax

    # legacy row-deletion diagnostic for quasi
    removed <- integer(0)
    if (feas$feasible) {
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
                            rows_used = mh$rows_used, n_used = length(mh$rows_used)),
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
                          rows_used = mh$rows_used, n_used = length(mh$rows_used)),
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

  # === Single-set path =========================================================
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
          missing_info = res$missing_info
        ))
      }

      return(list(
        type = "no separation problem",
        satisfied = FALSE,
        available_types = character(0),
        message = if (!is.null(res$message)) res$message else
          "No separation problem detected.",
        missing_info = res$missing_info
      ))
    }

    return(list(
      type = res$type,
      satisfied = TRUE,
      available_types = if (identical(res$type, "perfect separation"))
        c("perfect") else c("quasi"),
      removed = res$removed,
      message = res$message,
      missing_info = res$missing_info,
      diagnostics = res$diagnostics
    ))
  }

  # === Exhaustive enumeration ==================================================
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

  # === Minimal-subset search with pruning =====================================
  minimal <- list()
  minimal_idx <- list()
  contains_any_minimal <- function(cols) {
    if (!length(minimal_idx)) return(FALSE)
    for (m in minimal_idx) if (all(m %in% cols)) return(TRUE)
    FALSE
  }

  for (k in seq.int(min_vars, max_vars)) {
    cand <- combn(p, k, simplify = FALSE)
    cand <- Filter(function(cols) !contains_any_minimal(cols), cand)
    if (!length(cand)) next
    for (cols in cand) {
      out <- run_one(cols, mode)
      if (out$hit) {
        nm <- paste(var_names[cols], collapse = "_")
        minimal[[nm]] <- list(
          type = out$type,
          vars = var_names[cols],
          idx  = cols,
          removed = out$removed,
          missing_info = out$missing_info,
          diagnostics = out$diagnostics
        )
        minimal_idx[[length(minimal_idx) + 1L]] <- cols
        if (isTRUE(stop_at_first)) return(list(minimal_subsets = minimal))
      }
    }
  }

  list(minimal_subsets = minimal)
}
