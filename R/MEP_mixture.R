#' Severity-Adaptive MEP for Mixture (multi-predictor) Logistic
#'
#' Fits a multi-predictor logistic model with a **Multivariate Exponential Power (MEP)**
#' prior using a simple Random-Walk Metropolis–Hastings (RW–MH) sampler, where the *slope
#' prior scales* are **anchored by univariate DISCO severities**. A small grid over the
#' intercept prior mean, a global multiplier on slope scales, and the EP shape \eqn{\kappa}
#' is explored; one run is selected via acceptance window, posterior predictive agreement,
#' and—when available—closeness to GLM coefficient ratios relative to a reference predictor.
#'
#' **Design encoding.** Predictors in `X` are expanded with
#' `model.matrix(~ ., data = X)` (intercept dropped for slopes). Numeric columns remain
#' one column each. Factor columns are expanded to treatment-contrast dummies (baseline is
#' the first level). The RW–MH is run on the *standardized* encoded design (z-scored
#' columns). Optionally, coefficients can be back-transformed to the original (unscaled)
#' encoded columns via `transform_back = "logit"|"SAS"|"Long"`.
#'
#' @section What this function does:
#' \itemize{
#'   \item Applies a single missing-data policy to `(y, X)` **once** (complete-case or external imputation).
#'   \item For each *original* predictor in `X`, computes a univariate **severity**
#'         using `DISCO::uni_separation()` on the same rows used for modeling.
#'         Numeric predictors are **z-scored for the severity computation** (factors unchanged).
#'   \item Maps severity \eqn{s \in [0,1]} to an anchor slope prior scale \eqn{\sigma_j(s)}
#'         and an anchor \eqn{\kappa_j(s)}; the intercept uses a wide prior.
#'   \item Builds a small grid: intercept mean offsets, global multipliers on the \eqn{\sigma_j},
#'         and \eqn{\kappa} around the anchor average; runs RW–MH for each grid point.
#'   \item Selects one run using a score favoring acceptance in a target window, high posterior
#'         predictive agreement, and closeness of standardized slope ratios to **GLM** ratios
#'         (with a single *reference* denominator).
#' }
#'
#' @section Factor handling & column names:
#' \itemize{
#'   \item **Numeric predictors** appear as a single column with their original name
#'         (e.g., `X3`). No suffixes are added.
#'   \item **Factor predictors** are expanded by `model.matrix()` using treatment
#'         contrasts with the *first level as the baseline*. For a two-level factor
#'         `X3` with levels `A` and `B` (baseline = `A`), the encoded
#'         design includes a single dummy column `X3B`, which equals 1 when
#'         `X3 == "B"` and 0 when `X3 == "A"`. The coefficient for `X3B`
#'         is the log-odds difference *B vs A* (controlling for other predictors).
#'   \item The output `posterior$effects$Predictor` uses these encoded names: numeric
#'         predictors as `"Xk"`; factor dummies as `"FactorLevel"` (e.g., `X3B`).
#'   \item **Change the baseline** beforehand to alter dummy labels:
#' \preformatted{X$X3 <- stats::relevel(X$X3, ref = "B")  # baseline becomes B; dummy shows as X3A}
#'   \item If you truly intend to *treat a factor as numeric*, convert it yourself:
#' \preformatted{X$X3 <- as.numeric(X$X3)  # no dummy expansion; column remains 'X3'}
#' }
#'
#' @param y Numeric binary vector (0/1; logical or 2-level factor/character accepted; coerced to 0/1).
#' @param X Matrix or data.frame of predictors (no intercept). May include factors.
#'          Rows must align with `y`.
#' @param missing One of `"complete"` or `"impute"`; the same choice is applied
#'          to both modeling and severity diagnostics (shared rows/values). Default `"complete"`.
#' @param impute_args Optional list controlling simple external imputation when
#'          `missing = "impute"`. Supported:
#'          `numeric_method = "median"|"mean"` (default `"median"`),
#'          `factor_method = "mode"` (default `"mode"`).
#'
#' @param n_iter_grid Integer; MH iterations per grid point (including burn-in). Default `10000`.
#' @param burn_in_grid Integer; burn-in iterations per grid point. Default `1000`.
#' @param step_size Proposal standard deviation for RW–MH. Default `0.40`.
#'
#' @param mu_intercept_offsets Numeric vector of offsets added to \eqn{\mathrm{logit}(\bar{y})}
#'          for the intercept prior mean grid. Default `seq(-1, 1, by = 0.2)`.
#' @param sigma0_intercept Prior \eqn{\sigma_0} (sd) for the intercept (logit scale). Default `10`.
#' @param sigma_global_multipliers Numeric vector of global multipliers applied to all
#'          slope prior scales (after severity anchoring). Default `c(0.1, 0.5, 1, 2, 5, 10)`.
#' @param sigma_hi Slope prior sd under mild separation (\eqn{s=0}). Default `5`.
#' @param sigma_lo Slope prior sd under severe separation (\eqn{s=1}). Default `0.15`.
#' @param kappa_min,kappa_max EP shape at \eqn{s=0} and \eqn{s=1}. Defaults `1`, `2.5`.
#' @param kappa_delta Offsets added around the anchor-average \eqn{\kappa} to form the grid.
#'          Default `seq(-0.5, 0.5, by = 0.2)` (truncated to the interval \eqn{[0.5, 3]}).
#'
#' @param accept_window Numeric length-2 vector; acceptable MH acceptance interval.
#'          Default `c(0.30, 0.40)`.
#' @param accept_target Scalar acceptance target used as a fallback (closest is preferred)
#'          when no grid point falls inside `accept_window`. Default `0.35`.
#'
#' @param ref Predictor *name* or *index* (in the original `X`) to serve as
#'          ratio denominator. If `NULL` (default), the predictor with the highest
#'          univariate severity is used. If the reference is a factor, the denominator
#'          column is the *first* dummy generated for that factor.
#' @param transform_back One of `"none"`, `"logit"`, `"SAS"`, `"Long"`. Controls which
#'          back-transform is reported alongside standardized coefficients. Default `"none"`.
#'          \itemize{
#'            \item `"logit"`: per-**unit** effect on the original encoded column,
#'                  computed as \eqn{\beta^{\mathrm{std}}/s_x}.
#'            \item `"SAS"`: \eqn{\beta^{\mathrm{SAS}} = \beta^{\mathrm{logit}} \cdot \pi/\sqrt{3}}.
#'            \item `"Long"`: \eqn{\beta^{\mathrm{Long}} = \beta^{\mathrm{logit}} \cdot (\pi/\sqrt{3} + 1)}.
#'          }
#' @param ci_level Credible interval level in (0,1). Default `0.95`.
#' @param seed Optional integer; if provided, sets RNG seed for reproducibility.
#' @param return_draws Logical; if `TRUE`, return the post-burn MH chain for the
#'          selected grid run.
#'
#' @return A list with components:
#' \itemize{
#'   \item `ref_predictor`: list with `index` (1-based in original `X`) and `name`.
#'   \item `severity`: data.frame with per-*original* predictor severities used to anchor prior scales.
#'   \item `grid_summary`: data.frame summarizing all grid runs
#'         (`grid_id`, `mu`, `sigma_diag`, `kappa`, `acceptance_rate`,
#'          `prop_matched`, `posterior_ratio_std`).
#'   \item `best`: list describing the selected grid point
#'         (`grid_id`, `mu`, `sigma_diag`, `kappa`,
#'          `acceptance_rate`, `prop_matched`).
#'   \item `posterior`: list with
#'         \itemize{
#'           \item `means_std`: posterior means on the standardized encoded design
#'                 (length = intercept + encoded slopes),
#'           \item `effects`: data.frame with per encoded column:
#'                 `Predictor`, `Scaled`, `Scaled_CI_low`, `Scaled_CI_high`,
#'                 plus **one optional trio** depending on `transform_back`:
#'                 `b_logit_original`/`b_logit_CI_low`/`b_logit_CI_high` *or*
#'                 `b_SAS_original`/`b_SAS_CI_low`/`b_SAS_CI_high` *or*
#'                 `b_Long_original`/`b_Long_CI_low`/`b_Long_CI_high`.
#'         }
#'   \item `draws`: matrix of MH draws after burn-in for the selected run
#'         (returned only when `return_draws = TRUE`).
#' }
#'
#' @details
#' **Standardization & back-transforms.** The sampler runs on z-scored encoded columns.
#' `Scaled` slopes are posterior means on this working scale. Let \eqn{s_x} be the SD of the
#' unscaled encoded column: `logit` back-transform uses \eqn{\beta^{\mathrm{std}}/s_x};
#' `SAS` and `Long` multiply the `logit` transform by \eqn{\pi/\sqrt{3}} and \eqn{\pi/\sqrt{3}+1}
#' respectively.
#'
#' **Reference predictor & ratios.** The ratio denominator is chosen from the original
#' `X` columns (highest DISCO severity by default). If that predictor is a factor,
#' the denominator is the *first* dummy column generated for that factor by
#' `model.matrix()`. GLM ratios are computed on the standardized encoded design.
#'
#' **Shared missing handling.** The `missing` choice governs a single preprocessing
#' step applied to `(y, X)`. With `missing = "complete"` we drop rows with any NA in
#' `y` or `X`. With `missing = "impute"` we drop rows with NA in `y` and impute
#' NAs in `X` using `impute_args`. The same processed data are then used for both the
#' MEP fit and the DISCO severities (we call `DISCO::uni_separation(..., missing = "complete")`
#' because the data are already prepared). **Numeric predictors are z-scored when computing
#' DISCO severities** (factors unchanged) to align the anchoring scale with the modeling scale.
#'
#' @seealso `DISCO::uni_separation`, `DISCO::latent_separation`
#'
#' @importFrom DISCO uni_separation
#' @importFrom stats glm binomial coef plogis qlogis sd model.matrix terms median quantile
#' @examples
#' \donttest{
#' ## Toy data with one numeric factor and one 2-level factor
#' y <- c(0,0,0,0, 1,1,1,1)
#' X_toy <- data.frame(
#'   X1 = c(-1.86, -0.81,  1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
#'   X2 = c( 0.52,  -0.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09),
#'   X3 = factor(c(rep("A", 4), rep("B", 4)))
#' )
#'
#' ## 0) Univariate DISCO diagnostics (complete-case)
#' d3 <- DISCO::uni_separation(data.frame(y=y, X3=X_toy$X3), "X3", "y", "complete")
#' d1 <- DISCO::uni_separation(data.frame(y=y, X1=X_toy$X1), "X1", "y", "complete")
#' d2 <- DISCO::uni_separation(data.frame(y=y, X2=X_toy$X2), "X2", "y", "complete")
#' d3$separation_type; d1$separation_type; d2$separation_type
#'
#' ## 1) Standardized-only coefficients (includes CIs for Scaled slopes)
#' fit_std <- MEP_mixture(
#'   y, X_toy,
#'   n_iter_grid = 4000, burn_in_grid = 1000, seed = 42,
#'   transform_back = "none", ci_level = 0.95
#' )
#' head(fit_std$posterior$effects)
#' # Columns: Predictor, Scaled, Scaled_CI_low, Scaled_CI_high
#'
#' ## 2) Back-transform to LOGIT (per-encoded-unit effect) with CIs
#' fit_logit <- MEP_mixture(
#'   y, X_toy,
#'   n_iter_grid = 4000, burn_in_grid = 1000, seed = 42,
#'   transform_back = "logit"
#' )
#' subset(fit_logit$posterior$effects,
#'        select = c("Predictor","b_logit_original","b_logit_CI_low","b_logit_CI_high"))
#'
#' ## 3) Alternative effect scales with CIs (choose ONE per run)
#' fit_sas <- MEP_mixture(
#'   y, X_toy,
#'   n_iter_grid = 4000, burn_in_grid = 1000, seed = 42,
#'   transform_back = "SAS"
#' )
#' fit_long <- MEP_mixture(
#'   y, X_toy,
#'   n_iter_grid = 4000, burn_in_grid = 1000, seed = 42,
#'   transform_back = "Long"
#' )
#' subset(fit_sas$posterior$effects,
#'        select = c("Predictor","b_SAS_original","b_SAS_CI_low","b_SAS_CI_high"))
#' subset(fit_long$posterior$effects,
#'        select = c("Predictor","b_Long_original","b_Long_CI_low","b_Long_CI_high"))
#'
#' ## 4) Change baseline level to rename the dummy (now "X3A" = A vs B)
#' X_toy2 <- X_toy
#' X_toy2$X3 <- stats::relevel(X_toy2$X3, ref = "B")
#' fit_base <- MEP_mixture(
#'   y, X_toy2,
#'   n_iter_grid = 3000, burn_in_grid = 800, seed = 7,
#'   transform_back = "logit"
#' )
#' head(fit_base$posterior$effects)   # dummy appears as "X3A"
#'
#' ## 5) Shared missing handling
#' X_miss <- X_toy
#' X_miss$X1[c(2,6)] <- NA      # numeric NA
#' X_miss$X3[7]      <- NA      # factor NA
#'
#' # (a) complete-case: drops rows with any NA in X
#' fit_cc <- MEP_mixture(
#'   y, X_miss,
#'   missing = "complete",
#'   n_iter_grid = 3000, burn_in_grid = 800, seed = 9, transform_back = "Long"
#' )
#' head(fit_cc$posterior$effects)
#'
#' # (b) impute: imputes X (numeric=median; factor=mode by default),
#' #     then uses the same imputed data for DISCO + model
#' fit_im <- MEP_mixture(
#'   y, X_miss,
#'   missing = "impute",
#'   impute_args = list(numeric_method = "median"),
#'   n_iter_grid = 3000, burn_in_grid = 800, seed = 9, transform_back = "logit"
#' )
#' head(fit_im$posterior$effects)
#'
#' ## 6) Inspect selection & ratios
#' fit <- fit_logit
#' fit$severity
#' head(fit$grid_summary)
#' fit$best
#' }
#' @export
MEP_mixture <- function(
    y, X,
    missing = c("complete","impute"),
    impute_args = list(),
    n_iter_grid = 10000,
    burn_in_grid = 1000,
    step_size = 0.40,
    mu_intercept_offsets = seq(-1, 1, by = 0.2),
    sigma0_intercept = 10,
    sigma_global_multipliers = c(0.1, 0.5, 1, 2, 5, 10),
    sigma_hi = 5,
    sigma_lo = 0.15,
    kappa_min = 1,
    kappa_max = 2.5,
    kappa_delta = seq(-0.5, 0.5, by = 0.2),
    accept_window = c(0.30, 0.40),
    accept_target = 0.35,
    ref = NULL,
    transform_back = c("none","logit","SAS","Long"),
    ci_level = 0.95,
    seed = NULL,
    return_draws = FALSE
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  impute_frame <- function(df, args = list()) {
    num_method <- tolower(args$numeric_method %||% "median")
    fac_method <- tolower(args$factor_method  %||% "mode")
    for (nm in names(df)) {
      v <- df[[nm]]
      if (is.character(v) || is.logical(v)) v <- as.factor(v)
      if (is.numeric(v)) {
        if (anyNA(v)) {
          fill <- if (identical(num_method, "mean")) mean(v, na.rm = TRUE) else stats::median(v, na.rm = TRUE)
          if (!is.finite(fill)) fill <- 0
          v[is.na(v)] <- fill
        }
      } else if (is.factor(v)) {
        if (anyNA(v)) {
          tab <- table(v, useNA = "no")
          if (length(tab)) {
            lvl <- names(tab)[which.max(tab)]
            v[is.na(v)] <- factor(lvl, levels = levels(v))
          } else {
            v <- factor(v, levels = c(levels(v), "Missing"))
            v[is.na(v)] <- "Missing"
          }
        }
      } else {
        if (anyNA(v)) v[is.na(v)] <- 0
      }
      df[[nm]] <- v
    }
    df
  }

  scale_if_num <- function(v) if (is.numeric(v)) as.numeric(scale(v)) else v

  missing <- match.arg(missing)
  transform_back <- match.arg(transform_back)
  if (!is.null(seed)) set.seed(as.integer(seed))

  # normalize y
  if (!is.numeric(y)) {
    if (is.logical(y)) y <- as.integer(y)
    else if (is.factor(y) || is.character(y)) y <- as.integer(factor(y)) - 1L
    else stop("`y` must be numeric 0/1, logical, or 2-level factor/character.", call. = FALSE)
  }

  # data prep (shared for modeling & DISCO)
  X_df <- as.data.frame(X, stringsAsFactors = TRUE)
  stopifnot(nrow(X_df) == length(y))

  if (missing == "complete") {
    idx <- stats::complete.cases(data.frame(y = y, X_df, check.names = FALSE))
    y <- y[idx]; X_df <- X_df[idx, , drop = FALSE]
  } else { # "impute"
    keep <- !is.na(y)
    y <- y[keep]; X_df <- X_df[keep, , drop = FALSE]
    X_df <- impute_frame(X_df, impute_args)
  }
  if (ncol(X_df) < 1L) stop("X must have at least one predictor.", call. = FALSE)
  if (length(y) < 2L || length(unique(y)) < 2L) stop("Outcome must contain both 0 and 1 after missing handling.", call. = FALSE)

  # encoded design (no intercept), and mapping to original terms
  mm <- stats::model.matrix(~ . , data = X_df)  # includes intercept
  assign_vec <- attr(mm, "assign")
  term_labs  <- attr(stats::terms(~ . , data = X_df), "term.labels")
  X_mm <- mm[, -1, drop = FALSE]
  assign_mm <- assign_vec[-1]
  p_enc <- ncol(X_mm)
  p_terms <- length(term_labs)

  # univariate severities per original predictor (prepared data; numeric predictors scaled)
  sev_vec <- rep(0, p_terms)
  for (j in seq_len(p_terms)) {
    pred <- scale_if_num(X_df[[term_labs[j]]])
    out <- try(DISCO::uni_separation(
      data = data.frame(y = y, x = pred),
      predictor = "x", outcome = "y",
      missing = "complete"
    ), silent = TRUE)
    s <- 0
    if (!inherits(out, "try-error")) {
      if (!is.null(out$severity_score) && is.finite(out$severity_score)) s <- as.numeric(out$severity_score)
    }
    sev_vec[j] <- s
  }
  sev_df <- data.frame(Predictor = term_labs, Severity = sev_vec, stringsAsFactors = FALSE)

  # reference predictor (by highest severity unless specified)
  if (is.null(ref)) {
    ref_idx_term <- which.max(sev_vec)
  } else if (is.character(ref)) {
    ref_idx_term <- match(ref, term_labs); if (is.na(ref_idx_term)) stop("`ref` not found in colnames(X).", call. = FALSE)
  } else {
    ref_idx_term <- as.integer(ref); if (ref_idx_term < 1 || ref_idx_term > p_terms) stop("`ref` index out of range.", call. = FALSE)
  }
  ref_name <- term_labs[ref_idx_term]
  ref_enc_cols <- which(assign_mm == ref_idx_term)
  if (!length(ref_enc_cols)) stop("Internal: no encoded columns for selected `ref` predictor.", call. = FALSE)
  ref_pos_enc <- 1 + ref_enc_cols[1]

  # map severity -> anchor sigma/kappa
  map_uni_severity <- function(s, sigma_hi, sigma_lo, kappa_min, kappa_max) {
    s <- max(0, min(1, as.numeric(s)))
    sigma <- exp((1 - s) * log(sigma_hi) + s * log(sigma_lo))
    kappa_anchor <- kappa_min + s * (kappa_max - kappa_min)
    list(sigma = sigma, kappa_anchor = kappa_anchor)
  }
  anchors <- lapply(sev_vec, map_uni_severity,
                    sigma_hi = sigma_hi, sigma_lo = sigma_lo,
                    kappa_min = kappa_min, kappa_max = kappa_max)
  sigma_anchor_terms <- vapply(anchors, function(a) a$sigma, numeric(1))
  kappa_anchor_mean  <- mean(vapply(anchors, function(a) a$kappa_anchor, numeric(1)))
  sigma_anchor_enc <- sigma_anchor_terms[assign_mm]

  # grids
  mu0 <- stats::qlogis(pmin(pmax(mean(y), 1e-6), 1 - 1e-6))
  mu_grid <- lapply(mu_intercept_offsets, function(off) { v <- numeric(1 + p_enc); v[1] <- mu0 + off; v })
  build_sigma_diag <- function(global_mult = 1.0, sigma0_intercept = 10) {
    d <- c(sigma0_intercept, pmax(1e-6, sigma_anchor_enc * global_mult))
    diag(d, nrow = 1 + p_enc, ncol = 1 + p_enc)
  }
  sigma_grid <- lapply(sigma_global_multipliers, build_sigma_diag, sigma0_intercept = sigma0_intercept)
  kappa_grid <- pmax(0.5, pmin(3.0, kappa_anchor_mean + kappa_delta))

  # MH on standardized encoded design
  run_MH_sampler <- function(n_iter, init_beta, step_size, X_enc, y, mu, Sigma, kappa, burn_in = 1000, transform_back, ci_level) {
    Xstd <- scale(X_enc)
    X <- cbind(Intercept = 1, Xstd)
    Sigma_inv <- solve(Sigma)
    log_post <- function(beta) {
      eta <- as.vector(X %*% beta)
      pr  <- 1 / (1 + exp(-eta))
      loglik <- sum(y * log(pmax(pr, 1e-12)) + (1 - y) * log(pmax(1 - pr, 1e-12)))
      diff <- beta - mu
      qf   <- as.numeric(t(diff) %*% Sigma_inv %*% diff)
      logprior <- -0.5 * (qf^kappa)
      loglik + logprior
    }
    p_all <- ncol(X)
    chain <- matrix(0, n_iter, p_all); chain[1, ] <- rep(0.01, p_all)
    cur_lp <- log_post(chain[1, ]); acc <- 0L
    for (t in 2:n_iter) {
      prop <- chain[t - 1, ] + step_size * rnorm(p_all)
      prop_lp <- log_post(prop)
      if (log(runif(1)) < (prop_lp - cur_lp)) {
        chain[t, ] <- prop; cur_lp <- prop_lp; acc <- acc + 1L
      } else chain[t, ] <- chain[t - 1, ]
      if (t <= burn_in && (t %% 1000) == 0) {
        ar <- acc / t
        if (ar > 0.45) step_size <- step_size * 1.10
        if (ar < 0.20) step_size <- step_size * 0.90
      }
    }
    post <- chain[(burn_in + 1):n_iter, , drop = FALSE]
    pm   <- colMeans(post)

    # CIs
    qlo <- (1 - ci_level) / 2; qhi <- 1 - qlo
    qfun <- function(M) t(apply(M, 2, stats::quantile, probs = c(qlo, qhi), na.rm = TRUE))

    # standardized slopes (drop intercept)
    draws_scaled <- post[, -1, drop = FALSE]
    scaled_mean  <- colMeans(draws_scaled)
    scaled_ci    <- qfun(draws_scaled)

    # back-transforms for encoded columns (draw-by-draw)
    s_x <- apply(X_enc, 2, stats::sd); s_x[!is.finite(s_x) | s_x == 0] <- 1
    draws_logit <- sweep(draws_scaled, 2, s_x, "/")
    logit_sd <- pi / sqrt(3)
    draws_SAS  <- draws_logit * logit_sd
    draws_Long <- draws_logit * (logit_sd + 1)

    # per-encoded-column effects table (means + CIs)
    effects <- data.frame(
      Predictor      = colnames(X_enc),
      Scaled         = scaled_mean,
      Scaled_CI_low  = scaled_ci[, 1],
      Scaled_CI_high = scaled_ci[, 2],
      row.names = NULL, check.names = FALSE
    )
    if (transform_back == "logit") {
      ci_bt <- qfun(draws_logit)
      effects$b_logit_original <- colMeans(draws_logit)
      effects$b_logit_CI_low   <- ci_bt[, 1]
      effects$b_logit_CI_high  <- ci_bt[, 2]
    } else if (transform_back == "SAS") {
      ci_bt <- qfun(draws_SAS)
      effects$b_SAS_original <- colMeans(draws_SAS)
      effects$b_SAS_CI_low   <- ci_bt[, 1]
      effects$b_SAS_CI_high  <- ci_bt[, 2]
    } else if (transform_back == "Long") {
      ci_bt <- qfun(draws_Long)
      effects$b_Long_original <- colMeans(draws_Long)
      effects$b_Long_CI_low   <- ci_bt[, 1]
      effects$b_Long_CI_high  <- ci_bt[, 2]
    }

    # mean-based vectors for ratio reporting
    b_logit_mean <- colMeans(draws_logit)
    b_SAS_mean   <- colMeans(draws_SAS)
    b_Long_mean  <- colMeans(draws_Long)

    # posterior predictive check
    n_rep <- nrow(post); n_obs <- nrow(X)
    y_rep <- matrix(0, nrow = n_rep, ncol = n_obs)
    for (i in seq_len(n_rep)) {
      pr_i <- 1 / (1 + exp(-(X %*% post[i, ])))
      y_rep[i, ] <- rbinom(n_obs, 1, pr_i)
    }
    p_match <- colMeans(sweep(y_rep, 2, y, `==`))
    prop_matched <- mean(p_match >= 0.80, na.rm = TRUE)
    if (!is.finite(prop_matched)) prop_matched <- 0

    list(
      chain = post,
      pm = pm,
      effects = effects,
      bt = list(logit = b_logit_mean, SAS = b_SAS_mean, Long = b_Long_mean),
      prop_matched = prop_matched,
      acceptance_rate = acc / (n_iter - 1)
    )
  }

  # grid loop
  grid_summary <- data.frame(
    grid_id = integer(),
    mu = character(),
    sigma_diag = character(),
    kappa = numeric(),
    acceptance_rate = numeric(),
    prop_matched = numeric(),
    posterior_ratio_std = character(),
    stringsAsFactors = FALSE
  )
  runs <- list(); gid <- 1L

  for (mu in mu_grid) {
    mu_str <- paste(round(mu, 3), collapse = ", ")
    for (Sigma in sigma_grid) {
      sigma_str <- paste(round(diag(Sigma), 3), collapse = ", ")
      for (kappa in kappa_grid) {
        res <- run_MH_sampler(
          n_iter = n_iter_grid,
          init_beta = rep(0.01, 1 + p_enc),
          step_size = step_size,
          X_enc = X_mm,
          y = y,
          mu = mu,
          Sigma = Sigma,
          kappa = kappa,
          burn_in = burn_in_grid,
          transform_back = transform_back,
          ci_level = ci_level
        )
        vals <- res$pm
        num_idx <- setdiff(2:(1 + p_enc), ref_pos_enc)
        ratio_std_str <- if (length(num_idx) && is.finite(vals[ref_pos_enc]) && abs(vals[ref_pos_enc]) > 1e-12) {
          paste(round(vals[num_idx] / vals[ref_pos_enc], 3), collapse = ", ")
        } else NA_character_

        grid_summary <- rbind(
          grid_summary,
          data.frame(
            grid_id = gid,
            mu = mu_str,
            sigma_diag = sigma_str,
            kappa = kappa,
            acceptance_rate = res$acceptance_rate,
            prop_matched = res$prop_matched,
            posterior_ratio_std = ratio_std_str,
            stringsAsFactors = FALSE
          )
        )
        runs[[gid]] <- res
        gid <- gid + 1L
      }
    }
  }

  # GLM (on standardized encoded design)
  X_scaled <- scale(X_mm)
  glm_fit <- try(suppressWarnings(stats::glm(y ~ ., data = data.frame(y = y, X_scaled), family = stats::binomial())), silent = TRUE)
  glm_ok <- !(inherits(glm_fit, "try-error"))
  glm_coef <- if (glm_ok) stats::coef(glm_fit) else rep(NA_real_, 1 + p_enc)

  safe_ratio_from_glm <- function(coefs, ref_pos) {
    if (length(coefs) >= ref_pos && is.finite(coefs[ref_pos]) && abs(coefs[ref_pos]) > 0) {
      num <- coefs[-c(1, ref_pos)]
      if (length(num) == 0) return(NA_character_)
      paste(round(num / coefs[ref_pos], 3), collapse = ", ")
    } else NA_character_
  }
  parse_ratio <- function(s) {
    if (is.null(s) || is.na(s) || s == "") return(NA_real_)
    as.numeric(strsplit(s, ",\\s*")[[1]])
  }
  ref_ratio_str <- safe_ratio_from_glm(glm_coef, ref_pos_enc)
  ref_ratio_vec <- parse_ratio(ref_ratio_str)

  # select best
  z <- grid_summary
  z2 <- subset(z, is.finite(acceptance_rate) & acceptance_rate >= accept_window[1] & acceptance_rate <= accept_window[2])
  if (nrow(z2) == 0) {
    z2 <- z
    z2$acc_term <- abs(z2$acceptance_rate - accept_target)
  } else {
    z2$acc_term <- 0
  }
  z2$prop_term <- -z2$prop_matched
  z2$ratio_term <- NA_real_
  if (!all(is.na(ref_ratio_vec))) {
    z2$ratio_term <- vapply(z2$posterior_ratio_std, function(rs) {
      pv <- parse_ratio(rs)
      if (all(is.na(pv))) return(NA_real_)
      if (length(pv) != length(ref_ratio_vec)) return(NA_real_)
      mean(abs(pv - ref_ratio_vec))
    }, numeric(1))
  }
  if (any(is.finite(z2$ratio_term))) {
    z2$score <- z2$acc_term + z2$prop_term + z2$ratio_term
    z2 <- z2[order(z2$score, z2$acc_term, z2$prop_term, na.last = TRUE), , drop = FALSE]
  } else {
    z2$score <- z2$acc_term + z2$prop_term
    z2 <- z2[order(z2$score, z2$acc_term, na.last = TRUE), , drop = FALSE]
  }
  best_id <- z2$grid_id[1]
  best_run <- runs[[best_id]]
  best_row <- z[z$grid_id == best_id, , drop = FALSE]

  # outputs
  effects <- best_run$effects
  pm      <- best_run$pm
  #safe_div <- function(num, den) ifelse(is.finite(den) & abs(den) > 1e-12, num / den, NA_real_)

  glm_beta <- glm_coef[-1]
  #glm_ratio <- safe_div(glm_beta, glm_beta[ref_enc_cols[1]])

  mep_beta_std <- pm[-1]
  #mep_ratio_std <- safe_div(mep_beta_std, mep_beta_std[ref_enc_cols[1]])

  #ratios_list <- list(
   # GLM_beta = glm_beta,
   # GLM_ratio = glm_ratio,
   # MEP_beta_std = mep_beta_std,
  #  MEP_ratio_std = mep_ratio_std
  #)

  if (transform_back != "none") {
    bt_vec <- best_run$bt[[transform_back]]
    names(bt_vec) <- colnames(X_mm)
    beta_name <- paste0("MEP_beta_b_", transform_back, "_original")
    ratio_name <- paste0("MEP_ratio_b_", transform_back, "_original")
    denom_name <- colnames(X_mm)[ref_enc_cols[1]]

    #ratios_list[[beta_name]]  <- bt_vec
    #ratios_list[[ratio_name]] <- safe_div(bt_vec, bt_vec[denom_name])
  }

  list(
    ref_predictor = list(index = ref_idx_term, name = ref_name),
    severity = sev_df,
    grid_summary = grid_summary,
    best = list(
      grid_id = best_id,
      mu = best_row$mu,
      sigma_diag = best_row$sigma_diag,
      kappa = best_row$kappa,
      acceptance_rate = best_row$acceptance_rate,
      prop_matched = best_row$prop_matched
    ),
    posterior = list(
      means_std = pm,
      effects = effects
    ),
    #ratios = ratios_list,
    draws = if (isTRUE(return_draws)) best_run$chain else NULL
  )
}
