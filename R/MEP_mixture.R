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
#' @param n_iter Integer; MH iterations per grid point (including burn-in). Default `10000`.
#' @param burn_in Integer; burn-in iterations per grid point. Default `1000`.
#' @param init_beta Initial value(s) for the MH chain. Default `0.01`.
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
#' @param ess_threshold Numeric; minimum effective sample size (ESS) required (across all
#'        parameters) to declare convergence when `coda` is available. Default `150`.
#' @param geweke_z_threshold Numeric; maximum allowed absolute Geweke z-score (across all
#'        parameters) to declare convergence when `coda` is available. Default `2`.
#' @param n_chains_best Integer; number of independent MH chains to rerun for the selected
#'        best grid point only. Default `1`. When `n_chains_best >= 2` and package `coda`
#'        is available, multi-chain diagnostics (Gelman-Rubin R-hat and ESS) are computed.
#' @param chain_seeds_best Optional integer vector of length `n_chains_best`. Each value is used
#'        as the RNG seed for the corresponding best-point chain. If `NULL`, seeds are generated
#'        deterministically from `seed` (if provided) or from a random base seed.
#' @param combine_chains One of `"stack"` or `"none"`. Controls how the selected best-point
#'        chains are summarized:
#'        \itemize{
#'          \item `"stack"`: stack post burn-in draws across all best-point chains, then compute
#'                posterior means, credible intervals, and convergence diagnostics from the pooled draws.
#'          \item `"none"`: summarize using only the first best-point chain.
#'        }
#'        Default `"stack"`.
#'
#' @return A list with components:
#' \itemize{
#'   \item `ref_predictor`: list with `index` (1-based in original `X`) and `name`.
#'   \item `severity`: data.frame with per-*original* predictor severities used to anchor prior scales:
#'         columns `Predictor`, `Severity`.
#'
#'   \item `grid_summary`: data.frame summarizing all grid runs with columns:
#'         \itemize{
#'           \item `grid_id`: integer id for the grid point.
#'           \item `mu`: character representation of the prior mean vector (intercept entry varies; others are 0).
#'           \item `sigma_diag`: character representation of the diagonal of the prior scale matrix.
#'           \item `kappa`: EP shape parameter used for that run.
#'           \item `acceptance_rate`: MH acceptance rate.
#'           \item `prop_matched`: posterior predictive agreement summary (proportion of observations with
#'                 per-observation match probability >= 0.80, averaged).
#'           \item `posterior_ratio_std`: character representation of standardized posterior mean ratios
#'                 relative to the reference predictor (when defined).
#'           \item `converged`: logical; convergence flag for this run (requires package `coda`, otherwise `NA`).
#'           \item `ess_min`: numeric; minimum ESS across parameters for this run (requires `coda`, otherwise `NA`).
#'           \item `geweke_max_abs`: numeric; maximum absolute Geweke z-score across parameters for this run
#'                 (requires `coda`, otherwise `NA`).
#'         }
#'
#'   \item `best`: list describing the selected grid point with elements:
#'         `grid_id`, `mu`, `sigma_diag`, `kappa`, `acceptance_rate`, `prop_matched`,
#'         and diagnostics copied from the selected run:
#'         `converged`, `ess_min`, `geweke_max_abs`.
#'
#'   \item `posterior`: list with
#'         \itemize{
#'           \item `means_std`: posterior means on the working scale (intercept + standardized encoded slopes),
#'                 length = `1 + p_enc`.
#'           \item `effects`: data.frame with one row per encoded predictor column (slopes only) containing:
#'                 \itemize{
#'                   \item `Predictor`: encoded column name (from `model.matrix`, intercept excluded).
#'                   \item `Scaled`: posterior mean on standardized encoded design.
#'                   \item `Scaled_CI_low`, `Scaled_CI_high`: credible interval endpoints on the standardized scale.
#'                   \item `sig_scaled`: logical; `TRUE` if the `Scaled` credible interval excludes 0.
#'                   \item `star_scaled`: `"*"` if `sig_scaled` is `TRUE`, else `""`.
#'                   \item If `transform_back = "logit"`: `b_logit_original`, `b_logit_CI_low`, `b_logit_CI_high`,
#'                         plus `sig_original`, `star_original` computed from that interval.
#'                   \item If `transform_back = "SAS"`:  `b_SAS_original`, `b_SAS_CI_low`, `b_SAS_CI_high`,
#'                         plus `sig_original`, `star_original`.
#'                   \item If `transform_back = "Long"`: `b_Long_original`, `b_Long_CI_low`, `b_Long_CI_high`,
#'                         plus `sig_original`, `star_original`.
#'                 }
#'         }
#'
#'   \item `convergence`: list of diagnostics for the selected run (requires `coda`; otherwise entries are `NA`):
#'         \itemize{
#'           \item `ess`: numeric vector of ESS for each parameter (intercept + slopes).
#'           \item `geweke_z`: numeric vector of Geweke z-scores for each parameter.
#'           \item `ess_min`: minimum ESS across parameters.
#'           \item `geweke_max_abs`: maximum absolute Geweke z-score across parameters.
#'           \item `converged`: logical; `TRUE` if `ess_min >= ess_threshold` and
#'                 `geweke_max_abs <= geweke_z_threshold`.
#'         }
#'   \item `diagnostics_multi`: list of multi-chain diagnostics for the best grid point.
#'         Computed only when `n_chains_best >= 2` and package `coda` is available; otherwise NA.
#'         \itemize{
#'           \item `rhat`: numeric vector of Gelman-Rubin point estimates (one per parameter).
#'           \item `rhat_max`: maximum R-hat across parameters.
#'           \item `ess`: numeric vector of effective sample sizes from the multi-chain object.
#'           \item `ess_min`: minimum ESS across parameters.
#'         }
#'
#'   \item `best_chains`: summary of per-chain diagnostics and performance for the best grid point:
#'         \itemize{
#'           \item `diagnostics`: list of per-chain single-chain diagnostics (ESS, Geweke, etc).
#'           \item `acceptance_rate`: numeric vector, one per chain.
#'           \item `prop_matched`: numeric vector, one per chain.
#'           \item `chain_seeds`: integer vector of seeds used.
#'         }
#'
#'   \item `best$diagnostics_multi`: a copy of `diagnostics_multi` stored inside `best` for convenience.
#'
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
#' @importFrom stats glm binomial coef plogis qlogis sd model.matrix terms median quantile rbinom complete.cases
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
#'   n_iter = 4000, burn_in = 1000, seed = 42,
#'   transform_back = "none", ci_level = 0.95
#' )
#' head(fit_std$posterior$effects)
#' # Columns: Predictor, Scaled, Scaled_CI_low, Scaled_CI_high
#'
#' ## 2) Back-transform to LOGIT (per-encoded-unit effect) with CIs
#' fit_logit <- MEP_mixture(
#'   y, X_toy,
#'   n_iter = 4000, burn_in = 1000, seed = 42,
#'   transform_back = "logit"
#' )
#' subset(fit_logit$posterior$effects,
#'        select = c("Predictor","b_logit_original","b_logit_CI_low","b_logit_CI_high"))
#'
#' ## 3) Alternative effect scales with CIs (choose ONE per run)
#' fit_sas <- MEP_mixture(
#'   y, X_toy,
#'   n_iter = 4000, burn_in = 1000, seed = 42,
#'   transform_back = "SAS"
#' )
#' fit_long <- MEP_mixture(
#'   y, X_toy,
#'   n_iter = 4000, burn_in = 1000, seed = 42,
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
#'   n_iter = 3000, burn_in = 800, seed = 7,
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
#'   n_iter = 3000, burn_in = 800, seed = 9, transform_back = "Long"
#' )
#' head(fit_cc$posterior$effects)
#'
#' # (b) impute: imputes X (numeric=median; factor=mode by default),
#' #     then uses the same imputed data for DISCO + model
#' fit_im <- MEP_mixture(
#'   y, X_miss,
#'   missing = "impute",
#'   impute_args = list(numeric_method = "median"),
#'   n_iter = 3000, burn_in = 800, seed = 9, transform_back = "logit"
#' )
#' head(fit_im$posterior$effects)
#'
#' ## 6) Inspect selection & ratios
#' fit <- fit_logit
#' fit$severity
#' head(fit$grid_summary)
#' fit$best
#' }
#'
#' ## 7) Rerun the selected best grid point with multiple chains
#' fit_multi <- MEP_mixture(
#'   y, X_toy,
#'   n_iter = 4000, burn_in = 1000,
#'   n_chains_best = 4,
#'   chain_seeds_best = c(101, 102, 103, 104),
#'   combine_chains = "stack",
#'   transform_back = "none",
#'   seed = 9
#' )
#' fit_multi$diagnostics_multi
#' head(fit_multi$posterior$effects)
#' @export
MEP_mixture <- function(
    y, X,
    missing = c("complete","impute"),
    impute_args = list(),
    n_iter = 10000,
    burn_in = 1000,
    init_beta = 0.01,
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
    return_draws = FALSE,
    ess_threshold = 150,
    geweke_z_threshold = 2,
    n_chains_best = 1,
    chain_seeds_best = NULL,
    combine_chains = c("stack","none")
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

  missing <- base::match.arg(as.character(missing), choices = c("complete","impute"))
  transform_back <- base::match.arg(as.character(transform_back), choices = c("none","logit","SAS","Long"))
  combine_chains <- base::match.arg(as.character(combine_chains), choices = c("stack","none"))

  if (!is.null(seed)) set.seed(as.integer(seed))

  if (!is.numeric(n_chains_best) || length(n_chains_best) != 1 || n_chains_best < 1) {
    stop("`n_chains_best` must be a positive integer.", call. = FALSE)
  }
  n_chains_best <- as.integer(n_chains_best)

  if (!is.null(chain_seeds_best)) {
    if (length(chain_seeds_best) != n_chains_best) {
      stop("`chain_seeds_best` must have length `n_chains_best`.", call. = FALSE)
    }
    chain_seeds_best <- as.integer(chain_seeds_best)
  } else {
    base_seed <- if (!is.null(seed)) as.integer(seed) else sample.int(1e9, 1)
    chain_seeds_best <- base_seed + seq_len(n_chains_best)
  }

  # normalize y
  if (!is.numeric(y)) {
    if (is.logical(y)) y <- as.integer(y)
    else if (is.factor(y) || is.character(y)) y <- as.integer(factor(y)) - 1L
    else stop("`y` must be numeric 0/1, logical, or 2-level factor/character.", call. = FALSE)
  }
  y <- as.numeric(y)

  # data prep
  X_df <- as.data.frame(X, stringsAsFactors = TRUE)
  stopifnot(nrow(X_df) == length(y))

  if (missing == "complete") {
    idx <- stats::complete.cases(data.frame(y = y, X_df, check.names = FALSE))
    y <- y[idx]
    X_df <- X_df[idx, , drop = FALSE]
  } else {
    keep <- !is.na(y)
    y <- y[keep]
    X_df <- X_df[keep, , drop = FALSE]
    X_df <- impute_frame(X_df, impute_args)
  }

  if (ncol(X_df) < 1L) stop("X must have at least one predictor.", call. = FALSE)
  if (length(y) < 2L || length(unique(y)) < 2L) stop("Outcome must contain both 0 and 1 after missing handling.", call. = FALSE)

  # encoded design
  mm <- stats::model.matrix(~ . , data = X_df)
  assign_vec <- attr(mm, "assign")
  term_labs  <- attr(stats::terms(~ . , data = X_df), "term.labels")
  X_mm <- mm[, -1, drop = FALSE]
  assign_mm <- assign_vec[-1]
  p_enc <- ncol(X_mm)
  p_terms <- length(term_labs)

  # severities
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

  # reference predictor
  if (is.null(ref)) {
    ref_idx_term <- which.max(sev_vec)
  } else if (is.character(ref)) {
    ref_idx_term <- match(ref, term_labs)
    if (is.na(ref_idx_term)) stop("`ref` not found in colnames(X).", call. = FALSE)
  } else {
    ref_idx_term <- as.integer(ref)
    if (ref_idx_term < 1 || ref_idx_term > p_terms) stop("`ref` index out of range.", call. = FALSE)
  }

  ref_name <- term_labs[ref_idx_term]
  ref_enc_cols <- which(assign_mm == ref_idx_term)
  if (!length(ref_enc_cols)) stop("Internal: no encoded columns for selected `ref` predictor.", call. = FALSE)
  ref_pos_enc <- 1 + ref_enc_cols[1]

  map_uni_severity <- function(s, sigma_hi, sigma_lo, kappa_min, kappa_max) {
    s <- max(0, min(1, as.numeric(s)))
    sigma <- exp((1 - s) * log(sigma_hi) + s * log(sigma_lo))
    kappa_anchor <- kappa_min + s * (kappa_max - kappa_min)
    list(sigma = sigma, kappa_anchor = kappa_anchor)
  }

  anchors <- lapply(
    sev_vec, map_uni_severity,
    sigma_hi = sigma_hi, sigma_lo = sigma_lo,
    kappa_min = kappa_min, kappa_max = kappa_max
  )
  sigma_anchor_terms <- vapply(anchors, function(a) a$sigma, numeric(1))
  kappa_anchor_mean  <- mean(vapply(anchors, function(a) a$kappa_anchor, numeric(1)))
  sigma_anchor_enc <- sigma_anchor_terms[assign_mm]

  # grids
  mu0 <- stats::qlogis(pmin(pmax(mean(y), 1e-6), 1 - 1e-6))
  mu_grid <- lapply(mu_intercept_offsets, function(off) {
    v <- numeric(1 + p_enc)
    v[1] <- mu0 + off
    v
  })

  build_sigma_diag <- function(global_mult = 1.0, sigma0_intercept = 10) {
    d <- c(sigma0_intercept, pmax(1e-6, sigma_anchor_enc * global_mult))
    diag(d, nrow = 1 + p_enc, ncol = 1 + p_enc)
  }
  sigma_grid <- lapply(sigma_global_multipliers, build_sigma_diag, sigma0_intercept = sigma0_intercept)
  kappa_grid <- pmax(0.5, pmin(3.0, kappa_anchor_mean + kappa_delta))

  summarize_post <- function(post, X_enc, transform_back, ci_level, ess_threshold, geweke_z_threshold) {
    qlo <- (1 - ci_level) / 2
    qhi <- 1 - qlo
    qfun <- function(M) t(apply(M, 2, stats::quantile, probs = c(qlo, qhi), na.rm = TRUE))

    draws_scaled <- post[, -1, drop = FALSE]
    scaled_mean <- colMeans(draws_scaled)
    scaled_ci <- qfun(draws_scaled)

    s_x <- apply(X_enc, 2, stats::sd)
    s_x[!is.finite(s_x) | s_x == 0] <- 1
    draws_logit <- sweep(draws_scaled, 2, s_x, "/")

    logit_sd <- pi / sqrt(3)
    draws_SAS  <- draws_logit * logit_sd
    draws_Long <- draws_logit * (logit_sd + 1)

    effects <- data.frame(
      Predictor      = colnames(X_enc),
      Scaled         = scaled_mean,
      Scaled_CI_low  = scaled_ci[, 1],
      Scaled_CI_high = scaled_ci[, 2],
      row.names = NULL,
      check.names = FALSE
    )

    effects$sig_scaled <- with(
      effects,
      (Scaled_CI_low > 0 & Scaled_CI_high > 0) | (Scaled_CI_low < 0 & Scaled_CI_high < 0)
    )
    effects$star_scaled <- ifelse(effects$sig_scaled, "*", "")

    if (transform_back == "logit") {
      ci_bt <- qfun(draws_logit)
      effects$b_logit_original <- colMeans(draws_logit)
      effects$b_logit_CI_low   <- ci_bt[, 1]
      effects$b_logit_CI_high  <- ci_bt[, 2]
      effects$sig_original <- (effects$b_logit_CI_low > 0 & effects$b_logit_CI_high > 0) |
        (effects$b_logit_CI_low < 0 & effects$b_logit_CI_high < 0)
      effects$star_original <- ifelse(effects$sig_original, "*", "")
    } else if (transform_back == "SAS") {
      ci_bt <- qfun(draws_SAS)
      effects$b_SAS_original <- colMeans(draws_SAS)
      effects$b_SAS_CI_low   <- ci_bt[, 1]
      effects$b_SAS_CI_high  <- ci_bt[, 2]
      effects$sig_original <- (effects$b_SAS_CI_low > 0 & effects$b_SAS_CI_high > 0) |
        (effects$b_SAS_CI_low < 0 & effects$b_SAS_CI_high < 0)
      effects$star_original <- ifelse(effects$sig_original, "*", "")
    } else if (transform_back == "Long") {
      ci_bt <- qfun(draws_Long)
      effects$b_Long_original <- colMeans(draws_Long)
      effects$b_Long_CI_low   <- ci_bt[, 1]
      effects$b_Long_CI_high  <- ci_bt[, 2]
      effects$sig_original <- (effects$b_Long_CI_low > 0 & effects$b_Long_CI_high > 0) |
        (effects$b_Long_CI_low < 0 & effects$b_Long_CI_high < 0)
      effects$star_original <- ifelse(effects$sig_original, "*", "")
    }

    ess <- rep(NA_real_, ncol(post))
    geweke_z <- rep(NA_real_, ncol(post))
    converged_flag <- NA
    ess_min <- NA_real_
    geweke_max_abs <- NA_real_

    if (requireNamespace("coda", quietly = TRUE)) {
      mcmc_obj <- coda::mcmc(post)
      ess <- as.numeric(coda::effectiveSize(mcmc_obj))
      gz  <- coda::geweke.diag(mcmc_obj)$z
      geweke_z <- as.numeric(gz)

      ess_min <- suppressWarnings(min(ess, na.rm = TRUE))
      geweke_max_abs <- suppressWarnings(max(abs(geweke_z), na.rm = TRUE))

      ess_ok <- is.finite(ess_min) && ess_min >= ess_threshold
      geweke_ok <- is.finite(geweke_max_abs) && geweke_max_abs <= geweke_z_threshold
      converged_flag <- isTRUE(ess_ok && geweke_ok)
    }

    list(
      pm = colMeans(post),
      effects = effects,
      diagnostics = list(
        ess = ess,
        geweke_z = geweke_z,
        ess_min = ess_min,
        geweke_max_abs = geweke_max_abs,
        converged = converged_flag
      )
    )
  }

  run_MH_sampler <- function(
    n_iter, init_beta, step_size, X_enc, y, mu, Sigma, kappa,
    burn_in = 1000, transform_back, ci_level, ess_threshold, geweke_z_threshold,
    chain_seed = NULL
  ) {
    if (!is.null(chain_seed)) set.seed(as.integer(chain_seed))

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
    chain <- matrix(0, n_iter, p_all)

    if (length(init_beta) == 1L) {
      init_vec <- rep(as.numeric(init_beta), p_all)
    } else if (length(init_beta) == p_all) {
      init_vec <- as.numeric(init_beta)
    } else {
      stop("`init_beta` must be a scalar or a numeric vector of length (1 + p_enc).", call. = FALSE)
    }

    chain[1, ] <- init_vec
    cur_lp <- log_post(chain[1, ])
    acc <- 0L
    cur_step <- step_size

    for (t in 2:n_iter) {
      prop <- chain[t - 1, ] + cur_step * rnorm(p_all)
      prop_lp <- log_post(prop)

      if (log(runif(1)) < (prop_lp - cur_lp)) {
        chain[t, ] <- prop
        cur_lp <- prop_lp
        acc <- acc + 1L
      } else {
        chain[t, ] <- chain[t - 1, ]
      }

      if (t <= burn_in && (t %% 1000) == 0) {
        ar <- acc / t
        if (ar > 0.45) cur_step <- cur_step * 1.10
        if (ar < 0.20) cur_step <- cur_step * 0.90
      }
    }

    post <- chain[(burn_in + 1):n_iter, , drop = FALSE]

    # posterior predictive check
    n_rep <- nrow(post)
    n_obs <- nrow(X)
    y_rep <- matrix(0, nrow = n_rep, ncol = n_obs)
    for (i in seq_len(n_rep)) {
      pr_i <- 1 / (1 + exp(-(X %*% post[i, ])))
      y_rep[i, ] <- rbinom(n_obs, 1, pr_i)
    }
    p_match <- colMeans(sweep(y_rep, 2, y, `==`))
    prop_matched <- mean(p_match >= 0.80, na.rm = TRUE)
    if (!is.finite(prop_matched)) prop_matched <- 0

    sum_obj <- summarize_post(
      post = post,
      X_enc = X_enc,
      transform_back = transform_back,
      ci_level = ci_level,
      ess_threshold = ess_threshold,
      geweke_z_threshold = geweke_z_threshold
    )

    list(
      chain = post,
      pm = sum_obj$pm,
      effects = sum_obj$effects,
      prop_matched = prop_matched,
      acceptance_rate = acc / (n_iter - 1),
      diagnostics = sum_obj$diagnostics
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
    converged = logical(),
    ess_min = numeric(),
    geweke_max_abs = numeric(),
    stringsAsFactors = FALSE
  )

  runs <- list()
  gid <- 1L

  for (mu in mu_grid) {
    mu_str <- paste(round(mu, 3), collapse = ", ")
    for (Sigma in sigma_grid) {
      sigma_str <- paste(round(diag(Sigma), 3), collapse = ", ")
      for (kappa in kappa_grid) {

        grid_chain_seed <- if (!is.null(seed)) as.integer(seed) + gid else NULL

        res <- run_MH_sampler(
          n_iter = n_iter,
          init_beta = init_beta,
          step_size = step_size,
          X_enc = X_mm,
          y = y,
          mu = mu,
          Sigma = Sigma,
          kappa = kappa,
          burn_in = burn_in,
          transform_back = transform_back,
          ci_level = ci_level,
          ess_threshold = ess_threshold,
          geweke_z_threshold = geweke_z_threshold,
          chain_seed = grid_chain_seed
        )

        vals <- res$pm
        num_idx <- setdiff(2:(1 + p_enc), ref_pos_enc)

        ratio_std_str <- if (length(num_idx) &&
                             is.finite(vals[ref_pos_enc]) && abs(vals[ref_pos_enc]) > 1e-12) {
          paste(round(vals[num_idx] / vals[ref_pos_enc], 3), collapse = ", ")
        } else {
          NA_character_
        }

        conv_val <- if (!is.null(res$diagnostics$converged)) res$diagnostics$converged else NA
        ess_min_val <- if (!is.null(res$diagnostics$ess_min)) res$diagnostics$ess_min else NA_real_
        geweke_max_val <- if (!is.null(res$diagnostics$geweke_max_abs)) res$diagnostics$geweke_max_abs else NA_real_

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
            converged = conv_val,
            ess_min = ess_min_val,
            geweke_max_abs = geweke_max_val,
            stringsAsFactors = FALSE
          )
        )

        runs[[gid]] <- res
        gid <- gid + 1L
      }
    }
  }

  # GLM ratios
  X_scaled <- scale(X_mm)
  glm_fit <- try(suppressWarnings(stats::glm(
    y ~ .,
    data = data.frame(y = y, X_scaled),
    family = stats::binomial()
  )), silent = TRUE)

  glm_ok <- !(inherits(glm_fit, "try-error"))
  glm_coef <- if (glm_ok) stats::coef(glm_fit) else rep(NA_real_, 1 + p_enc)

  safe_ratio_from_glm <- function(coefs, ref_pos) {
    if (length(coefs) >= ref_pos && is.finite(coefs[ref_pos]) && abs(coefs[ref_pos]) > 0) {
      num <- coefs[-c(1, ref_pos)]
      if (length(num) == 0) return(NA_character_)
      paste(round(num / coefs[ref_pos], 3), collapse = ", ")
    } else {
      NA_character_
    }
  }

  parse_ratio <- function(s) {
    if (is.null(s) || is.na(s) || s == "") return(NA_real_)
    as.numeric(strsplit(s, ",\\s*")[[1]])
  }

  ref_ratio_str <- safe_ratio_from_glm(glm_coef, ref_pos_enc)
  ref_ratio_vec <- parse_ratio(ref_ratio_str)

  # select best
  z <- grid_summary
  z2 <- subset(z, is.finite(acceptance_rate) &
                 acceptance_rate >= accept_window[1] &
                 acceptance_rate <= accept_window[2])

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
  best_row <- z[z$grid_id == best_id, , drop = FALSE]

  # decode best grid indices
  n_mu <- length(mu_grid)
  n_sigma <- length(sigma_grid)
  n_kappa <- length(kappa_grid)

  id0 <- best_id - 1L
  i_mu <- (id0 %/% (n_sigma * n_kappa)) + 1L
  rem <- id0 %% (n_sigma * n_kappa)
  i_sigma <- (rem %/% n_kappa) + 1L
  i_kappa <- (rem %% n_kappa) + 1L

  mu_best <- mu_grid[[i_mu]]
  Sigma_best <- sigma_grid[[i_sigma]]
  kappa_best <- kappa_grid[[i_kappa]]

  # rerun best with multiple chains
  best_chains <- vector("list", n_chains_best)
  for (ch in seq_len(n_chains_best)) {
    best_chains[[ch]] <- run_MH_sampler(
      n_iter = n_iter,
      init_beta = init_beta,
      step_size = step_size,
      X_enc = X_mm,
      y = y,
      mu = mu_best,
      Sigma = Sigma_best,
      kappa = kappa_best,
      burn_in = burn_in,
      transform_back = transform_back,
      ci_level = ci_level,
      ess_threshold = ess_threshold,
      geweke_z_threshold = geweke_z_threshold,
      chain_seed = chain_seeds_best[ch]
    )
  }

  prop_matched_best <- mean(vapply(best_chains, function(x) x$prop_matched, numeric(1)), na.rm = TRUE)
  acc_best <- mean(vapply(best_chains, function(x) x$acceptance_rate, numeric(1)), na.rm = TRUE)

  # multi-chain diagnostics (Rhat + multi-chain ESS)
  diagnostics_multi <- list(
    rhat = rep(NA_real_, ncol(best_chains[[1]]$chain)),
    rhat_max = NA_real_,
    ess = rep(NA_real_, ncol(best_chains[[1]]$chain)),
    ess_min = NA_real_
  )

  if (requireNamespace("coda", quietly = TRUE) && n_chains_best >= 2) {
    mlist <- coda::mcmc.list(lapply(best_chains, function(x) coda::mcmc(x$chain)))

    gd <- coda::gelman.diag(mlist, autoburnin = FALSE, multivariate = FALSE)$psrf
    diagnostics_multi$rhat <- as.numeric(gd[, "Point est."])
    diagnostics_multi$rhat_max <- suppressWarnings(max(diagnostics_multi$rhat, na.rm = TRUE))

    ess_m <- coda::effectiveSize(mlist)
    diagnostics_multi$ess <- as.numeric(ess_m)
    diagnostics_multi$ess_min <- suppressWarnings(min(diagnostics_multi$ess, na.rm = TRUE))
  }

  # combine for final posterior summary
  if (combine_chains == "stack") {
    post_all <- do.call(rbind, lapply(best_chains, `[[`, "chain"))
    sum_all <- summarize_post(
      post = post_all,
      X_enc = X_mm,
      transform_back = transform_back,
      ci_level = ci_level,
      ess_threshold = ess_threshold,
      geweke_z_threshold = geweke_z_threshold
    )
    pm <- sum_all$pm
    effects <- sum_all$effects
    best_diag <- sum_all$diagnostics
    convergence_out <- best_diag
  } else {
    pm <- best_chains[[1]]$pm
    effects <- best_chains[[1]]$effects
    best_diag <- best_chains[[1]]$diagnostics
    convergence_out <- best_diag
  }

  draws_out <- NULL
  if (isTRUE(return_draws)) {
    if (n_chains_best == 1L) draws_out <- best_chains[[1]]$chain
    else draws_out <- lapply(best_chains, `[[`, "chain")
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
      acceptance_rate = acc_best,
      prop_matched = prop_matched_best,
      converged = best_diag$converged,
      ess_min = best_diag$ess_min,
      geweke_max_abs = best_diag$geweke_max_abs,
      n_chains_best = n_chains_best,
      chain_seeds_best = chain_seeds_best,
      combine_chains = combine_chains,
      diagnostics_multi = diagnostics_multi
    ),
    posterior = list(
      means_std = pm,
      effects = effects
    ),
    convergence = convergence_out,
    diagnostics_multi = diagnostics_multi,
    best_chains = list(
      diagnostics = lapply(best_chains, function(x) x$diagnostics),
      acceptance_rate = vapply(best_chains, function(x) x$acceptance_rate, numeric(1)),
      prop_matched = vapply(best_chains, function(x) x$prop_matched, numeric(1)),
      chain_seeds = chain_seeds_best
    ),
    draws = draws_out
  )
}
