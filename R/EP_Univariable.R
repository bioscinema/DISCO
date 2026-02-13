#' Univariate MEP Bayes with DISCO severity (EP_univariable)
#'
#' Runs a DISCO-based univariate separation diagnostic, constructs an adaptive
#' Multivariate Exponential Power (MEP) prior from the severity score, and fits
#' a univariate logistic model (intercept + one predictor) via random-walk
#' Metropolis–Hastings (RW–MH).
#'
#' @section Output conventions:
#' By default, the function returns the **standardized** coefficients:
#' \itemize{
#'   \item \code{beta0}: intercept in the standardized-\eqn{X} (z-scored) working scale.
#'   \item \code{beta1}: slope with respect to the standardized predictor.
#' }
#' Optionally, set \code{transform_beta1} to back-transform to the original
#' predictor scale, which adds:
#' \itemize{
#'   \item \code{beta0_orig}: intercept on the original predictor scale.
#'   \item \code{beta1_logit} \emph{or} \code{beta1_SAS} \emph{or} \code{beta1_Long}:
#'         slope per 1 unit change in the original predictor, on the chosen scale.
#' }
#'
#' @param data A data.frame containing \code{outcome} and \code{predictor}.
#' @param predictor String; name of the predictor column.
#' @param outcome String; name of the binary outcome column (default \code{"y"}).
#' @param missing How to handle missing data (applied \emph{once}, shared by severity
#'   and estimation): \code{"complete"} (drop rows with any NA in
#'   \code{predictor} or \code{outcome}) or \code{"impute"} (drop rows with NA in
#'   \code{outcome}, and impute \code{predictor} only). Default \code{"complete"}.
#' @param impute_args Optional list of imputation settings used only when
#'   \code{missing = "impute"}:
#'   \itemize{
#'     \item \code{numeric_method = "median"|"mean"} (default \code{"median"})
#'     \item \code{factor_method = "mode"} (default \code{"mode"})
#'   }
#'
#' @param n_iter Integer; total MCMC iterations (including burn-in).
#'   Default \code{20000}.
#' @param burn_in Integer; burn-in iterations discarded from the front.
#'   Default \code{5000}.
#' @param init_beta Initial value(s) for the MH chain on the standardized scale.
#'   Either a numeric vector of length 2 (intercept, slope) or a 2 by \code{n_chains}
#'   matrix for per-chain initial values. Default \code{c(0, 0)}.
#'
#' @param step_hi,step_lo RW–MH proposal s.d. blended by severity as
#'   \code{step = step_hi*(1 - severity) + step_lo*severity}.
#'   Defaults \code{0.30} and \code{0.12}.
#' @param ci_level Credible interval level in (0,1). Default \code{0.95}.
#' @param compare Logical; if \code{TRUE} (default) fit a GLM comparator on the same rows
#'   \emph{with standardized X} for reference.
#' @param return_draws Logical; if \code{TRUE}, include posterior draws on
#'   standardized and original scales. Default \code{FALSE}.
#' @param seed Optional integer; if provided, sets RNG seed for reproducibility.
#'
#' @param transform_beta1 One of \code{"none"}, \code{"logit"}, \code{"SAS"}, \code{"Long"}.
#'   If not \code{"none"}, also reports \code{beta0_orig} and the chosen slope on the
#'   original predictor scale. Default \code{"none"} (standardized-only output).
#'
#' @param sigma0 Prior sd for intercept (logit scale). Default \code{10}.
#' @param sigma1_hi Prior sd for slope under mild separation (\code{severity = 0}).
#'   Default \code{5}.
#' @param sigma1_lo Prior sd for slope under severe separation (\code{severity = 1}).
#'   Default \code{0.15}.
#' @param kappa_min,kappa_max Exponential-power shapes at \code{severity = 0} and
#'   \code{1}, blended linearly. Defaults \code{1} and \code{2.5}.
#'
#' @param tune_threshold_hi,tune_threshold_lo Burn-in acceptance thresholds for
#'   auto step-size tuning (increase if \code{> hi}; decrease if \code{< lo}).
#'   Defaults \code{0.45} and \code{0.20}.
#' @param tune_interval Iterations between tuning checks during burn-in.
#'   Default \code{500}.
#'
#' @param n_chains Integer; number of MH chains to run. Default \code{1}.
#' @param chain_seeds Optional integer vector of length \code{n_chains} giving
#'   per-chain RNG seeds. If \code{NULL}, seeds are generated deterministically
#'   from \code{seed}.
#' @param combine_chains How to combine chains for posterior summaries.
#'   \code{"stack"} binds post-burn draws across chains; \code{"none"} uses only
#'   the first chain. Default \code{"stack"}.
#'
#' @param ess_threshold Minimum effective sample size required (across parameters)
#'   to declare convergence when \pkg{coda} is available. Default \code{150}.
#' @param geweke_z_threshold Maximum allowed absolute Geweke z-score (across parameters)
#'   to declare convergence when \pkg{coda} is available. Default \code{2}.
#'
#' @param ci_levels_for_stars Numeric vector of nested credible interval levels
#'   used to assign “evidence stars” in \code{posterior}. Default
#'   \code{c(0.90, 0.95, 0.99)}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{predictor}, \code{outcome}
#'   \item \code{disco}: list with \code{separation_type}, \code{severity_score},
#'         \code{boundary_threshold}, \code{single_tie_boundary}, \code{missing_info}.
#'   \item \code{prior}: list with \code{mu}, \code{Sigma}, \code{kappa},
#'         \code{sigma0}, \code{sigma1} (slope prior scale implied by severity).
#'   \item \code{mcmc}: list with \code{acceptance_rate}, \code{step_size_used},
#'         \code{n_iter}, \code{burn_in}.
#'   \item \code{posterior}: data.frame of summaries with columns \code{Param},
#'         \code{Mean}, \code{SD}, \code{CI_low}, \code{CI_high}, plus
#'         \code{Sig_0} (CI excludes 0) and \code{Star} (based on nested credible
#'         intervals). Contains \code{beta0}, \code{beta1} (standardized) and, if
#'         requested via \code{transform_beta1}, \code{beta0_orig} plus one of
#'         \code{beta1_logit}/\code{beta1_SAS}/\code{beta1_Long}.
#'   \item \code{comparators}: list with a \code{glm} coefficient vector (fit on standardized X)
#'         when available.
#'   \item \code{rows_used}: integer indices of rows used after missing handling.
#'   \item \code{diagnostics_multi}: list with multi-chain diagnostics (when
#'         \code{n_chains >= 2} and \pkg{coda} is available): \code{rhat},
#'         \code{rhat_max}, \code{ess}, \code{ess_min}.
#'   \item \code{convergence}: list with single-posterior diagnostics for the
#'         final combined draws: \code{ess}, \code{geweke_z}, \code{ess_min},
#'         \code{geweke_max_abs}, \code{converged}.
#'   \item \code{mcmc}: list also includes per-chain acceptance rates when
#'         \code{n_chains > 1}.
#'   \item \code{draws} (optional): list with \code{chain_std} and \code{chain_orig}
#'         (included when \code{return_draws = TRUE}).
#' }
#'
#' @details
#' \strong{Shared missing handling.}
#' The \code{missing} choice is applied once to \code{(data[[outcome]], data[[predictor]])}.
#' With \code{missing = "complete"} we drop rows with any NA in those two columns.
#' With \code{missing = "impute"} we drop rows with NA in the outcome, then impute
#' NAs in the predictor using \code{impute_args}. The resulting rows/values are used
#' for both the DISCO severity and the model fit. We call
#' \code{DISCO::uni_separation(..., missing = "complete")} on the already-prepared data.
#'
#' \strong{Scaling for consistency.}
#' Numeric predictors are z-scored (\emph{standardized}) for: (i) the Bayesian fit,
#' (ii) the GLM comparator, and (iii) the DISCO severity computation.
#' Factors are left unchanged in the severity step; for estimation/comparators a
#' 2-level factor is converted to a 0/1 indicator (error if >2 levels).
#'
#' \strong{Adaptive prior from severity.}
#' Let \eqn{\bar{y}} be the sample mean on the analyzed rows and \eqn{s_x} the
#' predictor SD (on those rows). The MEP prior on \eqn{(\beta_0,\beta_1)} (working,
#' standardized-\eqn{X} scale) uses
#' \deqn{\mu = (\mathrm{logit}(\bar{y}),\, 0),\quad
#'       \Sigma = \mathrm{diag}(\sigma_0^2,\ \sigma_1^2),\quad
#'       \log \sigma_1 = (1-s)\log\sigma_{1,\mathrm{hi}} + s\log\sigma_{1,\mathrm{lo}},}
#' with \eqn{s} the severity score from \code{DISCO::uni_separation()} and
#' \eqn{\kappa = \kappa_{\min} + s(\kappa_{\max} - \kappa_{\min})}.
#'
#' \strong{Back-transform formulas.} If requested via \code{transform_beta1}:
#' \deqn{\beta_{1}^{\mathrm{logit}} = \beta_{1}^{\mathrm{std}}/s_x,\qquad
#'       \beta_{0}^{\mathrm{orig}} = \beta_{0}^{\mathrm{std}} - \beta_{1}^{\mathrm{std}}\cdot (\bar{x}/s_x),}
#' and we also report one of:
#' \deqn{\beta_{1}^{\mathrm{SAS}}  = \beta_{1}^{\mathrm{logit}} \cdot \pi/\sqrt{3},\qquad
#'       \beta_{1}^{\mathrm{Long}} = \beta_{1}^{\mathrm{logit}} \cdot (\pi/\sqrt{3} + 1).}
#'
#' \strong{Tuning.} During burn-in, the RW–MH step size is multiplicatively
#' adapted every \code{tune_interval} iterations to aim for an acceptance rate
#' between \code{tune_threshold_lo} and \code{tune_threshold_hi}.
#'
#' \strong{Initialization and multiple chains.}
#' The chain is initialized at \code{init_beta} on the standardized scale. If
#' \code{n_chains > 1}, chains are run with independent seeds and can be combined
#' via \code{combine_chains}.
#'
#' \strong{Convergence diagnostics.}
#' When \pkg{coda} is available, the function reports ESS and Geweke diagnostics
#' for the combined draws in \code{convergence}. For \code{n_chains >= 2}, it also
#' reports Gelman Rhat and multi-chain ESS in \code{diagnostics_multi}.
#'
#' \strong{Evidence stars.}
#' The \code{posterior} table includes \code{Sig_0} and \code{Star}, where stars
#' are assigned by whether nested credible intervals exclude 0.
#' @seealso
#' \code{\link[DISCO]{uni_separation}} for severity diagnostics.
#'
#' @importFrom DISCO uni_separation
#' @importFrom stats qlogis plogis quantile sd glm binomial coef median
#' @export
#'
#' @examples
#' \donttest{
#' ## Toy data
#' y <- c(0,0,0,0, 1,1,1,1)
#' x <- c(-0.52, -0.07, -0.60, -0.67, 1.39, 0.16, 1.40, 0.09)
#' df <- data.frame(y = y, x = x)
#'
#' ## 0) Detect Separation
#' detect <- DISCO::uni_separation(df, predictor = "x", outcome = "y")
#' detect$separation_type # e.g., "Perfect separation"
#'
#' ## 1) Default: STANDARDIZED coefficients (beta0, beta1)
#' fit_std <- EP_univariable(
#'   data = df, predictor = "x", outcome = "y",
#'   n_iter = 6000, burn_in = 2000, seed = 42
#' )
#' fit_std$posterior
#'
#' ## 2) Back-transform slope to ORIGINAL-x units on the LOGIT scale
#' fit_logit <- EP_univariable(
#'   data = df, predictor = "x", outcome = "y",
#'   transform_beta1 = "logit",
#'   n_iter = 6000, burn_in = 2000, seed = 42
#' )
#' fit_logit$posterior
#'
#' ## 3) Alternative effect scales on ORIGINAL-x units (choose ONE)
#' fit_sas <- EP_univariable(
#'   data = df, predictor = "x", outcome = "y",
#'   transform_beta1 = "SAS",
#'   n_iter = 6000, burn_in = 2000, seed = 42
#' )
#' fit_long <- EP_univariable(
#'   data = df, predictor = "x", outcome = "y",
#'   transform_beta1 = "Long",
#'   n_iter = 6000, burn_in = 2000, seed = 42
#' )
#' fit_sas$posterior   # contains beta1_SAS
#' fit_long$posterior  # contains beta1_Long
#'
#' ## 4) Shared missing handling
#' df2 <- df; df2$x[c(5, 8)] <- NA
#' fit_cc <- EP_univariable(df2, "x", "y", missing = "complete",
#'                          n_iter = 4000, burn_in = 1500, seed = 9)
#' fit_im <- EP_univariable(df2, "x", "y", missing = "impute",
#'                          impute_args = list(numeric_method = "median"),
#'                          n_iter = 4000, burn_in = 1500, seed = 9)
#'
#' ## 5) Multiple chains and Rhat (if coda installed)
#' fit_multi <- EP_univariable(
#'   data = df, predictor = "x", outcome = "y",
#'   n_iter = 20000, burn_in = 5000,
#'   init_beta = c(0, 0),
#'   n_chains = 4,
#'   chain_seeds = c(101, 102, 103, 104),
#'   combine_chains = "stack",
#'   seed = 9
#' )
#' fit_multi$diagnostics_multi$rhat
#' fit_multi$diagnostics_multi$rhat_max
#' fit_multi$convergence$converged

#' }
EP_univariable <- function(
    data,
    predictor,
    outcome = "y",
    missing = c("complete","impute"),
    impute_args = list(),
    n_iter = 20000,
    burn_in = 5000,
    init_beta = c(0, 0),
    step_hi = 0.30,
    step_lo = 0.12,
    ci_level = 0.95,
    n_chains = 1,
    chain_seeds = NULL,
    combine_chains = c("stack","none"),
    ess_threshold = 150,
    geweke_z_threshold = 2,
    compare = TRUE,
    return_draws = FALSE,
    seed = NULL,
    transform_beta1 = c("none","logit","SAS","Long"),
    sigma0 = 10,
    sigma1_hi = 5,
    sigma1_lo = 0.15,
    kappa_min = 1,
    kappa_max = 2.5,
    tune_threshold_hi = 0.45,
    tune_threshold_lo = 0.20,
    tune_interval = 500
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  missing <- match.arg(missing)
  transform_beta1 <- match.arg(transform_beta1)
  combine_chains <- match.arg(combine_chains)

  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  if (!predictor %in% names(data)) stop(sprintf("Predictor '%s' not found.", predictor), call. = FALSE)
  if (!outcome %in% names(data)) stop(sprintf("Outcome '%s' not found.", outcome), call. = FALSE)
  if (!is.null(seed)) set.seed(as.integer(seed))
  if (!is.numeric(n_iter) || length(n_iter) != 1 || n_iter <= 1) stop("`n_iter` must be > 1.", call. = FALSE)
  if (!is.numeric(burn_in) || length(burn_in) != 1 || burn_in < 0) stop("`burn_in` must be >= 0.", call. = FALSE)
  if (n_iter <= burn_in) stop("`n_iter` must be greater than `burn_in`.", call. = FALSE)

  if (!is.numeric(init_beta)) stop("`init_beta` must be numeric.", call. = FALSE)
  if (length(init_beta) == 1L) init_beta <- rep(as.numeric(init_beta), 2)
  if (length(init_beta) != 2L) stop("`init_beta` must be length 2 (intercept, slope).", call. = FALSE)
  init_beta <- as.numeric(init_beta)

  if (!is.numeric(n_chains) || length(n_chains) != 1 || n_chains < 1) stop("`n_chains` must be a positive integer.", call. = FALSE)
  n_chains <- as.integer(n_chains)

  if (!is.null(chain_seeds)) {
    if (length(chain_seeds) != n_chains) stop("`chain_seeds` must have length `n_chains`.", call. = FALSE)
    chain_seeds <- as.integer(chain_seeds)
  } else {
    base_seed <- if (!is.null(seed)) as.integer(seed) else sample.int(1e9, 1)
    chain_seeds <- base_seed + seq_len(n_chains)
  }

  impute_vec <- function(v, args) {
    if (is.numeric(v)) {
      if (anyNA(v)) {
        method <- tolower(args$numeric_method %||% "median")
        fill <- if (identical(method, "mean")) mean(v, na.rm = TRUE) else stats::median(v, na.rm = TRUE)
        if (!is.finite(fill)) fill <- 0
        v[is.na(v)] <- fill
      }
      v
    } else if (is.factor(v) || is.character(v) || is.logical(v)) {
      v <- as.factor(v)
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
      v
    } else {
      v <- suppressWarnings(as.numeric(v))
      v[is.na(v)] <- 0
      v
    }
  }

  # --- 1) Shared missing handling (applied once) ------------------------------
  y_full <- data[[outcome]]
  x_full <- data[[predictor]]

  if (missing == "complete") {
    keep_idx <- which(!is.na(y_full) & !is.na(x_full))
    if (!length(keep_idx)) stop("No complete cases for outcome & predictor.", call. = FALSE)
    y <- y_full[keep_idx]; x <- x_full[keep_idx]
  } else {
    keep_idx <- which(!is.na(y_full))
    if (!length(keep_idx)) stop("No rows with observed outcome.", call. = FALSE)
    y <- y_full[keep_idx]; x <- x_full[keep_idx]
    if (anyNA(x)) x <- impute_vec(x, impute_args)
  }

  # Normalize y to {0,1}
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y) || is.character(y)) y <- as.integer(factor(y)) - 1L
  y <- as.numeric(y)
  if (length(unique(y)) < 2L) stop("Outcome must contain both 0 and 1 after missing handling.", call. = FALSE)

  # Predictor for estimation: numeric vector (standardized later)
  if (!is.numeric(x)) {
    if (is.factor(x) || is.character(x) || is.logical(x)) {
      x <- as.factor(x)
      if (nlevels(x) != 2L) stop("Predictor has >2 levels; EP_univariable handles numeric or 2-level factors.", call. = FALSE)
      x <- as.integer(x) - 1L
    } else {
      x <- suppressWarnings(as.numeric(x))
    }
  }

  # --- 2) DISCO severity on the same rows; numeric predictor scaled -----------
  x_for_disco <- if (is.numeric(x)) as.numeric(scale(x)) else x
  res_disco <- DISCO::uni_separation(
    data = data.frame(y = y, x = x_for_disco),
    predictor = "x", outcome = "y",
    missing = "complete"
  )

  # --- 3) Standardize predictor for Bayesian and GLM comparator ---------------
  X_std <- scale(x)
  x_mean <- as.numeric(attr(X_std, "scaled:center"))
  x_sd   <- as.numeric(attr(X_std, "scaled:scale"))
  if (!is.finite(x_sd) || x_sd == 0) {
    return(list(
      predictor = predictor,
      outcome   = outcome,
      disco = list(
        separation_type = res_disco$separation_type,
        severity_score  = res_disco$severity_score,
        boundary_threshold = res_disco$boundary_threshold,
        single_tie_boundary = res_disco$single_tie_boundary,
        missing_info    = list(rows_used = keep_idx, policy = missing, imputed = (missing == "impute"))
      ),
      message = "Predictor has zero variance after missing-data handling; model not fit.",
      rows_used = keep_idx
    ))
  }
  X_std_num <- as.numeric(X_std)

  # --- 4) Adaptive prior from severity ----------------------------------------
  severity <- res_disco$severity_score
  if (is.null(severity) || is.na(severity)) severity <- 0
  y_bar <- mean(y)
  logit_clip <- function(p) stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))

  prior_from_severity <- function(severity_score, y_bar,
                                  sigma0, sigma1_hi, sigma1_lo,
                                  kappa_min, kappa_max) {
    mu0 <- logit_clip(y_bar)
    mu1 <- 0
    log_sigma1 <- (1 - severity_score) * log(sigma1_hi) + severity_score * log(sigma1_lo)
    sigma1 <- exp(log_sigma1)
    kappa  <- kappa_min + severity_score * (kappa_max - kappa_min)
    mu     <- c(mu0, mu1)
    Sigma  <- diag(c(sigma0^2, sigma1^2))
    list(mu = mu, Sigma = Sigma, kappa = kappa, sigma0 = sigma0, sigma1 = sigma1)
  }

  pr <- prior_from_severity(severity, y_bar, sigma0, sigma1_hi, sigma1_lo, kappa_min, kappa_max)
  mu <- pr$mu; Sigma <- pr$Sigma; kappa <- pr$kappa

  # Proposal step blended by severity
  step_size0 <- step_hi * (1 - severity) + step_lo * severity

  # --- 5) RW–MH sampler (one chain) -------------------------------------------
  log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))

  run_MH_uni_one <- function(n_iter, burn_in, init_beta, step_size, X_std, y, mu, Sigma, kappa,
                             tune_hi, tune_lo, tune_interval, chain_seed = NULL) {
    if (!is.null(chain_seed)) set.seed(as.integer(chain_seed))
    X <- cbind(Intercept = 1, X_std)
    Sigma_inv <- solve(Sigma)

    log_post <- function(beta) {
      eta <- as.vector(X %*% beta)
      ll <- sum(-y * log1pexp(-eta) - (1 - y) * log1pexp(eta))
      diff <- beta - mu
      d2 <- as.numeric(t(diff) %*% Sigma_inv %*% diff)
      ll - 0.5 * (d2^kappa)
    }

    chain <- matrix(0, nrow = n_iter, ncol = 2)
    colnames(chain) <- c("beta0", "beta1")
    chain[1, ] <- init_beta

    cur_lp <- log_post(chain[1, ])
    accept <- 0L
    ss <- step_size

    for (t in 2:n_iter) {
      prop <- chain[t - 1, ] + ss * rnorm(2)
      prop_lp <- log_post(prop)
      if (log(runif(1)) < (prop_lp - cur_lp)) {
        chain[t, ] <- prop
        cur_lp <- prop_lp
        accept <- accept + 1L
      } else {
        chain[t, ] <- chain[t - 1, ]
      }

      if (t <= burn_in && (t %% tune_interval) == 0) {
        acc <- accept / t
        if (acc > tune_hi) ss <- ss * 1.10
        if (acc < tune_lo) ss <- ss * 0.90
      }
    }

    post <- chain[(burn_in + 1):n_iter, , drop = FALSE]
    list(post = post, acceptance_rate = accept / (n_iter - 1), step_size_final = ss)
  }

  chains <- vector("list", n_chains)
  for (i in seq_len(n_chains)) {
    chains[[i]] <- run_MH_uni_one(
      n_iter = n_iter,
      burn_in = burn_in,
      init_beta = init_beta,
      step_size = step_size0,
      X_std = X_std_num,
      y = y,
      mu = mu,
      Sigma = Sigma,
      kappa = kappa,
      tune_hi = tune_threshold_hi,
      tune_lo = tune_threshold_lo,
      tune_interval = tune_interval,
      chain_seed = chain_seeds[i]
    )
  }

  # multi-chain diagnostics (Rhat and multi-chain ESS)
  diagnostics_multi <- list(
    rhat = c(beta0 = NA_real_, beta1 = NA_real_),
    rhat_max = NA_real_,
    ess = c(beta0 = NA_real_, beta1 = NA_real_),
    ess_min = NA_real_
  )

  if (requireNamespace("coda", quietly = TRUE) && n_chains >= 2) {
    mlist <- coda::mcmc.list(lapply(chains, function(z) coda::mcmc(z$post)))

    gd <- coda::gelman.diag(mlist, autoburnin = FALSE, multivariate = FALSE)$psrf
    # gd is a matrix with rows = params, columns include "Point est."
    rhat <- as.numeric(gd[, "Point est."])
    names(rhat) <- colnames(chains[[1]]$post)

    diagnostics_multi$rhat <- rhat
    diagnostics_multi$rhat_max <- suppressWarnings(max(rhat, na.rm = TRUE))

    ess_m <- coda::effectiveSize(mlist)
    diagnostics_multi$ess <- as.numeric(ess_m)
    names(diagnostics_multi$ess) <- names(rhat)
    diagnostics_multi$ess_min <- suppressWarnings(min(diagnostics_multi$ess, na.rm = TRUE))
  }

  # combine post-burn draws
  if (combine_chains == "stack") {
    post_std <- do.call(rbind, lapply(chains, `[[`, "post"))
  } else {
    post_std <- chains[[1]]$post
  }

  # --- 6) Back-transform chain to original predictor scale --------------------
  transform_chain_to_original <- function(chain_std, x_mean, x_sd) {
    b0_std <- chain_std[, 1]
    b1_std <- chain_std[, 2]
    beta1_logit <- b1_std / x_sd
    beta0_orig  <- b0_std - b1_std * (x_mean / x_sd)
    logit_sd <- pi / sqrt(3)
    beta1_SAS  <- beta1_logit * logit_sd
    beta1_Long <- beta1_logit * (logit_sd + 1)
    cbind(beta0_orig = beta0_orig,
          beta1_logit = beta1_logit,
          beta1_SAS = beta1_SAS,
          beta1_Long = beta1_Long)
  }
  post_orig <- transform_chain_to_original(post_std, x_mean, x_sd)

  # --- 7) Posterior summaries + significance stars ----------------------------
  ci_mat_for_level <- function(M, lvl) {
    qlo <- (1 - lvl) / 2
    qhi <- 1 - qlo
    t(apply(M, 2, stats::quantile, probs = c(qlo, qhi), na.rm = TRUE))
  }

  summarize_draws <- function(draws_mat, ci_level_main) {
    pm <- colMeans(draws_mat)
    sdv <- apply(draws_mat, 2, stats::sd)
    ci_main <- ci_mat_for_level(draws_mat, ci_level_main)

    ci90 <- ci_mat_for_level(draws_mat, 0.90)
    ci95 <- ci_mat_for_level(draws_mat, 0.95)
    ci99 <- ci_mat_for_level(draws_mat, 0.99)

    star <- rep("", ncol(draws_mat))
    star[(ci90[, 1] > 0) | (ci90[, 2] < 0)] <- "*"
    star[(ci95[, 1] > 0) | (ci95[, 2] < 0)] <- "**"
    star[(ci99[, 1] > 0) | (ci99[, 2] < 0)] <- "***"

    sig0 <- (ci_main[, 1] > 0) | (ci_main[, 2] < 0)

    data.frame(
      Mean = pm,
      SD = sdv,
      CI_low = ci_main[, 1],
      CI_high = ci_main[, 2],
      Sig_0 = as.logical(sig0),
      Star = star,
      check.names = FALSE
    )
  }

  # standardized params table
  tab_std <- summarize_draws(post_std, ci_level)
  posterior <- data.frame(
    Param = c("beta0", "beta1"),
    tab_std,
    row.names = NULL,
    check.names = FALSE
  )

  # optional original-scale rows
  if (transform_beta1 != "none") {
    tab_b0 <- summarize_draws(matrix(post_orig[, "beta0_orig"], ncol = 1), ci_level)
    row_b0 <- data.frame(Param = "beta0_orig", tab_b0, row.names = NULL, check.names = FALSE)

    if (transform_beta1 == "logit") {
      tab_b1 <- summarize_draws(matrix(post_orig[, "beta1_logit"], ncol = 1), ci_level)
      row_b1 <- data.frame(Param = "beta1_logit", tab_b1, row.names = NULL, check.names = FALSE)
    } else if (transform_beta1 == "SAS") {
      tab_b1 <- summarize_draws(matrix(post_orig[, "beta1_SAS"], ncol = 1), ci_level)
      row_b1 <- data.frame(Param = "beta1_SAS", tab_b1, row.names = NULL, check.names = FALSE)
    } else {
      tab_b1 <- summarize_draws(matrix(post_orig[, "beta1_Long"], ncol = 1), ci_level)
      row_b1 <- data.frame(Param = "beta1_Long", tab_b1, row.names = NULL, check.names = FALSE)
    }

    posterior <- rbind(posterior, row_b0, row_b1)
    rownames(posterior) <- NULL
  }

  # --- 8) Convergence diagnostics (coda if available) -------------------------
  convergence <- list(
    ess = rep(NA_real_, ncol(post_std)),
    geweke_z = rep(NA_real_, ncol(post_std)),
    ess_min = NA_real_,
    geweke_max_abs = NA_real_,
    converged = NA
  )

  if (requireNamespace("coda", quietly = TRUE)) {
    m <- coda::mcmc(post_std)
    ess <- as.numeric(coda::effectiveSize(m))
    gz  <- as.numeric(coda::geweke.diag(m)$z)

    convergence$ess <- ess
    convergence$geweke_z <- gz
    convergence$ess_min <- suppressWarnings(min(ess, na.rm = TRUE))
    convergence$geweke_max_abs <- suppressWarnings(max(abs(gz), na.rm = TRUE))

    ess_ok <- is.finite(convergence$ess_min) && convergence$ess_min >= ess_threshold
    gz_ok  <- is.finite(convergence$geweke_max_abs) && convergence$geweke_max_abs <= geweke_z_threshold
    convergence$converged <- isTRUE(ess_ok && gz_ok)
  }

  # --- 9) GLM comparator (standardized X) -------------------------------------
  comparators <- list(glm = NULL)
  if (isTRUE(compare)) {
    df_fit <- data.frame(y = y, X = X_std_num)
    glm_fit <- try(suppressWarnings(stats::glm(y ~ X, data = df_fit, family = stats::binomial())), silent = TRUE)
    if (!inherits(glm_fit, "try-error")) comparators$glm <- stats::coef(glm_fit)
  }

  # --- 10) Assemble and return ------------------------------------------------
  mcmc_info <- list(
    acceptance_rate_by_chain = vapply(chains, function(z) z$acceptance_rate, numeric(1)),
    step_size_final_by_chain = vapply(chains, function(z) z$step_size_final, numeric(1)),
    n_iter = n_iter,
    burn_in = burn_in,
    n_chains = n_chains,
    combine_chains = combine_chains,
    init_beta = init_beta,
    chain_seeds = chain_seeds
  )

  out <- list(
    predictor = predictor,
    outcome = outcome,
    disco = list(
      separation_type = res_disco$separation_type,
      severity_score = res_disco$severity_score,
      boundary_threshold = res_disco$boundary_threshold,
      single_tie_boundary = res_disco$single_tie_boundary,
      missing_info = list(rows_used = keep_idx, policy = missing, imputed = (missing == "impute"))
    ),
    prior = list(mu = mu, Sigma = Sigma, kappa = kappa,
                 sigma0 = pr$sigma0, sigma1 = pr$sigma1),
    mcmc = mcmc_info,
    posterior = posterior,
    convergence = convergence,
    diagnostics_multi = diagnostics_multi,
    comparators = comparators,
    rows_used = keep_idx
  )

  if (isTRUE(return_draws)) {
    out$draws <- list(
      chain_std = if (combine_chains == "stack") post_std else chains[[1]]$post,
      chain_orig = if (combine_chains == "stack") post_orig else transform_chain_to_original(chains[[1]]$post, x_mean, x_sd)
    )
  }

  out
}
