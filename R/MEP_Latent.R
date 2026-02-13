#' Severity-Adaptive MEP for Pure Latent (multi-predictor) Logistic — unified
#'
#' Fits a logistic model with a Multivariate Exponential Power (MEP) prior
#' via a random-walk Metropolis–Hastings (RW–MH) sampler and performs a small
#' grid search over prior settings \eqn{(\mu, \Sigma, \kappa)}. Predictors are
#' always z-scored internally for fitting (safe-scaling). Summaries and CIs are
#' reported on the working (logit/standardized) scale and back-transformed to
#' the original (encoded) predictor scale.
#'
#' @section What this function does:
#' \itemize{
#'   \item Handles missingness once (complete-case or imputation) on the
#'         \emph{raw} \code{X} and \code{y}, then uses the same rows for both
#'         fitting and the GLM reference.
#'   \item Encodes factors in \code{X} with \code{model.matrix(~ ., data = X)}
#'         (treatment contrasts, baseline = first level), drops the intercept,
#'         and fits on the numeric encoded design.
#'   \item Standardizes encoded predictors with a safe-scaler that guards against
#'         zero-variance columns (sets sd = 1 when sd = 0 or non-finite).
#'   \item For each grid setting, runs RW–MH with an EP prior and collects
#'         posterior means, credible intervals, and a posterior predictive
#'         “match” statistic.
#'   \item Selects a single run by (i) an acceptance-rate window, (ii) closeness
#'         to a GLM coefficient-ratio reference (computed on the same standardized
#'         working scale), and (iii) posterior predictive agreement.
#' }
#'
#' @section Factor handling & column names:
#' \itemize{
#'   \item \strong{Numeric predictors} appear as a single column with their
#'         original name (e.g., \code{X3}). No suffixes are added.
#'   \item \strong{Factor predictors} are expanded by \code{model.matrix()}
#'         using treatment contrasts with the \emph{first level as the baseline}.
#'         For a two-level factor \code{X3} with levels \code{A} and \code{B}
#'         (baseline = \code{A}), the encoded design includes a single dummy
#'         column \code{X3B}, which equals 1 when \code{X3 == "B"} and 0 when
#'         \code{X3 == "A"}. The reported effects are for these encoded columns.
#'   \item \strong{Change the baseline} beforehand to alter dummy labels:
#' \preformatted{X$X3 <- stats::relevel(X$X3, ref = "B")  # baseline becomes B; dummy shows as X3A}
#' }
#'
#' @param y Binary outcome (0/1, logical, or 2-level factor/character).
#' @param X Matrix or data.frame of predictors (no intercept). May include factors.
#'          Rows must align with \code{y}.
#'
#' @param missing How to handle missing data (applied once, shared by fitting and GLM):
#'   \code{"complete"} (drop rows with any NA in \code{y} or any column of \code{X}),
#'   or \code{"impute"} (drop rows with NA in \code{y}, then impute NAs in \code{X}
#'   \emph{before} factor encoding: numeric columns by mean/median; factor/character/logical
#'   columns by mode). Default \code{"complete"}.
#' @param impute_args Optional list of imputation settings used only when
#'   \code{missing = "impute"}:
#'   \itemize{
#'     \item \code{numeric_method = "median"|"mean"} (default \code{"median"})
#'     \item \code{factor_method  = "mode"}           (default \code{"mode"})
#'   }
#'
#' @param n_iter Integer; total MCMC iterations per grid point (including burn-in).
#'   Default \code{10000}.
#' @param burn_in Integer; burn-in iterations discarded from the front.
#'   Default \code{1000}.
#' @param init_beta Initial value(s) for the MH chain. Either a scalar (recycled)
#'   or a numeric vector of length \eqn{p} (intercept + encoded slopes).
#'   Default \code{0.01}.
#' @param step_size Proposal standard deviation for RW–MH. Default \code{0.40}.
#'
#' @param mu_vals Numeric vector; each value is repeated to length \eqn{p}
#'   (intercept + predictors) to form \eqn{\mu} in the prior grid.
#'
#' @param sigma0_intercept Prior scale for the intercept entry of \eqn{\Sigma}
#'   (logit scale). Default \code{10}.
#' @param sigma_global_multipliers Numeric vector of global multipliers applied
#'   to all \emph{slope} prior scales (intercept held at \code{sigma0_intercept}).
#'   For each multiplier \code{gm}, the prior scale matrix is
#'   \code{diag(c(sigma0_intercept, rep(gm, p-1)))}. Default
#'   \code{c(0.1, 0.5, 1, 2, 5, 10)}.
#'
#' @param kappa_vals Numeric vector of positive \eqn{\kappa} values (EP shape)
#'   to include in the grid. Default \code{c(0.5, 1, 2)}.
#'
#' @param accept_window Numeric length-2 vector giving the acceptable
#'   MH acceptance-rate window. Default \code{c(0.30, 0.40)}.
#' @param accept_target Numeric; target acceptance used to pick the closest run
#'   if no grid point falls in \code{accept_window}. Default \code{0.35}.
#'
#' @param ci_level Credible interval level in \eqn{(0,1)}. Default \code{0.95}.
#' @param ppc_threshold Posterior predictive match threshold; the returned
#'   \code{prop_matched} is the fraction of observations whose replicated
#'   outcomes match the observed at least this proportion across posterior
#'   draws. Default \code{0.80}.
#'
#' @param tune_threshold_hi,tune_threshold_lo Burn-in acceptance thresholds for
#'   multiplicative step-size tuning (increase if \code{> hi}; decrease if
#'   \code{< lo}). Defaults \code{0.45} and \code{0.20}.
#' @param tune_interval Iterations between tuning checks during burn-in.
#'   Default \code{500}.
#' @param verbose Logical; print brief progress messages. Default \code{FALSE}.
#' @param n_chains_best Integer; number of MH chains to rerun for the selected
#'   best grid point. Default \code{1}.
#' @param chain_seeds_best Optional integer vector of length \code{n_chains_best}
#'   giving per-chain RNG seeds for the best-point reruns. If \code{NULL}, seeds
#'   are generated deterministically from \code{seed}.
#' @param combine_chains How to combine the best-point chains for the final
#'   posterior summaries. \code{"stack"} binds post-burn draws across chains;
#'   \code{"none"} uses only the first chain. Default \code{"stack"}.
#' @param return_draws Logical; if \code{TRUE}, return post-burn draws for the
#'   selected best grid point (a matrix if one chain, else a list of matrices).
#' @param ess_threshold Minimum effective sample size required (across parameters)
#'   to declare convergence when \pkg{coda} is available. Default \code{150}.
#' @param geweke_z_threshold Maximum allowed absolute Geweke z-score (across parameters)
#'   to declare convergence when \pkg{coda} is available. Default \code{2}.
#' @param seed Optional integer; if provided, sets RNG seed for reproducibility.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{best_settings}: list with chosen \code{mu} (as a string),
#'         \code{Sigma_diag} (prior diagonal as a string), and \code{kappa}.
#'   \item \code{best_acceptance}: MH acceptance rate of the selected run.
#'   \item \code{best_prop_matched}: posterior predictive “match” statistic for the selected run.
#'   \item \code{posterior_means}: posterior means (length \eqn{p}) for the selected run (working scale).
#'   \item \code{standardized_coefs_back}: data.frame (per \emph{encoded} column) with means and CIs for:
#'         \code{Scaled} (working space; slope w.r.t. standardized predictor),
#'         and back-transformed effects on
#'         \code{b_A_original}, \code{b_SAS_original}, \code{b_Long_original}.
#'   \item \code{scaled_summary}: data.frame (including Intercept) with columns
#'         \code{Mean}, \code{SD}, \code{CI_low}, \code{CI_high} in the working
#'         (logit/standardized) space, plus \code{Sig} and \code{Star} where
#'         stars are computed from nested 90/95/99 credible intervals.
#'   \item \code{ci_level}: credible interval level used.
#'   \item \code{rows_used}: integer indices of rows kept after missing handling.
#'   \item \code{missing_info}: list with \code{policy} and \code{imputed} flags.
#'   \item \code{diagnostics_multi}: list with multi-chain diagnostics (when \code{n_chains_best >= 2}
#'         and \pkg{coda} is available): \code{rhat}, \code{rhat_max}, \code{ess}, \code{ess_min}.
#'   \item \code{convergence}: list with single-posterior diagnostics for the final combined draws:
#'         \code{ess}, \code{geweke_z}, \code{ess_min}, \code{geweke_max_abs}, \code{converged}.
#'   \item \code{best}: list including \code{n_chains_best}, \code{chain_seeds_best}, and \code{combine_chains}.

#' }
#'
#' @details
#' \strong{Missing handling.}
#' With \code{missing="complete"}, rows with any NA in \code{y} or \code{X} are dropped.
#' With \code{missing="impute"}, rows with NA in \code{y} are dropped, then each column of
#' \code{X} is imputed on the \emph{raw scale} before encoding: numeric columns by
#' \code{median} (default) or \code{mean}; factor/character/logical columns by their
#' \emph{mode} (most frequent level). After imputation, factors are encoded with
#' \code{model.matrix()}, and the sampler runs on the standardized encoded design.
#'
#' \strong{Initialization.} The sampler starts at the prior mean \code{mu} for each grid
#' point (no user-specified initial values). Each MH chain starts at \code{init_beta} (scalar recycled)
#' or a user-supplied vector of length \eqn{p}. When \code{n_chains_best > 1}, the
#' selected grid point is rerun with independent seeds.
#'
#' \strong{Sampler & prior.}
#' The (unnormalized) MEP prior is
#' \deqn{\pi(\beta \mid \mu,\Sigma,\kappa) \propto
#'       \exp\left\{-\tfrac12\left[(\beta-\mu)^\top \Sigma^{-1}(\beta-\mu)\right]^\kappa\right\}.}
#' We evaluate the quadratic form via a Cholesky of \eqn{\Sigma} and use stable
#' \code{log1p/exp} algebra for the Bernoulli log-likelihood.
#'
#' \strong{Standardization and back-transforms (encoded scale).}
#' Let \eqn{s_x} be the SD of an encoded column (numeric or 0/1 dummy) on the encoded
#' \code{X} scale, and \eqn{\beta^{\mathrm{std}}} the slope in the standardized design.
#' We report:
#' \deqn{b_{\mathrm{A}} = \beta^{\mathrm{std}}/s_x, \quad
#'       b_{\mathrm{SAS}} = b_{\mathrm{A}} \cdot \pi/\sqrt{3}, \quad
#'       b_{\mathrm{Long}} = b_{\mathrm{A}} \cdot (\pi/\sqrt{3} + 1).}
#' For a 0/1 dummy with prevalence \eqn{p}, \eqn{s_x = \sqrt{p(1-p)}}.
#'
#' \strong{Selection.}
#' Among all grid runs, candidates are those with acceptance in \code{accept_window}
#' (or closest to \code{accept_target} if none). When a GLM ratio is available, we prefer
#' small mean absolute deviation from that ratio; ties are broken by higher posterior
#' predictive agreement.
#'
#' \strong{Evidence flags (Sig and stars).}
#' The returned \code{scaled_summary} includes \code{Sig}, which is \code{TRUE}
#' when the (main) credible interval in the table excludes 0. In addition,
#' \code{Star} is computed from nested credible intervals using posterior draws:
#' \code{*} if the 90\% interval excludes 0, \code{**} if the 95\% interval excludes 0,
#' and \code{***} if the 99\% interval excludes 0.
#' These are descriptive summaries of posterior uncertainty, not frequentist p-values.
#'
#' @seealso
#' \code{\link{MEP_mixture}} for severity-anchored mixed predictors;
#' \code{\link{EP_univariable}} for a univariate EP Bayes fit with DISCO severity.
#'
#' @importFrom stats quantile rbinom glm binomial coef sd plogis median
#'
#' @examples
#' \donttest{
#' ## Numeric example
#' y <- c(0,0,0,0, 1,1,1,1)
#' X <- data.frame(
#'   X1 = c(-1.86, -0.81,  1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
#'   X2 = c( 0.52,  -0.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09)
#' )
#'
#' ## Latent (automatic): inclusion-minimal separating subsets of {X1, X2}
#' lat <- DISCO::latent_separation(
#'   y = y,
#'   X = X,
#'   find_minimal  = TRUE,
#'   mode          = "either",
#'   missing       = "complete",
#'   scale_X       = FALSE
#' )
#' names(lat$minimal_subsets)              # e.g., "X1_X2"
#' lat$minimal_subsets[[1]]$type           # "Perfect separation"
#' lat$minimal_subsets[[1]]$vars           # c("X1","X2")
#'
#' ## Complete-case
#' gs_cc <- MEP_latent(
#'   y = y, X = X,
#'   n_iter = 10000, burn_in = 1000, step_size = 0.4,
#'   mu_vals = seq(-1, 1, by = 0.2),
#'   sigma0_intercept = 10,
#'   sigma_global_multipliers = c(0.1, 0.5, 1, 2, 5),
#'   kappa_vals = c(0.5, 1, 2),
#'   missing = "complete",
#'   seed = 42
#' )
#' gs_cc$best_settings
#' gs_cc$scaled_summary
#' gs_cc$standardized_coefs_back
#'
#' ## Factor example (treatment coding; baseline = first level)
#' Xf <- data.frame(
#'   X1 = X$X1,
#'   G  = factor(c("A","A","B","B","A","B","A","B"), levels = c("A","B"))
#' )
#' fit_f <- MEP_latent(
#'   y = y, X = Xf,
#'   n_iter = 8000, burn_in = 800,
#'   mu_vals = seq(-0.5, 0.5, by = 0.25),
#'   sigma_global_multipliers = c(0.5, 1, 2),
#'   kappa_vals = c(1, 2),
#'   missing = "complete",
#'   seed = 99
#' )
#' fit_f$standardized_coefs_back     # contains column "GB" (= B vs A)
#'
#' ## Impute: add an NA to a factor and a numeric column
#' Xfi <- Xf
#' Xfi$X1[1] <- NA
#' Xfi$G[4]  <- NA
#' fit_i <- MEP_latent(
#'   y = y, X = Xfi,
#'   n_iter = 6000, burn_in = 600,
#'   mu_vals = seq(-0.5, 0.5, by = 0.25),
#'   sigma_global_multipliers = c(0.5, 1, 2),
#'   kappa_vals = c(1, 2),
#'   missing = "impute",
#'   impute_args = list(numeric_method = "median", factor_method = "mode"),
#'   seed = 77
#' )
#' fit_i$rows_used
#' head(fit_i$scaled_summary)
#' }
#' @export
MEP_latent <- function(
    y, X,
    missing = c("complete","impute"),
    impute_args = list(),
    n_iter = 10000, burn_in = 1000,
    init_beta = 0.01,
    step_size = 0.4,
    mu_vals = seq(-1, 1, by = 0.1),
    sigma0_intercept = 10,
    sigma_global_multipliers = c(0.1, 0.5, 1, 2, 5, 10),
    kappa_vals = c(0.5, 1, 2),
    accept_window = c(0.3, 0.4),
    accept_target = 0.35,
    ci_level = 0.95,
    ppc_threshold = 0.80,
    tune_threshold_hi = 0.45,
    tune_threshold_lo = 0.20,
    tune_interval = 500,
    verbose = FALSE,
    seed = NULL,
    n_chains_best = 1,
    chain_seeds_best = NULL,
    combine_chains = c("stack","none"),
    return_draws = FALSE,
    ess_threshold = 150,
    geweke_z_threshold = 2
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  missing <- base::match.arg(as.character(missing), choices = c("complete","impute"))
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

  impute_numeric <- function(v, method = "median") {
    if (!anyNA(v)) return(v)
    fill <- if (identical(tolower(method), "mean")) mean(v, na.rm = TRUE) else stats::median(v, na.rm = TRUE)
    if (!is.finite(fill)) fill <- 0
    v[is.na(v)] <- fill
    v
  }

  impute_factor_mode <- function(v) {
    v <- as.factor(v)
    if (!anyNA(v)) return(v)
    tab <- table(v, useNA = "no")
    if (length(tab)) {
      lvl <- names(tab)[which.max(tab)]
      v[is.na(v)] <- factor(lvl, levels = levels(v))
    } else {
      v <- factor(v, levels = c(levels(v), "Missing"))
      v[is.na(v)] <- "Missing"
    }
    v
  }

  safe_scale <- function(M) {
    cen <- suppressWarnings(colMeans(M))
    Xc  <- sweep(M, 2, cen, FUN = "-", check.margin = FALSE)
    sc  <- suppressWarnings(apply(M, 2, stats::sd))
    sc[!is.finite(sc) | sc == 0] <- 1
    Xs  <- sweep(Xc, 2, sc, FUN = "/", check.margin = FALSE)
    list(Xstd = Xs, center = cen, scale = sc)
  }

  log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))

  qfun_mat <- function(M, ci_level) {
    qlo <- (1 - ci_level) / 2
    qhi <- 1 - qlo
    t(apply(M, 2, stats::quantile, probs = c(qlo, qhi), na.rm = TRUE))
  }

  summarize_post <- function(post, X_orig, ci_level, ppc_threshold, y, Xw) {
    pm  <- colMeans(post)
    se  <- apply(post, 2, stats::sd)
    ci_scaled <- qfun_mat(post, ci_level)

    scaled_summary <- data.frame(
      Param   = colnames(Xw),
      Mean    = pm,
      SD      = se,
      CI_low  = ci_scaled[, 1],
      CI_high = ci_scaled[, 2],
      row.names = NULL,
      check.names = FALSE
    )

    # stars based on nested credible intervals from posterior draws
    ci_mat_for_level <- function(M, lvl) {
      qlo <- (1 - lvl) / 2
      qhi <- 1 - qlo
      t(apply(M, 2, stats::quantile, probs = c(qlo, qhi), na.rm = TRUE))
    }

    ci90 <- ci_mat_for_level(post, 0.90)
    ci95 <- ci_mat_for_level(post, 0.95)
    ci99 <- ci_mat_for_level(post, 0.99)

    star <- rep("", ncol(post))
    star[(ci90[, 1] > 0) | (ci90[, 2] < 0)] <- "*"
    star[(ci95[, 1] > 0) | (ci95[, 2] < 0)] <- "**"
    star[(ci99[, 1] > 0) | (ci99[, 2] < 0)] <- "***"

    # decide what "Sig" should mean
    # Option A: match the main CI level used in the table
    scaled_summary$Sig <- (scaled_summary$CI_low > 0) | (scaled_summary$CI_high < 0)

    # Option B: always use 95% for Sig (then rename it)
    # scaled_summary$Sig_95 <- (ci95[, 1] > 0) | (ci95[, 2] < 0)

    scaled_summary$Star <- star

    # back-transform slopes to original encoded predictor units
    X_orig_m <- as.matrix(X_orig)
    s_x <- apply(X_orig_m, 2, stats::sd)
    s_x[!is.finite(s_x) | s_x == 0] <- 1

    samp_scaled <- post[, -1, drop = FALSE]
    A_samp <- sweep(samp_scaled, 2, s_x, "/")

    logit_sd <- pi / sqrt(3)
    SAS_samp  <- A_samp * logit_sd
    Long_samp <- A_samp * (logit_sd + 1)

    summarise_mat <- function(M) {
      ci <- qfun_mat(M, ci_level)
      data.frame(Mean = colMeans(M), CI_low = ci[, 1], CI_high = ci[, 2])
    }

    out_scaled <- summarise_mat(samp_scaled)
    out_A      <- summarise_mat(A_samp)
    out_SAS    <- summarise_mat(SAS_samp)
    out_Long   <- summarise_mat(Long_samp)

    std_back <- data.frame(
      Predictor       = colnames(X_orig),
      Scaled          = out_scaled$Mean,
      Scaled_CI_low   = out_scaled$CI_low,
      Scaled_CI_high  = out_scaled$CI_high,
      b_A_original    = out_A$Mean,
      b_A_CI_low      = out_A$CI_low,
      b_A_CI_high     = out_A$CI_high,
      b_SAS_original  = out_SAS$Mean,
      b_SAS_CI_low    = out_SAS$CI_low,
      b_SAS_CI_high   = out_SAS$CI_high,
      b_Long_original = out_Long$Mean,
      b_Long_CI_low   = out_Long$CI_low,
      b_Long_CI_high  = out_Long$CI_high,
      row.names = NULL,
      check.names = FALSE
    )

    # posterior predictive check
    n_rep <- nrow(post)
    n_obs <- nrow(Xw)
    y_rep <- matrix(0L, nrow = n_rep, ncol = n_obs)
    for (i in seq_len(n_rep)) {
      pr_i <- stats::plogis(Xw %*% post[i, ])
      y_rep[i, ] <- stats::rbinom(n_obs, 1, pr_i)
    }
    p_match <- colMeans(sweep(y_rep, 2, y, `==`))
    prop_matched <- mean(p_match >= ppc_threshold, na.rm = TRUE)
    if (!is.finite(prop_matched)) prop_matched <- 0

    list(
      posterior_means = pm,
      scaled_summary = scaled_summary,
      standardized_coefs_back = std_back,
      prop_matched = prop_matched
    )
  }


  compute_convergence <- function(post, ess_threshold, geweke_z_threshold) {
    out <- list(
      ess = rep(NA_real_, ncol(post)),
      geweke_z = rep(NA_real_, ncol(post)),
      ess_min = NA_real_,
      geweke_max_abs = NA_real_,
      converged = NA
    )
    if (requireNamespace("coda", quietly = TRUE)) {
      m <- coda::mcmc(post)
      ess <- as.numeric(coda::effectiveSize(m))
      gz  <- as.numeric(coda::geweke.diag(m)$z)

      out$ess <- ess
      out$geweke_z <- gz
      out$ess_min <- suppressWarnings(min(ess, na.rm = TRUE))
      out$geweke_max_abs <- suppressWarnings(max(abs(gz), na.rm = TRUE))

      ess_ok <- is.finite(out$ess_min) && out$ess_min >= ess_threshold
      gz_ok  <- is.finite(out$geweke_max_abs) && out$geweke_max_abs <= geweke_z_threshold
      out$converged <- isTRUE(ess_ok && gz_ok)
    }
    out
  }

  # normalize y -> {0,1}
  y_full <- y
  if (!is.numeric(y_full)) {
    if (is.logical(y_full)) y_full <- as.integer(y_full)
    else if (is.factor(y_full) || is.character(y_full)) y_full <- as.integer(factor(y_full)) - 1L
    else stop("`y` must be numeric 0/1, logical, or 2-level factor/character.", call. = FALSE)
  }
  y_full <- as.numeric(y_full)

  X_raw <- as.data.frame(X, stringsAsFactors = FALSE)
  if (nrow(X_raw) != length(y_full)) stop("Rows of X must match length of y.", call. = FALSE)
  if (ncol(X_raw) < 1L) stop("X must have at least one predictor.", call. = FALSE)

  # shared missing handling on raw (before encoding)
  if (missing == "complete") {
    keep_idx <- which(stats::complete.cases(data.frame(y = y_full, X_raw, check.names = FALSE)))
    if (!length(keep_idx)) stop("No complete rows for y and X.", call. = FALSE)
    y_used <- y_full[keep_idx]
    X_used <- X_raw[keep_idx, , drop = FALSE]
  } else {
    keep_idx <- which(!is.na(y_full))
    if (!length(keep_idx)) stop("No rows with observed y.", call. = FALSE)
    y_used <- y_full[keep_idx]
    X_used <- X_raw[keep_idx, , drop = FALSE]

    num_method <- impute_args$numeric_method %||% "median"
    fac_method <- impute_args$factor_method  %||% "mode"

    for (j in seq_len(ncol(X_used))) {
      v <- X_used[[j]]
      if (is.numeric(v)) {
        X_used[[j]] <- impute_numeric(v, method = num_method)
      } else if (is.factor(v) || is.character(v) || is.logical(v)) {
        if (!identical(tolower(fac_method), "mode")) {
          warning("Only factor_method = 'mode' is supported here; using mode.")
        }
        X_used[[j]] <- impute_factor_mode(v)
      } else {
        v2 <- suppressWarnings(as.numeric(v))
        X_used[[j]] <- impute_numeric(v2, method = num_method)
      }
    }
  }

  if (length(y_used) < 2L || length(unique(y_used)) < 2L) {
    stop("Outcome must contain both 0 and 1 after missing handling.", call. = FALSE)
  }

  # encode factors, drop intercept
  mm <- stats::model.matrix(~ ., data = X_used)
  X_enc <- mm[, -1, drop = FALSE]
  X_mat <- as.matrix(X_enc)
  p_enc <- ncol(X_mat)
  p_all <- 1L + p_enc

  # GLM ratio reference on standardized encoded design (same rows)
  S_ref <- safe_scale(X_mat)
  X_ref <- S_ref$Xstd
  df_ref <- data.frame(y = y_used, X_ref)
  glm_fit <- try(suppressWarnings(stats::glm(y ~ ., data = df_ref, family = stats::binomial())), silent = TRUE)
  glm_ok <- !(inherits(glm_fit, "try-error"))
  glm_coefs <- if (glm_ok) stats::coef(glm_fit) else rep(NA_real_, p_all)

  parse_ratio <- function(x) {
    if (is.null(x) || is.na(x) || !nzchar(x)) return(NA_real_)
    as.numeric(strsplit(x, ",\\s*")[[1]])
  }

  Ref_ratio <- if (length(glm_coefs) >= 3 && is.finite(glm_coefs[2]) && abs(glm_coefs[2]) > 0) {
    paste(round(glm_coefs[-c(1, 2)] / glm_coefs[2], 3), collapse = ", ")
  } else {
    NA_character_
  }
  ref_ratio_vec <- parse_ratio(Ref_ratio)

  # build grids
  mu_grid <- lapply(mu_vals, function(m) rep(m, p_all))
  build_sigma <- function(gm) diag(c(sigma0_intercept, rep(gm, p_all - 1L)), nrow = p_all, ncol = p_all)
  Sigma_list <- lapply(sigma_global_multipliers, build_sigma)
  kappa_grid <- kappa_vals

  # sampler for a single chain at one grid point
  run_chain_one <- function(n_iter, burn_in, init_beta, step_size,
                            X_orig, y, mu, Sigma, kappa,
                            ci_level, ppc_threshold,
                            tune_threshold_hi, tune_threshold_lo, tune_interval,
                            verbose, chain_seed = NULL) {

    if (!is.null(chain_seed)) set.seed(as.integer(chain_seed))

    S <- safe_scale(X_orig)
    X_work <- S$Xstd
    Xw <- cbind(Intercept = 1, X_work)
    colnames(Xw) <- c("Intercept", colnames(X_orig))

    cholSigma <- try(chol(Sigma), silent = TRUE)
    if (inherits(cholSigma, "try-error")) stop("Sigma is not positive-definite.", call. = FALSE)

    qform <- function(v) {
      z <- backsolve(cholSigma, v, transpose = TRUE)
      sum(z * z)
    }

    # init vector
    if (length(init_beta) == 1L) init_vec <- rep(as.numeric(init_beta), p_all)
    else if (length(init_beta) == p_all) init_vec <- as.numeric(init_beta)
    else stop("`init_beta` must be a scalar or a numeric vector of length (1 + p_enc).", call. = FALSE)

    log_post <- function(beta) {
      eta <- as.vector(Xw %*% beta)
      ll  <- sum(-y * log1pexp(-eta) - (1 - y) * log1pexp(eta))
      diff <- beta - mu
      ll - 0.5 * (qform(diff)^kappa)
    }

    chain <- matrix(0, nrow = n_iter, ncol = p_all)
    colnames(chain) <- colnames(Xw)
    chain[1, ] <- init_vec
    cur_lp <- log_post(chain[1, ])
    accept <- 0L
    ss <- step_size

    if (isTRUE(verbose)) {
      message("RW-MH: ", n_iter, " iters; burn-in ", burn_in)
    }

    for (it in 2:n_iter) {
      prop <- chain[it - 1, ] + ss * rnorm(p_all)
      prop_lp <- log_post(prop)

      if (log(runif(1)) < (prop_lp - cur_lp)) {
        chain[it, ] <- prop
        cur_lp <- prop_lp
        accept <- accept + 1L
      } else {
        chain[it, ] <- chain[it - 1, ]
      }

      if (it <= burn_in && (it %% tune_interval) == 0) {
        ar <- accept / it
        if (ar > tune_threshold_hi) ss <- ss * 1.10
        if (ar < tune_threshold_lo) ss <- ss * 0.90
      }
    }

    post <- chain[(burn_in + 1):n_iter, , drop = FALSE]
    list(
      post = post,
      acceptance_rate = accept / (n_iter - 1),
      scaling_info = list(center = S$center, scale = S$scale),
      Xw = Xw
    )
  }

  # grid search (single chain per grid point)
  results_df <- data.frame(
    grid_id = integer(),
    mu = character(),
    sigma_diag = character(),
    kappa = numeric(),
    acceptance_rate = numeric(),
    prop_matched = numeric(),
    posterior_ratio_scaled = character(),
    stringsAsFactors = FALSE
  )
  runs <- list()
  gid <- 1L

  for (mu in mu_grid) {
    mu_str <- paste(round(mu, 6), collapse = ", ")
    for (Sigma in Sigma_list) {
      sigma_diag_str <- paste(round(diag(Sigma), 6), collapse = ", ")
      for (kappa in kappa_grid) {

        grid_chain_seed <- if (!is.null(seed)) as.integer(seed) + gid else NULL

        ch <- run_chain_one(
          n_iter = n_iter,
          burn_in = burn_in,
          init_beta = init_beta,
          step_size = step_size,
          X_orig = X_mat,
          y = y_used,
          mu = mu,
          Sigma = Sigma,
          kappa = kappa,
          ci_level = ci_level,
          ppc_threshold = ppc_threshold,
          tune_threshold_hi = tune_threshold_hi,
          tune_threshold_lo = tune_threshold_lo,
          tune_interval = tune_interval,
          verbose = verbose,
          chain_seed = grid_chain_seed
        )

        sum_one <- summarize_post(
          post = ch$post,
          X_orig = X_mat,
          ci_level = ci_level,
          ppc_threshold = ppc_threshold,
          y = y_used,
          Xw = ch$Xw
        )

        vals <- sum_one$posterior_means
        ratio_str <- if (length(vals) >= 3 && is.finite(vals[2]) && abs(vals[2]) > 0) {
          paste(round(vals[-c(1, 2)] / vals[2], 3), collapse = ", ")
        } else {
          NA_character_
        }

        results_df <- rbind(
          results_df,
          data.frame(
            grid_id = gid,
            mu = mu_str,
            sigma_diag = sigma_diag_str,
            kappa = kappa,
            acceptance_rate = ch$acceptance_rate,
            prop_matched = sum_one$prop_matched,
            posterior_ratio_scaled = ratio_str,
            stringsAsFactors = FALSE
          )
        )

        runs[[gid]] <- list(
          grid_id = gid,
          mu = mu,
          mu_str = mu_str,
          Sigma = Sigma,
          sigma_diag_str = sigma_diag_str,
          kappa = kappa,
          acceptance_rate = ch$acceptance_rate,
          prop_matched = sum_one$prop_matched,
          posterior_means = sum_one$posterior_means,
          posterior_ratio_scaled = ratio_str
        )

        gid <- gid + 1L
      }
    }
  }

  # selection
  lo <- accept_window[1]
  hi <- accept_window[2]
  cand <- results_df[
    is.finite(results_df$acceptance_rate) &
      results_df$acceptance_rate >= lo &
      results_df$acceptance_rate <= hi,
    , drop = FALSE
  ]

  if (nrow(cand) == 0) {
    idx <- which.min(abs(results_df$acceptance_rate - accept_target))
    cand <- results_df[idx, , drop = FALSE]
  }

  cand$ratio_term <- NA_real_
  if (!all(is.na(ref_ratio_vec))) {
    cand$ratio_term <- vapply(cand$posterior_ratio_scaled, function(r) {
      pv <- parse_ratio(r)
      if (all(is.na(pv))) return(NA_real_)
      if (length(pv) != length(ref_ratio_vec)) return(NA_real_)
      mean(abs(pv - ref_ratio_vec))
    }, numeric(1))
  }

  # score: prefer lower ratio_term if defined, else prefer higher prop_matched
  if (any(is.finite(cand$ratio_term))) {
    cand <- cand[order(cand$ratio_term, -cand$prop_matched), , drop = FALSE]
  } else {
    cand <- cand[order(-cand$prop_matched), , drop = FALSE]
  }

  best_id <- cand$grid_id[1]
  best_run <- runs[[best_id]]
  mu_best <- best_run$mu
  Sigma_best <- best_run$Sigma
  kappa_best <- best_run$kappa

  # rerun best with multiple chains
  best_chains <- vector("list", n_chains_best)
  for (ch_i in seq_len(n_chains_best)) {
    ch <- run_chain_one(
      n_iter = n_iter,
      burn_in = burn_in,
      init_beta = init_beta,
      step_size = step_size,
      X_orig = X_mat,
      y = y_used,
      mu = mu_best,
      Sigma = Sigma_best,
      kappa = kappa_best,
      ci_level = ci_level,
      ppc_threshold = ppc_threshold,
      tune_threshold_hi = tune_threshold_hi,
      tune_threshold_lo = tune_threshold_lo,
      tune_interval = tune_interval,
      verbose = verbose,
      chain_seed = chain_seeds_best[ch_i]
    )

    sum_one <- summarize_post(
      post = ch$post,
      X_orig = X_mat,
      ci_level = ci_level,
      ppc_threshold = ppc_threshold,
      y = y_used,
      Xw = ch$Xw
    )

    best_chains[[ch_i]] <- list(
      post = ch$post,
      acceptance_rate = ch$acceptance_rate,
      prop_matched = sum_one$prop_matched
    )
  }

  best_acceptance <- mean(vapply(best_chains, function(x) x$acceptance_rate, numeric(1)), na.rm = TRUE)
  best_prop_matched <- mean(vapply(best_chains, function(x) x$prop_matched, numeric(1)), na.rm = TRUE)

  # multi-chain diagnostics (Rhat and multi-chain ESS)
  diagnostics_multi <- list(
    rhat = rep(NA_real_, ncol(best_chains[[1]]$post)),
    rhat_max = NA_real_,
    ess = rep(NA_real_, ncol(best_chains[[1]]$post)),
    ess_min = NA_real_
  )
  if (requireNamespace("coda", quietly = TRUE) && n_chains_best >= 2) {
    mlist <- coda::mcmc.list(lapply(best_chains, function(x) coda::mcmc(x$post)))
    gd <- coda::gelman.diag(mlist, autoburnin = FALSE, multivariate = FALSE)$psrf
    diagnostics_multi$rhat <- as.numeric(gd[, "Point est."])
    diagnostics_multi$rhat_max <- suppressWarnings(max(diagnostics_multi$rhat, na.rm = TRUE))
    ess_m <- coda::effectiveSize(mlist)
    diagnostics_multi$ess <- as.numeric(ess_m)
    diagnostics_multi$ess_min <- suppressWarnings(min(diagnostics_multi$ess, na.rm = TRUE))
  }

  # combine chains for final summaries
  if (combine_chains == "stack") {
    post_all <- do.call(rbind, lapply(best_chains, `[[`, "post"))
  } else {
    post_all <- best_chains[[1]]$post
  }

  # Rebuild Xw for summaries
  S_final <- safe_scale(X_mat)
  X_work_final <- S_final$Xstd
  Xw_final <- cbind(Intercept = 1, X_work_final)
  colnames(Xw_final) <- c("Intercept", colnames(X_mat))

  sum_final <- summarize_post(
    post = post_all,
    X_orig = X_mat,
    ci_level = ci_level,
    ppc_threshold = ppc_threshold,
    y = y_used,
    Xw = Xw_final
  )

  convergence <- compute_convergence(
    post = post_all,
    ess_threshold = ess_threshold,
    geweke_z_threshold = geweke_z_threshold
  )

  draws_out <- NULL
  if (isTRUE(return_draws)) {
    if (n_chains_best == 1L) draws_out <- best_chains[[1]]$post
    else draws_out <- lapply(best_chains, `[[`, "post")
  }

  list(
    best_settings = list(
      mu = best_run$mu_str,
      Sigma_diag = best_run$sigma_diag_str,
      kappa = kappa_best
    ),
    best_acceptance = best_acceptance,
    best_prop_matched = best_prop_matched,
    posterior_means = sum_final$posterior_means,
    standardized_coefs_back = sum_final$standardized_coefs_back,
    scaled_summary = sum_final$scaled_summary,
    convergence = convergence,
    diagnostics_multi = diagnostics_multi,
    rows_used = keep_idx,
    missing_info = list(policy = missing, imputed = identical(missing, "impute")),
    best = list(
      grid_id = best_id,
      acceptance_rate = best_acceptance,
      prop_matched = best_prop_matched,
      n_chains_best = n_chains_best,
      chain_seeds_best = chain_seeds_best,
      combine_chains = combine_chains,
      converged = convergence$converged,
      ess_min = convergence$ess_min,
      geweke_max_abs = convergence$geweke_max_abs
    ),
    draws = draws_out
  )
}
