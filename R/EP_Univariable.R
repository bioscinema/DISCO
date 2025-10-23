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
#' @param missing How to handle missing data: \code{"complete"} (drop rows with any
#'   NA in \code{predictor} or \code{outcome}) or \code{"impute"} (impute the
#'   predictor only; outcome NA is always dropped). Default \code{"complete"}.
#' @param impute_args Optional list of imputation settings passed to
#'   \code{DISCO::uni_separation()} (e.g., \code{list(numeric_method="median")}).
#'
#' @param n_iter Integer; total MCMC iterations (including burn-in).
#'   Default \code{20000}.
#' @param burn_in Integer; burn-in iterations discarded from the front.
#'   Default \code{5000}.
#' @param step_hi,step_lo RW–MH proposal s.d. blended by severity as
#'   \code{step = step_hi*(1 - severity) + step_lo*severity}.
#'   Defaults \code{0.30} and \code{0.12}.
#' @param ci_level Credible interval level in (0,1). Default \code{0.95}.
#' @param compare Logical; if \code{TRUE} (default) fit GLM and Firth logistic
#'   (\pkg{logistf}) on the same rows for reference.
#' @param return_draws Logical; if \code{TRUE}, include posterior draws on
#'   standardized and original scales. Default \code{FALSE}.
#' @param seed For reproducibility, default 2025.
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
#'         \code{Mean}, \code{SD}, \code{CI_low}, \code{CI_high}, containing
#'         \code{beta0}, \code{beta1} (standardized) and, if requested via
#'         \code{transform_beta1}, \code{beta0_orig} plus one of
#'         \code{beta1_logit}/\code{beta1_SAS}/\code{beta1_Long}.
#'   \item \code{comparators}: list with \code{glm} and \code{firth} coefficient
#'         vectors when available.
#'   \item \code{rows_used}: integer indices of rows used after missing handling.
#'   \item \code{draws} (optional): list with \code{chain_std} and \code{chain_orig}
#'         (included when \code{return_draws = TRUE}).
#' }
#'
#' @details
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
#' @section Interpretation:
#' \itemize{
#'   \item \code{beta1} (default) — change in log-odds per 1 SD increase in the predictor (standardized scale).
#'   \item \code{beta1_logit} — change in log-odds per 1 unit increase in the original predictor.
#'   \item \code{beta1_SAS} and \code{beta1_Long} — alternative effect scalings derived from \code{beta1_logit}.
#'   \item \code{beta0_orig} is the intercept consistent with original predictor units.
#' }
#'
#' @seealso
#' \code{\link[DISCO]{uni_separation}} for severity diagnostics.
#'
#' @importFrom DISCO uni_separation
#' @importFrom stats qlogis plogis quantile sd glm binomial coef
#' @export
#'
#' @examples
#' \donttest{
#' ## Toy data: y ~ Bernoulli(logit^{-1}(-0.2 + 1.0 * x)) on ORIGINAL x units
#' set.seed(1)
#' n  <- 60
#' x  <- rnorm(n, mean = 0.3, sd = 1)
#' eta <- -0.2 + 1.0 * x
#' y  <- rbinom(n, size = 1, prob = stats::plogis(eta))
#' df <- data.frame(y = y, x = x)
#'
#' ## 1) Default: STANDARDIZED coefficients (beta0, beta1)
#' fit_std <- EP_univariable(
#'   data = df, predictor = "x", outcome = "y",
#'   n_iter = 6000, burn_in = 2000, seed = 42  # reduced iters for example speed
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
#' }
#'
EP_univariable <- function(
    data,
    predictor,
    outcome = "y",
    missing = c("complete","impute"),
    impute_args = list(),
    n_iter = 20000,
    burn_in = 5000,
    step_hi = 0.30,
    step_lo = 0.12,
    ci_level = 0.95,
    compare = TRUE,
    return_draws = FALSE,
    seed = 2025,
    # back-transform choice for beta1 (and corresponding intercept)
    transform_beta1 = c("none","logit","SAS","Long"),
    # prior defaults (aligned with your SLURM script)
    sigma0 = 10,
    sigma1_hi = 5,
    sigma1_lo = 0.15,
    kappa_min = 1,
    kappa_max = 2.5,
    # tuning defaults
    tune_threshold_hi = 0.45,
    tune_threshold_lo = 0.20,
    tune_interval = 500
) {
  missing <- match.arg(missing)
  transform_beta1 <- match.arg(transform_beta1)

  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  if (!predictor %in% names(data)) stop(sprintf("Predictor '%s' not found.", predictor), call. = FALSE)
  if (!outcome %in% names(data)) stop(sprintf("Outcome '%s' not found.", outcome), call. = FALSE)
  if (!is.null(seed)) set.seed(as.integer(seed))
  if (n_iter <= burn_in) stop("`n_iter` must be greater than `burn_in`.", call. = FALSE)

  # --- 1) DISCO severity + rows used -----------------------------------------
  res_disco <- DISCO::uni_separation(
    data = data, predictor = predictor, outcome = outcome,
    missing = missing, impute_args = impute_args
  )
  rows_used <- res_disco$missing_info$rows_used
  if (is.null(rows_used) || length(rows_used) == 0) {
    stop("No rows available after missing-data handling.", call. = FALSE)
  }

  df_used <- data[rows_used, , drop = FALSE]
  y <- df_used[[outcome]]
  x <- df_used[[predictor]]

  # Normalize outcome to {0,1}
  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y) || is.character(y)) y <- as.integer(factor(y)) - 1L
  y <- as.numeric(y)
  # Coerce predictor to numeric if needed
  if (!is.numeric(x)) x <- as.numeric(x)

  if (length(unique(y)) < 2L) {
    return(list(
      predictor = predictor,
      outcome   = outcome,
      disco = list(
        separation_type = res_disco$separation_type,
        severity_score  = res_disco$severity_score,
        boundary_threshold = res_disco$boundary_threshold,
        single_tie_boundary = res_disco$single_tie_boundary,
        missing_info    = res_disco$missing_info
      ),
      message = "Fewer than two outcome classes after missing-data handling; model not fit.",
      rows_used = rows_used
    ))
  }

  # --- 2) Standardize predictor ----------------------------------------------
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
        missing_info    = res_disco$missing_info
      ),
      message = "Predictor has zero variance after missing-data handling; model not fit.",
      rows_used = rows_used
    ))
  }

  # --- 3) Adaptive prior from severity ----------------------------------------
  severity <- res_disco$severity_score
  if (is.null(severity) || is.na(severity)) severity <- 0
  y_bar <- mean(y)

  logit_clip <- function(p) stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))

  prior_from_severity <- function(severity_score, y_bar,
                                  sigma0, sigma1_hi, sigma1_lo,
                                  kappa_min, kappa_max) {
    mu0 <- logit_clip(y_bar); mu1 <- 0
    log_sigma1 <- (1 - severity_score) * log(sigma1_hi) + severity_score * log(sigma1_lo)
    sigma1 <- exp(log_sigma1)
    kappa  <- kappa_min + severity_score * (kappa_max - kappa_min)
    mu     <- c(mu0, mu1)
    Sigma  <- diag(c(sigma0^2, sigma1^2))
    list(mu = mu, Sigma = Sigma, kappa = kappa,
         sigma0 = sigma0, sigma1 = sigma1)
  }

  pr <- prior_from_severity(severity, y_bar,
                            sigma0 = sigma0, sigma1_hi = sigma1_hi, sigma1_lo = sigma1_lo,
                            kappa_min = kappa_min, kappa_max = kappa_max)
  mu <- pr$mu; Sigma <- pr$Sigma; kappa <- pr$kappa

  # Proposal step blended by severity (milder → larger step)
  step_size <- step_hi * (1 - severity) + step_lo * severity

  # --- 4) RW–MH sampler (intercept + standardized slope) ----------------------
  run_MH_sampler_uni <- function(n_iter, init_beta, step_size, X_std, y, mu, Sigma, kappa, burn_in,
                                 tune_hi, tune_lo, tune_interval) {
    X <- cbind(Intercept = 1, X_std); p <- ncol(X)
    Sigma_inv <- solve(Sigma)
    log_post <- function(beta) {
      eta <- as.vector(X %*% beta)
      pi_x <- 1 / (1 + exp(-eta))
      loglik <- sum(y * log(pi_x + 1e-12) + (1 - y) * log(1 - pi_x + 1e-12))
      diff <- beta - mu
      d2 <- as.numeric(t(diff) %*% Sigma_inv %*% diff)
      logprior <- -0.5 * (d2^kappa)
      loglik + logprior
    }
    chain <- matrix(0, n_iter, p); chain[1, ] <- init_beta
    cur_lp <- log_post(chain[1, ]); accept <- 0L
    for (t in 2:n_iter) {
      prop <- chain[t - 1, ] + step_size * rnorm(p)
      prop_lp <- log_post(prop)
      if (log(runif(1)) < (prop_lp - cur_lp)) {
        chain[t, ] <- prop; cur_lp <- prop_lp; accept <- accept + 1L
      } else {
        chain[t, ] <- chain[t - 1, ]
      }
      if (t <= burn_in && (t %% tune_interval) == 0) {
        acc <- accept / t
        if (acc > tune_hi) step_size <- step_size * 1.10
        if (acc < tune_lo) step_size <- step_size * 0.90
      }
    }
    list(chain = chain[(burn_in + 1):n_iter, , drop = FALSE],
         acceptance_rate = accept / (n_iter - 1))
  }

  fit <- run_MH_sampler_uni(
    n_iter    = n_iter,
    init_beta = c(mu[1], 0),
    step_size = step_size,
    X_std     = as.numeric(X_std),
    y         = y,
    mu        = mu,
    Sigma     = Sigma,
    kappa     = kappa,
    burn_in   = burn_in,
    tune_hi   = tune_threshold_hi,
    tune_lo   = tune_threshold_lo,
    tune_interval = tune_interval
  )
  chain_std <- fit$chain
  acc_rate  <- fit$acceptance_rate

  # --- 5) Back-transform to original predictor scale --------------------------
  # beta1_logit = b1_std / s_x ;  beta0_orig = b0_std - b1_std * (mean_x / s_x)
  # SAS/Long are scalar multiples of beta1_logit
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
  chain_orig <- transform_chain_to_original(chain_std, x_mean, x_sd)

  # --- 6) Posterior summaries --------------------------------------------------
  qlo <- (1 - ci_level) / 2; qhi <- 1 - qlo
  summ_one <- function(v) {
    c(Mean = mean(v), SD = stats::sd(v),
      CI_low = as.numeric(stats::quantile(v, qlo)),
      CI_high = as.numeric(stats::quantile(v, qhi)))
  }

  # Always include standardized coefficients as beta0, beta1
  rows <- list(
    beta0 = chain_std[, 1],
    beta1 = chain_std[, 2]
  )

  # Optionally add the original-scale intercept and the chosen slope scale
  if (transform_beta1 != "none") {
    rows$beta0_orig <- chain_orig[, "beta0_orig"]
    if (transform_beta1 == "logit") {
      rows$beta1_logit <- chain_orig[, "beta1_logit"]
    } else if (transform_beta1 == "SAS") {
      rows$beta1_SAS <- chain_orig[, "beta1_SAS"]
    } else if (transform_beta1 == "Long") {
      rows$beta1_Long <- chain_orig[, "beta1_Long"]
    }
  }

  posterior <- do.call(rbind, lapply(names(rows), function(nm) {
    z <- summ_one(rows[[nm]])
    data.frame(Param = nm, t(z), row.names = NULL, check.names = FALSE)
  }))

  # --- 7) Frequentist comparators (optional) ----------------------------------
  comparators <- list(glm = NULL, firth = NULL)
  if (isTRUE(compare)) {
    df_fit <- data.frame(y = y, X = as.numeric(x))
    glm_fit <- try(suppressWarnings(stats::glm(y ~ X, data = df_fit, family = stats::binomial())), silent = TRUE)
    if (!inherits(glm_fit, "try-error")) comparators$glm <- stats::coef(glm_fit)
    firth_fit <- try(logistf::logistf(y ~ X, data = df_fit), silent = TRUE)
    if (!inherits(firth_fit, "try-error")) comparators$firth <- stats::coef(firth_fit)
  }

  # --- 8) Assemble and return --------------------------------------------------
  out <- list(
    predictor = predictor,
    outcome   = outcome,
    disco = list(
      separation_type = res_disco$separation_type,
      severity_score  = severity,
      boundary_threshold = res_disco$boundary_threshold,
      single_tie_boundary = res_disco$single_tie_boundary,
      missing_info    = res_disco$missing_info
    ),
    prior = list(mu = mu, Sigma = Sigma, kappa = kappa,
                 sigma0 = pr$sigma0, sigma1 = pr$sigma1),
    mcmc  = list(acceptance_rate = acc_rate, step_size_used = step_size,
                 n_iter = n_iter, burn_in = burn_in),
    posterior = posterior,
    comparators = comparators,
    rows_used = rows_used
  )
  if (isTRUE(return_draws)) out$draws <- list(chain_std = chain_std, chain_orig = chain_orig)
  out
}
