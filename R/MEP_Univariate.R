#' Univariate Multiple adaptive Exponential power procedure Bayes with DISCO severity (MEP_Univariate)
#'
#' Runs a DISCO-based univariate separation diagnostic, constructs a
#' Multiple adaptive Exponential power procedure prior from the severity score,
#' and fits a univariate logistic model (intercept plus one predictor) via
#' random-walk Metropolis-Hastings (RW-MH).
#'
#' Missing handling is complete-case only: rows with any NA in the outcome or
#' predictor are dropped once and reused for both the DISCO severity diagnostic
#' and Bayesian estimation.
#'
#' @section Output conventions:
#' The function reports posterior summaries for the slope coefficient of the predictor:
#' \itemize{
#'   \item \code{beta1}: slope per 1 SD increase in the predictor (standardized scale).
#' }
#' Optionally, set \code{transform_beta} to also report the slope on the original
#' predictor scale:
#' \itemize{
#'   \item One of \code{beta1_logit}, \code{beta1_SAS}, or \code{beta1_Long}.
#' }
#'
#' @param data A data.frame containing \code{outcome} and \code{predictor}.
#' @param predictor String; name of the predictor column.
#' @param outcome String; name of the binary outcome column (default \code{"y"}).
#'
#' @param burn_in Integer; number of burn-in iterations per chain (discarded). Default \code{5000}.
#' @param n_iter Integer; number of post-burn-in MCMC iterations (posterior draws) per chain.
#' Default \code{20000}.
#' @param init_beta Numeric vector of length 2 giving initial values \code{c(beta0, beta1)}
#' on the standardized scale. Default \code{c(0, 0)}.
#'
#' @param step_hi,step_lo RW-MH proposal standard deviation blended by severity:
#' \code{step = step_hi*(1 - severity) + step_lo*severity}. Defaults \code{0.30} and \code{0.12}.
#' @param tune_threshold_hi,tune_threshold_lo Burn-in acceptance thresholds for auto step-size tuning.
#' Increase step size if acceptance \code{> tune_threshold_hi}; decrease if \code{< tune_threshold_lo}.
#' Defaults \code{0.45} and \code{0.20}.
#' @param tune_interval Iterations between tuning checks during burn-in. Default \code{500}.
#'
#' @param ci_level Credible interval level in (0,1). Default \code{0.95}.
#' @param ci_levels_for_stars Numeric vector of nested credible interval levels used to assign
#' evidence stars in \code{posterior}. Default \code{c(0.90, 0.95, 0.99)}.
#'
#' @param n_chains Integer; number of MH chains to run. Default \code{1}.
#' @param chain_seeds Optional integer vector of length \code{n_chains} giving per-chain RNG seeds.
#' If \code{NULL}, random seeds are generated.
#' @param combine_chains How to combine chains for posterior summaries.
#' \code{"stack"} binds post-burn draws across chains; \code{"none"} uses only the first chain.
#' Default \code{"stack"}.
#'
#' @param ess_threshold Minimum effective sample size required (across parameters) to declare
#' convergence for the single-chain case when \pkg{coda} is available. Default \code{150}.
#' @param geweke_z_threshold Maximum allowed absolute Geweke z-score (across parameters) to declare
#' convergence for the single-chain case when \pkg{coda} is available. Default \code{2}.
#'
#' @param compare Logical; if \code{TRUE} (default) fit a GLM comparator on the same rows
#' with standardized \eqn{X}.
#' @param return_draws Logical; if \code{TRUE}, include post-burn posterior draws on
#' standardized and original scales. Default \code{FALSE}.
#'
#' @param transform_beta One of \code{"none"}, \code{"logit"}, \code{"SAS"}, \code{"Long"}.
#' If not \code{"none"}, also reports the chosen slope on the original predictor scale.
#' Default \code{"none"}.
#'
#' @param sigma0 Prior sd for intercept (logit scale). Default \code{10}.
#' @param sigma1_hi Prior sd for slope under mild separation (severity near 0). Default \code{5}.
#' @param sigma1_lo Prior sd for slope under severe separation (severity near 1). Default \code{0.15}.
#' @param kappa_min,kappa_max Exponential-power shapes at severity 0 and 1, blended linearly.
#' Defaults \code{1} and \code{2.5}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{predictor}, \code{outcome}.
#'   \item \code{disco}: list with \code{separation_type}, \code{severity_score},
#'         \code{boundary_threshold}, \code{single_tie_boundary}, and \code{missing_info}.
#'   \item \code{prior}: list with \code{mu}, \code{Sigma}, \code{kappa}, \code{sigma0}, \code{sigma1}.
#'   \item \code{mcmc}: per-chain acceptance rates, step sizes, tuning traces, plus
#'         \code{burn_in}, \code{n_iter} (post-burn), \code{n_total} (burn_in + n_iter),
#'         \code{n_chains}, and \code{combine_chains}.
#'   \item \code{posterior}: data.frame with \code{Param}, \code{Mean}, \code{SD},
#'         \code{CI_low}, \code{CI_high}, \code{Sig_0}, and \code{Star}.
#'   \item \code{comparators}: list with GLM coefficients when \code{compare = TRUE}.
#'   \item \code{rows_used}: integer indices of rows used after missing handling.
#'   \item \code{diagnostics_multi}: included only when \code{n_chains >= 2} and \pkg{coda} is available.
#'   \item \code{diagnostics_single}: included only when \code{n_chains == 1} and \pkg{coda} is available.
#'   \item \code{draws} (optional): included when \code{return_draws = TRUE}.
#' }
#'
#' @importFrom stats qlogis quantile sd glm binomial coef median rnorm runif
#' @export
#'
#' @examples
#' \donttest{
#' y <- c(0,0,0,0, 1,1,1,1)
#' x <- c(-0.52, -0.07, -0.60, -0.67, 1.39, 0.16, 1.40, 0.09)
#' df <- data.frame(y = y, x = x)
#'
#'## 1) Default: STANDARDIZED coefficients for predictor
#' fit_std <- MEP_Univariate(
#'   data = df, predictor = "x", outcome = "y"
#' )
#' fit_std$posterior
#' fit_std$diagnostics_single
#'
#'## 2) Back-transform slope to ORIGINAL-x units on the LOGIT scale
#' fit_logit <- MEP_Univariate(
#'   data = df, predictor = "x", outcome = "y",
#'   transform_beta = "logit"
#' )
#' fit_logit$posterior
#' fit_logit$diagnostics_single
#'
#'## 3) Multiple chains
#' fit_multi <- MEP_Univariate(
#'   data = df, predictor = "x", outcome = "y",
#'   n_chains = 4,
#'   chain_seeds = c(101, 102, 103, 104),
#'   combine_chains = "stack"
#' )
#' fit_multi$posterior
#' fit_multi$diagnostics_multi
#' }
MEP_Univariate <- function(
    data,
    predictor,
    outcome = "y",
    burn_in = 5000,          # burn-in iterations per chain (discarded)
    n_iter = 15000,          # post-burn draws per chain
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
    return_draws = TRUE,
    transform_beta = "none",
    sigma0 = 10,
    sigma1_hi = 5,
    sigma1_lo = 0.15,
    kappa_min = 1,
    kappa_max = 2.5,
    tune_threshold_hi = 0.45,
    tune_threshold_lo = 0.20,
    tune_interval = 500,
    ci_levels_for_stars = c(0.90, 0.95, 0.99)
) {
  transform_beta <- match.arg(transform_beta, choices = c("none","logit","SAS","Long"))
  combine_chains <- match.arg(combine_chains)

  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  if (!predictor %in% names(data)) stop(sprintf("Predictor '%s' not found.", predictor), call. = FALSE)
  if (!outcome %in% names(data)) stop(sprintf("Outcome '%s' not found.", outcome), call. = FALSE)

  if (!is.numeric(n_iter) || length(n_iter) != 1 || n_iter < 1) stop("`n_iter` must be >= 1.", call. = FALSE)
  if (!is.numeric(burn_in) || length(burn_in) != 1 || burn_in < 0) stop("`burn_in` must be >= 0.", call. = FALSE)
  n_iter <- as.integer(n_iter)
  burn_in <- as.integer(burn_in)
  n_total <- burn_in + n_iter
  if (n_total <= 1) stop("`burn_in + n_iter` must be > 1.", call. = FALSE)

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
    chain_seeds <- sample.int(1e9, n_chains)
  }

  if (!is.numeric(ci_level) || length(ci_level) != 1 || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be in (0,1).", call. = FALSE)
  }
  ci_levels_for_stars <- sort(unique(as.numeric(ci_levels_for_stars)))
  if (any(!is.finite(ci_levels_for_stars)) || any(ci_levels_for_stars <= 0) || any(ci_levels_for_stars >= 1)) {
    stop("`ci_levels_for_stars` must be in (0,1).", call. = FALSE)
  }

  y_full <- data[[outcome]]
  x_full <- data[[predictor]]

  keep_idx <- which(!is.na(y_full) & !is.na(x_full))
  if (!length(keep_idx)) stop("No complete cases for outcome and predictor.", call. = FALSE)

  y <- y_full[keep_idx]
  x <- x_full[keep_idx]

  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y) || is.character(y)) y <- as.integer(factor(y)) - 1L
  y <- as.numeric(y)
  if (!all(y %in% c(0, 1))) stop("Outcome must be binary (0/1) after missing handling.", call. = FALSE)
  if (length(unique(y)) < 2L) stop("Outcome must contain both 0 and 1 after missing handling.", call. = FALSE)

  if (!is.numeric(x)) {
    if (is.factor(x) || is.character(x) || is.logical(x)) {
      x <- as.factor(x)
      if (nlevels(x) != 2L) stop("Predictor has >2 levels; MEP_Univariate supports numeric or 2-level factors.", call. = FALSE)
      x <- as.integer(x) - 1L
    } else {
      x <- suppressWarnings(as.numeric(x))
      if (anyNA(x)) stop("Predictor could not be converted to numeric without NA.", call. = FALSE)
    }
  }

  x_for_disco <- as.numeric(scale(x))
  res_disco <- uni_separation(
    data = data.frame(y = y, x = x_for_disco),
    predictor = "x", outcome = "y",
    missing = "complete"
  )

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
        missing_info = list(rows_used = keep_idx, policy = "complete")
      ),
      message = "Predictor has zero variance after missing-data handling; model not fit.",
      rows_used = keep_idx
    ))
  }
  X_std_num <- as.numeric(X_std)

  severity <- res_disco$severity_score
  if (is.null(severity) || is.na(severity) || !is.finite(severity)) severity <- 0
  severity <- max(min(severity, 1), 0)

  y_bar <- mean(y)
  logit_clip <- function(p) stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))

  mu0 <- logit_clip(y_bar)
  mu1 <- 0
  log_sigma1 <- (1 - severity) * log(sigma1_hi) + severity * log(sigma1_lo)
  sigma1 <- exp(log_sigma1)
  kappa  <- kappa_min + severity * (kappa_max - kappa_min)

  mu <- c(mu0, mu1)
  Sigma <- diag(c(sigma0^2, sigma1^2))

  step_size0 <- step_hi * (1 - severity) + step_lo * severity
  log1pexp <- function(z) ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))

  run_MH_uni_one <- function(n_total, burn_in, init_beta, step_size, X_std, y, mu, Sigma, kappa,
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

    chain <- matrix(0, nrow = n_total, ncol = 2)
    colnames(chain) <- c("beta0", "beta1")
    chain[1, ] <- init_beta

    cur_lp <- log_post(chain[1, ])
    accept <- 0L
    ss <- step_size

    tune_pts <- if (burn_in >= tune_interval && tune_interval > 0) floor(burn_in / tune_interval) else 0L
    ss_trace <- if (tune_pts > 0L) numeric(tune_pts) else numeric(0)
    it_trace <- if (tune_pts > 0L) integer(tune_pts) else integer(0)
    ar_trace <- if (tune_pts > 0L) numeric(tune_pts) else numeric(0)
    idx_t <- 0L

    for (t in 2:n_total) {
      prop <- chain[t - 1, ] + ss * stats::rnorm(2)
      prop_lp <- log_post(prop)
      if (log(stats::runif(1)) < (prop_lp - cur_lp)) {
        chain[t, ] <- prop
        cur_lp <- prop_lp
        accept <- accept + 1L
      } else {
        chain[t, ] <- chain[t - 1, ]
      }

      if (t <= burn_in && tune_interval > 0 && (t %% tune_interval) == 0) {
        ar <- accept / t
        if (ar > tune_hi) ss <- ss * 1.10
        if (ar < tune_lo) ss <- ss * 0.90

        idx_t <- idx_t + 1L
        if (idx_t <= length(ss_trace)) {
          it_trace[idx_t] <- t
          ar_trace[idx_t] <- ar
          ss_trace[idx_t] <- ss
        }
      }
    }

    if (burn_in >= n_total) stop("`burn_in` must be < `burn_in + n_iter`.", call. = FALSE)

    post <- chain[(burn_in + 1):n_total, , drop = FALSE]
    list(
      post = post,
      acceptance_rate = accept / (n_total - 1),
      step_size_final = ss,
      burnin_step_trace = data.frame(iter = it_trace, acceptance = ar_trace, step_size = ss_trace)
    )
  }

  chains <- vector("list", n_chains)
  for (i in seq_len(n_chains)) {
    chains[[i]] <- run_MH_uni_one(
      n_total = n_total,
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

  if (combine_chains == "stack") {
    post_std <- do.call(rbind, lapply(chains, `[[`, "post"))
  } else {
    post_std <- chains[[1]]$post
  }

  diagnostics_multi <- list(
    rhat = c(beta0 = NA_real_, beta1 = NA_real_),
    rhat_max = NA_real_,
    ess = c(beta0 = NA_real_, beta1 = NA_real_),
    ess_min = NA_real_
  )

  diagnostics_single <- list(
    ess = rep(NA_real_, ncol(post_std)),
    geweke_z = rep(NA_real_, ncol(post_std)),
    ess_min = NA_real_,
    geweke_max_abs = NA_real_,
    converged = NA
  )

  if (requireNamespace("coda", quietly = TRUE)) {
    if (n_chains >= 2) {
      mlist <- coda::mcmc.list(lapply(chains, function(z) coda::mcmc(z$post)))

      gd <- coda::gelman.diag(mlist, autoburnin = FALSE, multivariate = FALSE)$psrf
      rhat <- as.numeric(gd[, "Point est."])
      names(rhat) <- colnames(chains[[1]]$post)

      diagnostics_multi$rhat <- rhat
      diagnostics_multi$rhat_max <- suppressWarnings(max(rhat, na.rm = TRUE))

      ess_m <- coda::effectiveSize(mlist)
      diagnostics_multi$ess <- as.numeric(ess_m)
      names(diagnostics_multi$ess) <- names(rhat)
      diagnostics_multi$ess_min <- suppressWarnings(min(diagnostics_multi$ess, na.rm = TRUE))
    } else {
      m <- coda::mcmc(post_std)
      ess <- as.numeric(coda::effectiveSize(m))
      gz  <- as.numeric(coda::geweke.diag(m)$z)

      diagnostics_single$ess <- ess
      diagnostics_single$geweke_z <- gz
      diagnostics_single$ess_min <- suppressWarnings(min(ess, na.rm = TRUE))
      diagnostics_single$geweke_max_abs <- suppressWarnings(max(abs(gz), na.rm = TRUE))

      ess_ok <- is.finite(diagnostics_single$ess_min) && diagnostics_single$ess_min >= ess_threshold
      gz_ok  <- is.finite(diagnostics_single$geweke_max_abs) && diagnostics_single$geweke_max_abs <= geweke_z_threshold
      diagnostics_single$converged <- isTRUE(ess_ok && gz_ok)
    }
  }

  # Back-transform (kept for transformed beta1 reporting and optional draws)
  transform_chain_to_original <- function(chain_std, x_mean, x_sd) {
    b0_std <- chain_std[, 1]
    b1_std <- chain_std[, 2]
    beta1_logit <- b1_std / x_sd
    beta0_orig  <- b0_std - b1_std * (x_mean / x_sd)
    logit_sd <- pi / sqrt(3)
    beta1_SAS  <- beta1_logit * logit_sd
    beta1_Long <- beta1_logit * (logit_sd + 1)
    cbind(beta0_orig = beta0_orig, beta1_logit = beta1_logit, beta1_SAS = beta1_SAS, beta1_Long = beta1_Long)
  }
  post_orig <- transform_chain_to_original(post_std, x_mean, x_sd)

  ci_mat_for_level <- function(M, lvl) {
    qlo <- (1 - lvl) / 2
    qhi <- 1 - qlo
    t(apply(M, 2, stats::quantile, probs = c(qlo, qhi), na.rm = TRUE))
  }

  summarize_draws <- function(draws_mat, ci_level_main, star_levels) {
    pm <- colMeans(draws_mat)
    sdv <- apply(draws_mat, 2, stats::sd)
    ci_main <- ci_mat_for_level(draws_mat, ci_level_main)

    star_levels <- sort(unique(star_levels))
    star <- rep("", ncol(draws_mat))
    for (j in seq_along(star_levels)) {
      ci_j <- ci_mat_for_level(draws_mat, star_levels[j])
      ok_j <- (ci_j[, 1] > 0) | (ci_j[, 2] < 0)
      star[ok_j] <- strrep("*", j)
    }

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

  # Posterior summaries: report only beta1 (and only the chosen transformed beta1 if requested)
  beta1_std <- matrix(post_std[, "beta1"], ncol = 1)
  colnames(beta1_std) <- "beta1"

  tab_b1 <- summarize_draws(beta1_std, ci_level, ci_levels_for_stars)
  posterior <- data.frame(
    Param = "beta1",
    tab_b1,
    row.names = NULL,
    check.names = FALSE
  )

  if (transform_beta != "none") {
    key <- if (transform_beta == "logit") "beta1_logit" else if (transform_beta == "SAS") "beta1_SAS" else "beta1_Long"
    b1_tr <- matrix(post_orig[, key], ncol = 1)
    colnames(b1_tr) <- key
    tab_tr <- summarize_draws(b1_tr, ci_level, ci_levels_for_stars)
    posterior <- rbind(
      posterior,
      data.frame(Param = key, tab_tr, row.names = NULL, check.names = FALSE)
    )
    rownames(posterior) <- NULL
  }

  comparators <- list(glm = NULL)
  if (isTRUE(compare)) {
    df_fit <- data.frame(y = y, X = X_std_num)
    glm_fit <- try(suppressWarnings(stats::glm(y ~ X, data = df_fit, family = stats::binomial())), silent = TRUE)
    if (!inherits(glm_fit, "try-error")) comparators$glm <- stats::coef(glm_fit)
  }

  mcmc_info <- list(
    acceptance_rate_by_chain = vapply(chains, function(z) z$acceptance_rate, numeric(1)),
    step_size_final_by_chain = vapply(chains, function(z) z$step_size_final, numeric(1)),
    burnin_step_trace_by_chain = lapply(chains, `[[`, "burnin_step_trace"),
    n_iter = n_iter,
    burn_in = burn_in,
    n_total = n_total,
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
      missing_info = list(rows_used = keep_idx, policy = "complete")
    ),
    prior = list(mu = mu, Sigma = Sigma, kappa = kappa, sigma0 = sigma0, sigma1 = sigma1),
    mcmc = mcmc_info,
    posterior = posterior,
    comparators = comparators,
    rows_used = keep_idx
  )

  if (n_chains >= 2) {
    out$diagnostics_multi <- diagnostics_multi
  } else {
    out$diagnostics_single <- diagnostics_single
  }

  if (isTRUE(return_draws)) {
    chain_std_by_chain <- lapply(chains, function(z) z$post)
    chain_orig_by_chain <- lapply(chains, function(z) transform_chain_to_original(z$post, x_mean, x_sd))

    out$draws <- list(
      chain_std = post_std,          # combined (stacked or first chain)
      chain_orig = post_orig,        # combined (stacked or first chain)
      chain_std_by_chain = chain_std_by_chain,    # list of length n_chains
      chain_orig_by_chain = chain_orig_by_chain   # list of length n_chains
    )
  }

  out
}
