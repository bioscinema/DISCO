#' Severity-Adaptive MEP for Pure Latent (multi-predictor) Logistic
#'
#' Fits a logistic regression model with a Multivariate Exponential Power (MEP) prior
#' using a random-walk Metropolis-Hastings (RW-MH) sampler and performs a small grid
#' search over prior settings \eqn{(\mu, \Sigma, \kappa)}. Predictors are encoded
#' (via \code{model.matrix}) and then z-scored internally for fitting using a safe scaler.
#' Summaries and credible intervals are reported on the working (logit, standardized-design)
#' scale and also back-transformed to the original encoded predictor scale.
#'
#' This function assumes the user inputs complete data. If any missing values are found
#' in \code{y} or \code{X}, the function stops with an error.
#'
#' @section What this function does:
#' \itemize{
#'   \item Validates that \code{y} and \code{X} contain no missing values and have compatible sizes.
#'   \item Encodes factors in \code{X} with \code{model.matrix(~ ., data = X)} using treatment contrasts
#'         (baseline is the first level), drops the intercept, and fits on the numeric encoded design.
#'   \item Standardizes encoded predictors with a safe scaler (sets sd = 1 when sd is 0 or non-finite).
#'   \item Runs one RW-MH chain for each grid setting and computes acceptance rate, posterior summaries,
#'         and a posterior predictive match statistic.
#'   \item Selects one grid setting using: (i) an acceptance-rate window (or closest to a target),
#'         (ii) closeness to a GLM coefficient-ratio reference (when available), and (iii) posterior
#'         predictive agreement.
#'   \item Reruns the selected grid point using \code{n_chains} chains and returns final summaries.
#' }
#'
#' @section Factor handling and column names:
#' \itemize{
#'   \item Numeric predictors remain a single encoded column with their original name.
#'   \item Factor predictors are expanded by \code{model.matrix()} using treatment contrasts with the first
#'         level as baseline. For a 2-level factor \code{X3} with levels \code{A} and \code{B} (baseline \code{A}),
#'         the encoded column is \code{X3B}, representing \code{B} versus \code{A}.
#'   \item To change the baseline, relevel before calling:
#' \preformatted{X$X3 <- stats::relevel(X$X3, ref = "B")}
#' }
#'
#' @param y Binary outcome (0/1, logical, or 2-level factor/character). Must contain both classes.
#' @param X Matrix or data.frame of predictors (no intercept). May include factors. Must contain no missing values.
#'
#' @param burn_in Integer; number of burn-in iterations per chain (discarded). Default \code{1000}.
#' @param n_iter Integer; number of post-burn-in iterations (posterior draws) per chain. Default \code{9000}.
#'
#' @param init_beta Initial value(s) for the MH chain. Either a scalar (recycled) or a numeric vector of length
#'   \eqn{p} (intercept plus encoded slopes). Default \code{0.01}.
#' @param step_size Proposal standard deviation for RW-MH. Default \code{0.40}.
#'
#' @param mu_vals Numeric vector; each value is repeated to length \eqn{p} to form \eqn{\mu} in the prior grid.
#' @param sigma0_intercept Prior scale for the intercept entry of \eqn{\Sigma} (logit scale). Default \code{10}.
#' @param sigma_global_multipliers Numeric vector of global multipliers applied to all slope prior scales
#'   (intercept held at \code{sigma0_intercept}). Default \code{c(0.1, 0.5, 1, 2, 5, 10)}.
#'
#' @param kappa_mode How to set \eqn{\kappa}. \code{"auto"} uses \code{kappa_vals} as a grid. \code{"fixed"}
#'   uses \code{kappa_fixed} for all runs. Default \code{"auto"}.
#' @param kappa_fixed Single positive value used when \code{kappa_mode="fixed"}. Default \code{1}.
#' @param kappa_vals Numeric vector of positive \eqn{\kappa} values used when \code{kappa_mode="auto"}.
#'   Default \code{c(0.5, 1, 2)}.
#'
#' @param accept_window Numeric length-2 vector giving the acceptable MH acceptance-rate window.
#'   Default \code{c(0.30, 0.40)}.
#' @param accept_target Numeric; target acceptance used if no grid point falls in \code{accept_window}.
#'   Default \code{0.35}.
#'
#' @param ci_level Credible interval level in \eqn{(0,1)}. Default \code{0.95}.
#' @param ppc_threshold Posterior predictive match threshold. Default \code{0.80}.
#'
#' @param tune_threshold_hi,tune_threshold_lo Burn-in acceptance thresholds for multiplicative step-size tuning.
#'   Increase step size if acceptance \code{> hi}; decrease if \code{< lo}. Defaults \code{0.45} and \code{0.20}.
#' @param tune_interval Iterations between tuning checks during burn-in. Default \code{500}.
#' @param verbose Logical; print brief progress messages. Default \code{FALSE}.
#'
#' @param n_chains Integer; number of MH chains to rerun for the selected best grid point. Default \code{1}.
#' @param chain_seeds Optional integer vector of length \code{n_chains} giving per-chain RNG seeds for
#'   the best-point reruns. If \code{NULL}, random seeds are generated.
#' @param combine_chains How to combine the best-point chains for final summaries. \code{"stack"} binds post-burn
#'   draws across chains; \code{"none"} uses only the first chain. Default \code{"stack"}.
#'
#' @param return_draws Logical; if \code{TRUE}, return post-burn draws for the selected best grid point
#'   (a matrix if one chain, else a list of matrices). Default \code{FALSE}.
#'
#' @param ess_threshold Minimum effective sample size required (across parameters) to declare convergence in the
#'   single-chain case when \pkg{coda} is available. Default \code{150}.
#' @param geweke_z_threshold Maximum allowed absolute Geweke z-score (across parameters) to declare convergence in
#'   the single-chain case when \pkg{coda} is available. Default \code{2}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{best_settings}: list with chosen \code{mu} (string), \code{Sigma_diag} (string), chosen \code{kappa},
#'         , \code{kappa_mode}, \code{acceptance_rate}, and \code{prop_matched}.
#'   \item \code{posterior_means}: posterior means (length \eqn{p}) for the selected run (working scale).
#'   \item \code{standardized_coefs_back}: data.frame per encoded column with means and CIs for standardized slopes
#'         and back-transformed effects (\code{b_A_original}, \code{b_SAS_original}, \code{b_Long_original}).
#'   \item \code{scaled_summary}: data.frame including Intercept with \code{Mean}, \code{SD}, \code{CI_low}, \code{CI_high},
#'         plus \code{Sig} and \code{Star}.
#'   \item \code{diagnostics_single}: returned only when \code{n_chains == 1} and \pkg{coda} is available.
#'   \item \code{diagnostics_multiple}: returned only when \code{n_chains >= 2} and \pkg{coda} is available.
#'   \item \code{burnin_step_trace_best}: list of burn-in tuning checkpoints per best-point chain.
#'   \item \code{step_size_final_best}: numeric vector of final tuned step sizes per best-point chain.
#'   \item \code{draws}: optional draws for the best-point reruns (if \code{return_draws = TRUE}).
#' }
#'
#' @importFrom stats quantile rbinom glm binomial coef sd plogis median
#' @export
#'
#' @examples
#' \donttest{
#' y <- c(0,0,0,0, 1,1,1,1)
#' X <- data.frame(
#'   X1 = c(-1.86, -0.81,  1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
#'   X2 = c( 0.52,  -0.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09)
#' )
#'
#' ## Single Chain
#' fit_single <- MEP_latent(
#'   y, X,
#'   n_chains = 1,
#'   chain_seeds = 9
#' )
#' fit_single$scaled_summary
#' fit_single$diagnostics_single
#'
#' ## Multiple Chains
#' fit_multi <- MEP_latent(
#'   y, X,
#'   n_chains = 4,
#'   chain_seeds = c(101, 102, 103, 104),
#'   combine_chains = "stack",
#'   return_draws = TRUE
#' )
#' fit_multi$scaled_summary
#' fit_multi$diagnostics_multiple
#' }
MEP_latent <- function(
    y, X,
    burn_in = 1000,
    n_iter = 9000,
    init_beta = 0.01,
    step_size = 0.4,
    mu_vals = seq(-1, 1, by = 0.1),
    sigma0_intercept = 10,
    sigma_global_multipliers = c(0.1, 0.5, 1, 2, 5, 10),
    kappa_mode = c("auto","fixed"),
    kappa_vals = c(0.5, 1, 2),
    kappa_fixed = 1,
    accept_window = c(0.3, 0.4),
    accept_target = 0.35,
    ci_level = 0.95,
    ppc_threshold = 0.80,
    tune_threshold_hi = 0.45,
    tune_threshold_lo = 0.20,
    tune_interval = 500,
    verbose = FALSE,
    n_chains = 1,
    chain_seeds = NULL,
    combine_chains = c("stack","none"),
    return_draws = FALSE,
    ess_threshold = 150,
    geweke_z_threshold = 2
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  combine_chains <- base::match.arg(as.character(combine_chains), choices = c("stack","none"))
  kappa_mode <- base::match.arg(as.character(kappa_mode), choices = c("auto","fixed"))

  if (!is.numeric(burn_in) || length(burn_in) != 1L || burn_in < 0) stop("`burn_in` must be >= 0.", call. = FALSE)
  if (!is.numeric(n_iter) || length(n_iter) != 1L || n_iter < 1) stop("`n_iter` must be >= 1 (post-burn draws).", call. = FALSE)
  burn_in <- as.integer(burn_in)
  n_iter <- as.integer(n_iter)
  n_total <- burn_in + n_iter
  if (n_total <= 1L) stop("`burn_in + n_iter` must be > 1.", call. = FALSE)

  if (!is.numeric(n_chains) || length(n_chains) != 1L || n_chains < 1) {
    stop("`n_chains` must be a positive integer.", call. = FALSE)
  }
  n_chains <- as.integer(n_chains)

  if (!is.null(chain_seeds)) {
    if (length(chain_seeds) != n_chains) stop("`chain_seeds` must have length `n_chains`.", call. = FALSE)
    chain_seeds <- as.integer(chain_seeds)
  } else {
    chain_seeds <- sample.int(1e9, n_chains)
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

  qfun_mat <- function(M, lvl) {
    qlo <- (1 - lvl) / 2
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

    ci90 <- qfun_mat(post, 0.90)
    ci95 <- qfun_mat(post, 0.95)
    ci99 <- qfun_mat(post, 0.99)

    star <- rep("", ncol(post))
    star[(ci90[, 1] > 0) | (ci90[, 2] < 0)] <- "*"
    star[(ci95[, 1] > 0) | (ci95[, 2] < 0)] <- "**"
    star[(ci99[, 1] > 0) | (ci99[, 2] < 0)] <- "***"

    scaled_summary$Sig <- (scaled_summary$CI_low > 0) | (scaled_summary$CI_high < 0)
    scaled_summary$Star <- star

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

  compute_diagnostics_single <- function(post, ess_threshold, geweke_z_threshold) {
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

  # y -> {0,1}
  y_full <- y
  if (!is.numeric(y_full)) {
    if (is.logical(y_full)) y_full <- as.integer(y_full)
    else if (is.factor(y_full) || is.character(y_full)) y_full <- as.integer(factor(y_full)) - 1L
    else stop("`y` must be numeric 0/1, logical, or 2-level factor/character.", call. = FALSE)
  }
  y_full <- as.numeric(y_full)
  if (anyNA(y_full)) stop("Missing values in `y` are not allowed. Please provide complete data.", call. = FALSE)
  if (!all(y_full %in% c(0, 1))) stop("`y` must be binary (0/1) after coercion.", call. = FALSE)
  if (length(unique(y_full)) < 2L) stop("`y` must contain both 0 and 1.", call. = FALSE)

  X_raw <- as.data.frame(X, stringsAsFactors = FALSE)
  if (nrow(X_raw) != length(y_full)) stop("Rows of X must match length of y.", call. = FALSE)
  if (ncol(X_raw) < 1L) stop("X must have at least one predictor.", call. = FALSE)
  if (anyNA(X_raw)) stop("Missing values in `X` are not allowed. Please provide complete data.", call. = FALSE)

  # encode factors, drop intercept
  mm <- stats::model.matrix(~ ., data = X_raw)
  X_enc <- mm[, -1, drop = FALSE]
  X_mat <- as.matrix(X_enc)
  p_enc <- ncol(X_mat)
  p_all <- 1L + p_enc

  # GLM ratio reference on standardized encoded design
  S_ref <- safe_scale(X_mat)
  X_ref <- S_ref$Xstd
  df_ref <- data.frame(y = y_full, X_ref)
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

  # grids
  mu_grid <- lapply(mu_vals, function(m) rep(m, p_all))
  build_sigma <- function(gm) diag(c(sigma0_intercept, rep(gm, p_all - 1L)), nrow = p_all, ncol = p_all)
  Sigma_list <- lapply(sigma_global_multipliers, build_sigma)

  if (kappa_mode == "fixed") {
    if (!is.numeric(kappa_fixed) || length(kappa_fixed) != 1L || !is.finite(kappa_fixed) || kappa_fixed <= 0) {
      stop("When kappa_mode = 'fixed', `kappa_fixed` must be a single positive finite number.", call. = FALSE)
    }
    kappa_grid <- as.numeric(kappa_fixed)
  } else {
    if (!is.numeric(kappa_vals) || length(kappa_vals) < 1L || any(!is.finite(kappa_vals)) || any(kappa_vals <= 0)) {
      stop("When kappa_mode = 'auto', `kappa_vals` must be a numeric vector of positive finite values.", call. = FALSE)
    }
    kappa_grid <- kappa_vals
  }

  run_chain_one <- function(n_total, burn_in, init_beta, step_size,
                            X_orig, y, mu, Sigma, kappa,
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

    if (length(init_beta) == 1L) init_vec <- rep(as.numeric(init_beta), p_all)
    else if (length(init_beta) == p_all) init_vec <- as.numeric(init_beta)
    else stop("`init_beta` must be a scalar or a numeric vector of length (1 + p_enc).", call. = FALSE)

    log_post <- function(beta) {
      eta <- as.vector(Xw %*% beta)
      ll  <- sum(-y * log1pexp(-eta) - (1 - y) * log1pexp(eta))
      diff <- beta - mu
      ll - 0.5 * (qform(diff)^kappa)
    }

    chain <- matrix(0, nrow = n_total, ncol = p_all)
    colnames(chain) <- colnames(Xw)
    chain[1, ] <- init_vec
    cur_lp <- log_post(chain[1, ])
    accept <- 0L
    ss <- step_size

    tune_pts <- if (burn_in >= tune_interval && tune_interval > 0) floor(burn_in / tune_interval) else 0L
    ss_trace <- if (tune_pts > 0L) numeric(tune_pts) else numeric(0)
    it_trace <- if (tune_pts > 0L) integer(tune_pts) else integer(0)
    ar_trace <- if (tune_pts > 0L) numeric(tune_pts) else numeric(0)
    idx_t <- 0L

    if (isTRUE(verbose)) message("RW-MH: total iters ", n_total, "; burn-in ", burn_in)

    for (it in 2:n_total) {
      prop <- chain[it - 1, ] + ss * stats::rnorm(p_all)
      prop_lp <- log_post(prop)

      if (log(stats::runif(1)) < (prop_lp - cur_lp)) {
        chain[it, ] <- prop
        cur_lp <- prop_lp
        accept <- accept + 1L
      } else {
        chain[it, ] <- chain[it - 1, ]
      }

      if (it <= burn_in && tune_interval > 0 && (it %% tune_interval) == 0) {
        ar <- accept / it
        if (ar > tune_threshold_hi) ss <- ss * 1.10
        if (ar < tune_threshold_lo) ss <- ss * 0.90

        idx_t <- idx_t + 1L
        if (idx_t <= length(ss_trace)) {
          it_trace[idx_t] <- it
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
      Xw = Xw,
      step_size_final = ss,
      burnin_step_trace = data.frame(
        iter = it_trace,
        acceptance = ar_trace,
        step_size = ss_trace,
        row.names = NULL
      )
    )
  }

  # grid search: one chain per grid point
  results_df <- data.frame(
    grid_id = integer(),
    mu = character(),
    sigma_diag = character(),
    kappa = numeric(),
    acceptance_rate = numeric(),
    prop_matched = numeric(),
    posterior_ratio_scaled = character(),
    step_size_final = numeric(),
    stringsAsFactors = FALSE
  )
  runs <- list()
  gid <- 1L

  for (mu in mu_grid) {
    mu_str <- paste(round(mu, 6), collapse = ", ")
    for (Sigma in Sigma_list) {
      sigma_diag_str <- paste(round(diag(Sigma), 6), collapse = ", ")
      for (kappa in kappa_grid) {
        ch <- run_chain_one(
          n_total = n_total,
          burn_in = burn_in,
          init_beta = init_beta,
          step_size = step_size,
          X_orig = X_mat,
          y = y_full,
          mu = mu,
          Sigma = Sigma,
          kappa = kappa,
          tune_threshold_hi = tune_threshold_hi,
          tune_threshold_lo = tune_threshold_lo,
          tune_interval = tune_interval,
          verbose = verbose,
          chain_seed = NULL
        )

        sum_one <- summarize_post(
          post = ch$post,
          X_orig = X_mat,
          ci_level = ci_level,
          ppc_threshold = ppc_threshold,
          y = y_full,
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
            step_size_final = ch$step_size_final,
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
          posterior_ratio_scaled = ratio_str,
          step_size_final = ch$step_size_final
        )

        gid <- gid + 1L
      }
    }
  }

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

  if (any(is.finite(cand$ratio_term))) cand <- cand[order(cand$ratio_term, -cand$prop_matched), , drop = FALSE]
  else cand <- cand[order(-cand$prop_matched), , drop = FALSE]

  best_id <- cand$grid_id[1]
  best_run <- runs[[best_id]]
  mu_best <- best_run$mu
  Sigma_best <- best_run$Sigma
  kappa_best <- best_run$kappa

  # rerun best point with n_chains chains
  best_chains <- vector("list", n_chains)
  for (ch_i in seq_len(n_chains)) {
    ch <- run_chain_one(
      n_total = n_total,
      burn_in = burn_in,
      init_beta = init_beta,
      step_size = step_size,
      X_orig = X_mat,
      y = y_full,
      mu = mu_best,
      Sigma = Sigma_best,
      kappa = kappa_best,
      tune_threshold_hi = tune_threshold_hi,
      tune_threshold_lo = tune_threshold_lo,
      tune_interval = tune_interval,
      verbose = verbose,
      chain_seed = chain_seeds[ch_i]
    )

    sum_one <- summarize_post(
      post = ch$post,
      X_orig = X_mat,
      ci_level = ci_level,
      ppc_threshold = ppc_threshold,
      y = y_full,
      Xw = ch$Xw
    )

    best_chains[[ch_i]] <- list(
      post = ch$post,
      acceptance_rate = ch$acceptance_rate,
      prop_matched = sum_one$prop_matched,
      step_size_final = ch$step_size_final,
      burnin_step_trace = ch$burnin_step_trace
    )
  }

  best_acceptance <- mean(vapply(best_chains, function(x) x$acceptance_rate, numeric(1)), na.rm = TRUE)
  best_prop_matched <- mean(vapply(best_chains, function(x) x$prop_matched, numeric(1)), na.rm = TRUE)

  # diagnostics for multiple chains only
  diagnostics_multiple <- list(
    rhat = rep(NA_real_, ncol(best_chains[[1]]$post)),
    rhat_max = NA_real_,
    ess = rep(NA_real_, ncol(best_chains[[1]]$post)),
    ess_min = NA_real_
  )
  if (requireNamespace("coda", quietly = TRUE) && n_chains >= 2) {
    mlist <- coda::mcmc.list(lapply(best_chains, function(x) coda::mcmc(x$post)))
    gd <- coda::gelman.diag(mlist, autoburnin = FALSE, multivariate = FALSE)$psrf
    diagnostics_multiple$rhat <- as.numeric(gd[, "Point est."])
    diagnostics_multiple$rhat_max <- suppressWarnings(max(diagnostics_multiple$rhat, na.rm = TRUE))
    ess_m <- coda::effectiveSize(mlist)
    diagnostics_multiple$ess <- as.numeric(ess_m)
    diagnostics_multiple$ess_min <- suppressWarnings(min(diagnostics_multiple$ess, na.rm = TRUE))
  }

  # combine post-burn draws
  if (combine_chains == "stack") post_all <- do.call(rbind, lapply(best_chains, `[[`, "post"))
  else post_all <- best_chains[[1]]$post

  # build Xw for summary (consistent columns)
  S_final <- safe_scale(X_mat)
  X_work_final <- S_final$Xstd
  Xw_final <- cbind(Intercept = 1, X_work_final)
  colnames(Xw_final) <- c("Intercept", colnames(X_mat))

  sum_final <- summarize_post(
    post = post_all,
    X_orig = X_mat,
    ci_level = ci_level,
    ppc_threshold = ppc_threshold,
    y = y_full,
    Xw = Xw_final
  )

  draws_out <- NULL
  if (isTRUE(return_draws)) {
    if (n_chains == 1L) draws_out <- best_chains[[1]]$post
    else draws_out <- lapply(best_chains, `[[`, "post")
  }

  out <- list(
    best_settings = list(
      mu = best_run$mu_str,
      Sigma_diag = best_run$sigma_diag_str,
      kappa = kappa_best,
      kappa_mode = kappa_mode,
      acceptance_rate = best_acceptance,
      prop_matched = best_prop_matched
    ),
    posterior_means = sum_final$posterior_means,
    standardized_coefs_back = sum_final$standardized_coefs_back,
    scaled_summary = sum_final$scaled_summary,
    burnin_step_trace_best = lapply(best_chains, `[[`, "burnin_step_trace"),
    step_size_final_best = vapply(best_chains, `[[`, numeric(1), "step_size_final"),
    draws = draws_out
  )

  if (requireNamespace("coda", quietly = TRUE)) {
    if (n_chains == 1L) {
      out$diagnostics_single <- compute_diagnostics_single(
        post = post_all,
        ess_threshold = ess_threshold,
        geweke_z_threshold = geweke_z_threshold
      )
    } else {
      out$diagnostics_multiple <- diagnostics_multiple
    }
  }

  out
}
