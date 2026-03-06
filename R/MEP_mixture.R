#' Severity-Adaptive MEP for Mixture (multi-predictor) Logistic
#'
#' Fits a multi-predictor logistic model with a Multivariate Exponential Power (MEP)
#' prior using a Random-Walk Metropolis-Hastings (RW-MH) sampler. Slope prior scales
#' are anchored by univariate DISCO severities computed on the modeling data.
#' A small grid over intercept prior mean offsets, a global multiplier on slope scales,
#' and the EP shape kappa is explored. One run is selected using an acceptance-rate window,
#' posterior predictive agreement, and (when available) closeness to GLM coefficient ratios
#' relative to a reference predictor.
#'
#' Complete-case only. This function assumes inputs are complete; it stops if any NA
#' appears in y or X.
#'
#' Design encoding. Predictors in X are expanded with model.matrix(~ ., data = X)
#' (intercept dropped for slopes). Numeric columns remain one column each. Factor columns
#' are expanded to treatment-contrast dummies (baseline is the first level). RW-MH is run
#' on the standardized encoded design (z-scored columns). Back-transforms are reported
#' to the original (unscaled) encoded columns as b_A_original, b_SAS_original, b_Long_original.
#'
#' @param y Numeric binary vector (0/1; logical or 2-level factor/character accepted; coerced to 0/1).
#' @param X Matrix or data.frame of predictors (no intercept). May include factors. Rows must align with y.
#'
#' @param burn_in Integer; number of burn-in iterations per chain (discarded). Default 1000.
#' @param n_iter Integer; number of post-burn-in MCMC iterations (posterior draws) per chain. Default 9000.
#' @param init_beta Initial value(s) for the MH chain. Scalar (recycled) or numeric vector of length (1 + p_enc).
#'   Default 0.01.
#' @param step_size Proposal standard deviation for RW-MH. Default 0.40.
#'
#' @param mu_intercept_offsets Numeric vector of offsets added to logit(mean(y)) for intercept prior mean grid.
#'   Default seq(-1, 1, by = 0.2).
#' @param sigma0_intercept Prior sd for the intercept (logit scale). Default 10.
#' @param sigma_global_multipliers Numeric vector of global multipliers applied to all slope prior scales.
#'   Default c(0.1, 0.5, 1, 2, 5, 10).
#' @param sigma_hi Slope prior sd under mild separation (s=0). Default 5.
#' @param sigma_lo Slope prior sd under severe separation (s=1). Default 0.15.
#' @param kappa_min,kappa_max EP shape at s=0 and s=1. Defaults 1 and 2.5.
#' @param kappa_delta Offsets around anchor-average kappa to form the grid.
#'   Default seq(-0.5, 0.5, by = 0.2), truncated to \code{[0.5, 3]}.
#'
#' @param accept_window Numeric length-2 vector; acceptable MH acceptance interval. Default c(0.30, 0.40).
#' @param accept_target Scalar acceptance target used if no grid point is inside accept_window. Default 0.35.
#'
#' @param ref Predictor name or index (in original X) to serve as ratio denominator. If NULL, uses max severity.
#' @param transform_back One of "none","logit","SAS","Long". Default "none".
#' @param ci_level Credible interval level in (0,1). Default 0.95.
#' @param ppc_threshold Posterior predictive match threshold. Default \code{0.80}.
#'
#' @param return_draws Logical; if TRUE return post-burn draws for selected run.
#' @param n_chains Integer; number of MH chains to rerun for the selected best grid point. Default 1.
#' @param chain_seeds Optional integer vector of length n_chains for best-point reruns. If NULL, random seeds are generated.
#' @param combine_chains One of "stack" or "none". Default "stack".
#'
#' @param tune_threshold_hi,tune_threshold_lo Burn-in acceptance thresholds for multiplicative tuning
#'   (increase if > hi; decrease if < lo). Defaults 0.45 and 0.20.
#' @param tune_interval Iterations between tuning checks during burn-in. Default 1000.
#'
#' @param ess_threshold Minimum ESS required for single-chain diagnostics (requires coda). Default 150.
#' @param geweke_z_threshold Maximum abs Geweke z allowed for single-chain diagnostics (requires coda). Default 2.
#'
#' @return A list with:
#'   \item \code{best_settings}: list with chosen \code{mu} (string), \code{Sigma_diag} (string), chosen \code{kappa},
#'         , \code{kappa_mode}, \code{acceptance_rate}, and \code{prop_matched}.
#'   - posterior_means, scaled_summary, standardized_coefs_back
#'   - diagnostics_single (if n_chains==1) or diagnostics_multiple (if n_chains>=2) when coda is available
#'   - burnin_step_trace_best, step_size_final_best
#'   - draws (optional)
#'   - plus: ref_predictor, severity, grid_summary
#'
#' @export
#' @examples
#' \donttest{
#' y <- c(0,0,0,0, 1,1,1,1)
#' X <- data.frame(
#'   X1 = c(-1.86,-0.81, 1.32,-0.40, 0.91, 2.49, 0.34, 0.25),
#'   X2 = c( 0.52,-0.07, 0.60, 0.67,-1.39, 0.16,-1.40,-0.09),
#'   X3 = factor(c(rep("A",4), rep("B",4)))
#' )
#'
#' ## Single chain
#' fit_single <- MEP_mixture(
#'   y, X,
#'   n_chains = 1,
#'   chain_seeds = 9
#' )
#' fit_single$scaled_summary
#' fit_single$diagnostics_single
#'
#' ## Multiple chains
#' fit_multi <- MEP_mixture(
#'   y, X,
#'   n_chains = 4,
#'   chain_seeds = c(101, 102, 103, 104),
#'   combine_chains = "stack"
#' )
#' fit_multi$scaled_summary
#' fit_multi$diagnostics_multiple
#' }
MEP_mixture <- function(
    y, X,
    burn_in = 1000,
    n_iter = 9000,
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
    ppc_threshold = 0.80,
    return_draws = FALSE,
    n_chains = 1,
    chain_seeds = NULL,
    combine_chains = c("stack","none"),
    tune_threshold_hi = 0.45,
    tune_threshold_lo = 0.20,
    tune_interval = 1000,
    ess_threshold = 150,
    geweke_z_threshold = 2
) {

  transform_back <- base::match.arg(as.character(transform_back), choices = c("none","logit","SAS","Long"))
  combine_chains <- base::match.arg(as.character(combine_chains), choices = c("stack","none"))

  if (!is.numeric(burn_in) || length(burn_in) != 1 || burn_in < 0) stop("`burn_in` must be >= 0.", call. = FALSE)
  if (!is.numeric(n_iter) || length(n_iter) != 1 || n_iter < 1) stop("`n_iter` must be >= 1.", call. = FALSE)
  burn_in <- as.integer(burn_in)
  n_iter <- as.integer(n_iter)
  n_total <- burn_in + n_iter
  if (n_total <= 1) stop("`burn_in + n_iter` must be > 1.", call. = FALSE)

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
  # ci_levels_for_stars <- sort(unique(as.numeric(ci_levels_for_stars)))
  # if (any(!is.finite(ci_levels_for_stars)) || any(ci_levels_for_stars <= 0) || any(ci_levels_for_stars >= 1)) {
  #   stop("`ci_levels_for_stars` must be in (0,1).", call. = FALSE)
  # }

  # normalize y to {0,1}
  if (!is.numeric(y)) {
    if (is.logical(y)) y <- as.integer(y)
    else if (is.factor(y) || is.character(y)) y <- as.integer(factor(y)) - 1L
    else stop("`y` must be numeric 0/1, logical, or 2-level factor/character.", call. = FALSE)
  }
  y <- as.numeric(y)
  if (anyNA(y)) stop("`y` contains NA. This function assumes complete inputs.", call. = FALSE)
  if (!all(y %in% c(0, 1))) stop("`y` must be binary (0/1) after coercion.", call. = FALSE)
  if (length(unique(y)) < 2L) stop("`y` must contain both 0 and 1.", call. = FALSE)

  X_df <- as.data.frame(X, stringsAsFactors = TRUE)
  if (nrow(X_df) != length(y)) stop("Rows of X must match length of y.", call. = FALSE)
  if (ncol(X_df) < 1L) stop("X must have at least one predictor.", call. = FALSE)
  if (anyNA(X_df)) stop("`X` contains NA. This function assumes complete inputs.", call. = FALSE)

  # encoded design
  mm <- stats::model.matrix(~ . , data = X_df)
  assign_vec <- attr(mm, "assign")
  term_labs  <- attr(stats::terms(~ . , data = X_df), "term.labels")
  X_mm <- mm[, -1, drop = FALSE]
  assign_mm <- assign_vec[-1]
  p_enc <- ncol(X_mm)
  p_terms <- length(term_labs)
  p_all <- 1L + p_enc

  # severities per original term (anchoring)
  scale_if_num <- function(v) if (is.numeric(v)) as.numeric(scale(v)) else v
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
    sev_vec[j] <- max(0, min(1, s))
  }
  severity_df <- data.frame(Predictor = term_labs, Severity = sev_vec, stringsAsFactors = FALSE)

  # reference predictor selection
  if (is.null(ref)) {
    ref_idx_term <- which.max(sev_vec)
  } else if (is.character(ref)) {
    ref_idx_term <- match(ref, term_labs)
    if (is.na(ref_idx_term)) stop("`ref` not found in predictors.", call. = FALSE)
  } else {
    ref_idx_term <- as.integer(ref)
    if (ref_idx_term < 1 || ref_idx_term > p_terms) stop("`ref` index out of range.", call. = FALSE)
  }
  ref_name <- term_labs[ref_idx_term]
  ref_enc_cols <- which(assign_mm == ref_idx_term)
  if (!length(ref_enc_cols)) stop("Internal: no encoded columns for selected `ref` predictor.", call. = FALSE)
  ref_pos_enc <- 1 + ref_enc_cols[1]  # position in full beta

  map_uni_severity <- function(s) {
    sigma <- exp((1 - s) * log(sigma_hi) + s * log(sigma_lo))
    kappa_anchor <- kappa_min + s * (kappa_max - kappa_min)
    list(sigma = sigma, kappa_anchor = kappa_anchor)
  }

  anchors <- lapply(sev_vec, map_uni_severity)
  sigma_anchor_terms <- vapply(anchors, function(a) a$sigma, numeric(1))
  kappa_anchor_mean  <- mean(vapply(anchors, function(a) a$kappa_anchor, numeric(1)))
  sigma_anchor_enc <- sigma_anchor_terms[assign_mm]

  # grids
  mu0 <- stats::qlogis(pmin(pmax(mean(y), 1e-6), 1 - 1e-6))
  mu_grid <- lapply(mu_intercept_offsets, function(off) {
    v <- numeric(p_all)
    v[1] <- mu0 + off
    v
  })

  build_sigma <- function(global_mult) {
    d <- c(sigma0_intercept, pmax(1e-6, sigma_anchor_enc * global_mult))
    diag(d, nrow = p_all, ncol = p_all)
  }
  Sigma_list <- lapply(sigma_global_multipliers, build_sigma)
  kappa_grid <- pmax(0.5, pmin(3.0, kappa_anchor_mean + kappa_delta))

  # helpers
  safe_scale <- function(M) {
    cen <- suppressWarnings(colMeans(M))
    Xc  <- sweep(M, 2, cen, FUN = "-", check.margin = FALSE)
    sc  <- suppressWarnings(apply(M, 2, stats::sd))
    sc[!is.finite(sc) | sc == 0] <- 1
    Xs  <- sweep(Xc, 2, sc, FUN = "/", check.margin = FALSE)
    list(Xstd = Xs, center = cen, scale = sc)
  }

  log1pexp <- function(z) ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))

  qfun_mat <- function(M, lvl) {
    qlo <- (1 - lvl) / 2
    qhi <- 1 - qlo
    t(apply(M, 2, stats::quantile, probs = c(qlo, qhi), na.rm = TRUE))
  }

  summarize_post <- function(post, X_enc_mat, ci_level#, ci_levels_for_stars
                             ) {
    pm <- colMeans(post)
    sdv <- apply(post, 2, stats::sd)
    ci_main <- qfun_mat(post, ci_level)

    scaled_summary <- data.frame(
      Param = colnames(post),
      Mean = pm,
      SD = sdv,
      CI_low = ci_main[, 1],
      CI_high = ci_main[, 2],
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

    # back-transforms (encoded scale)
    s_x <- apply(X_enc_mat, 2, stats::sd)
    s_x[!is.finite(s_x) | s_x == 0] <- 1

    slopes_scaled <- post[, -1, drop = FALSE]
    A_samp <- sweep(slopes_scaled, 2, s_x, "/")
    logit_sd <- pi / sqrt(3)
    SAS_samp <- A_samp * logit_sd
    Long_samp <- A_samp * (logit_sd + 1)

    summarise_slopes <- function(M) {
      ci <- qfun_mat(M, ci_level)
      data.frame(Mean = colMeans(M), CI_low = ci[, 1], CI_high = ci[, 2])
    }

    out_scaled <- summarise_slopes(slopes_scaled)
    out_A <- summarise_slopes(A_samp)
    out_SAS <- summarise_slopes(SAS_samp)
    out_Long <- summarise_slopes(Long_samp)

    standardized_coefs_back <- data.frame(
      Predictor = colnames(X_enc_mat),
      Scaled = out_scaled$Mean,
      Scaled_CI_low = out_scaled$CI_low,
      Scaled_CI_high = out_scaled$CI_high,
      b_A_original = out_A$Mean,
      b_A_CI_low = out_A$CI_low,
      b_A_CI_high = out_A$CI_high,
      b_SAS_original = out_SAS$Mean,
      b_SAS_CI_low = out_SAS$CI_low,
      b_SAS_CI_high = out_SAS$CI_high,
      b_Long_original = out_Long$Mean,
      b_Long_CI_low = out_Long$CI_low,
      b_Long_CI_high = out_Long$CI_high,
      row.names = NULL,
      check.names = FALSE
    )

    list(
      posterior_means = pm,
      scaled_summary = scaled_summary,
      standardized_coefs_back = standardized_coefs_back
    )
  }

  run_chain_one <- function(mu, Sigma, kappa, chain_seed) {
    if (!is.null(chain_seed)) set.seed(as.integer(chain_seed))

    S <- safe_scale(as.matrix(X_mm))
    Xstd <- S$Xstd
    Xmat <- cbind(Intercept = 1, Xstd)
    colnames(Xmat) <- c("Intercept", colnames(X_mm))

    Sigma_inv <- solve(Sigma)

    log_post <- function(beta) {
      eta <- as.vector(Xmat %*% beta)
      ll <- sum(-y * log1pexp(-eta) - (1 - y) * log1pexp(eta))
      diff <- beta - mu
      qf <- as.numeric(t(diff) %*% Sigma_inv %*% diff)
      ll - 0.5 * (qf^kappa)
    }

    if (length(init_beta) == 1L) init_vec <- rep(as.numeric(init_beta), p_all)
    else if (length(init_beta) == p_all) init_vec <- as.numeric(init_beta)
    else stop("`init_beta` must be a scalar or a numeric vector of length (1 + p_enc).", call. = FALSE)

    chain <- matrix(0, nrow = n_total, ncol = p_all)
    colnames(chain) <- colnames(Xmat)
    chain[1, ] <- init_vec

    cur_lp <- log_post(chain[1, ])
    acc <- 0L
    ss <- step_size

    tune_pts <- if (burn_in >= tune_interval && tune_interval > 0) floor(burn_in / tune_interval) else 0L
    ss_trace <- if (tune_pts > 0L) numeric(tune_pts) else numeric(0)
    it_trace <- if (tune_pts > 0L) integer(tune_pts) else integer(0)
    ar_trace <- if (tune_pts > 0L) numeric(tune_pts) else numeric(0)
    idx_t <- 0L

    for (t in 2:n_total) {
      prop <- chain[t - 1, ] + ss * stats::rnorm(p_all)
      prop_lp <- log_post(prop)

      if (log(stats::runif(1)) < (prop_lp - cur_lp)) {
        chain[t, ] <- prop
        cur_lp <- prop_lp
        acc <- acc + 1L
      } else {
        chain[t, ] <- chain[t - 1, ]
      }

      if (t <= burn_in && tune_interval > 0 && (t %% tune_interval) == 0) {
        ar <- acc / t
        if (ar > tune_threshold_hi) ss <- ss * 1.10
        if (ar < tune_threshold_lo) ss <- ss * 0.90

        idx_t <- idx_t + 1L
        if (idx_t <= length(ss_trace)) {
          it_trace[idx_t] <- t
          ar_trace[idx_t] <- ar
          ss_trace[idx_t] <- ss
        }
      }
    }

    post <- chain[(burn_in + 1):n_total, , drop = FALSE]

    # posterior predictive match
    n_rep <- nrow(post)
    n_obs <- nrow(Xmat)
    y_rep <- matrix(0L, nrow = n_rep, ncol = n_obs)
    for (i in seq_len(n_rep)) {
      pr_i <- stats::plogis(Xmat %*% post[i, ])
      y_rep[i, ] <- stats::rbinom(n_obs, 1, pr_i)
    }
    p_match <- colMeans(sweep(y_rep, 2, y, `==`))
    prop_matched <- mean(p_match >= ppc_threshold, na.rm = TRUE)
    if (!is.finite(prop_matched)) prop_matched <- 0

    list(
      post = post,
      acceptance_rate = acc / (n_total - 1),
      prop_matched = prop_matched,
      step_size_final = ss,
      burnin_step_trace = data.frame(iter = it_trace, acceptance = ar_trace, step_size = ss_trace, row.names = NULL)
    )
  }

  # GLM ratio reference on standardized encoded design
  glm_fit <- try(suppressWarnings(stats::glm(
    y ~ .,
    data = data.frame(y = y, scale(X_mm)),
    family = stats::binomial()
  )), silent = TRUE)
  glm_ok <- !(inherits(glm_fit, "try-error"))
  glm_coef <- if (glm_ok) stats::coef(glm_fit) else rep(NA_real_, p_all)

  parse_ratio <- function(s) {
    if (is.null(s) || is.na(s) || !nzchar(s)) return(NA_real_)
    as.numeric(strsplit(s, ",\\s*")[[1]])
  }

  ref_ratio_str <- if (length(glm_coef) >= ref_pos_enc && is.finite(glm_coef[ref_pos_enc]) && abs(glm_coef[ref_pos_enc]) > 0) {
    num <- glm_coef[-c(1, ref_pos_enc)]
    if (!length(num)) NA_character_ else paste(round(num / glm_coef[ref_pos_enc], 3), collapse = ", ")
  } else {
    NA_character_
  }
  ref_ratio_vec <- parse_ratio(ref_ratio_str)

  # grid search: one chain per grid point
  grid_summary <- data.frame(
    grid_id = integer(),
    mu = character(),
    sigma_diag = character(),
    kappa = numeric(),
    acceptance_rate = numeric(),
    prop_matched = numeric(),
    posterior_ratio_std = character(),
    step_size_final = numeric(),
    stringsAsFactors = FALSE
  )
  runs <- list()
  gid <- 1L

  for (mu in mu_grid) {
    mu_str <- paste(round(mu, 6), collapse = ", ")
    for (Sigma in Sigma_list) {
      sigma_str <- paste(round(diag(Sigma), 6), collapse = ", ")
      for (kappa in kappa_grid) {

        ch <- run_chain_one(mu = mu, Sigma = Sigma, kappa = kappa, chain_seed = NULL)

        pm_tmp <- colMeans(ch$post)
        num_idx <- setdiff(2:(1 + p_enc), ref_pos_enc)
        ratio_std_str <- if (length(num_idx) &&
                             is.finite(pm_tmp[ref_pos_enc]) && abs(pm_tmp[ref_pos_enc]) > 1e-12) {
          paste(round(pm_tmp[num_idx] / pm_tmp[ref_pos_enc], 3), collapse = ", ")
        } else {
          NA_character_
        }

        grid_summary <- rbind(
          grid_summary,
          data.frame(
            grid_id = gid,
            mu = mu_str,
            sigma_diag = sigma_str,
            kappa = kappa,
            acceptance_rate = ch$acceptance_rate,
            prop_matched = ch$prop_matched,
            posterior_ratio_std = ratio_std_str,
            step_size_final = ch$step_size_final,
            stringsAsFactors = FALSE
          )
        )

        runs[[gid]] <- list(mu = mu, mu_str = mu_str, Sigma = Sigma, sigma_str = sigma_str, kappa = kappa, chain = ch)
        gid <- gid + 1L
      }
    }
  }

  # selection
  cand <- grid_summary[
    is.finite(grid_summary$acceptance_rate) &
      grid_summary$acceptance_rate >= accept_window[1] &
      grid_summary$acceptance_rate <= accept_window[2],
    , drop = FALSE
  ]
  if (nrow(cand) == 0) {
    idx <- which.min(abs(grid_summary$acceptance_rate - accept_target))
    cand <- grid_summary[idx, , drop = FALSE]
  }

  cand$ratio_term <- NA_real_
  if (!all(is.na(ref_ratio_vec))) {
    cand$ratio_term <- vapply(cand$posterior_ratio_std, function(r) {
      pv <- parse_ratio(r)
      if (all(is.na(pv))) return(NA_real_)
      if (length(pv) != length(ref_ratio_vec)) return(NA_real_)
      mean(abs(pv - ref_ratio_vec))
    }, numeric(1))
  }

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

  # rerun best grid point with n_chains
  best_chains <- vector("list", n_chains)
  for (i in seq_len(n_chains)) {
    best_chains[[i]] <- run_chain_one(mu = mu_best, Sigma = Sigma_best, kappa = kappa_best, chain_seed = chain_seeds[i])
  }

  best_acceptance <- mean(vapply(best_chains, function(x) x$acceptance_rate, numeric(1)), na.rm = TRUE)
  best_prop_matched <- mean(vapply(best_chains, function(x) x$prop_matched, numeric(1)), na.rm = TRUE)

  # combine draws
  if (combine_chains == "stack") {
    post_all <- do.call(rbind, lapply(best_chains, `[[`, "post"))
  } else {
    post_all <- best_chains[[1]]$post
  }

  sum_final <- summarize_post(
    post = post_all,
    X_enc_mat = as.matrix(X_mm),
    ci_level = ci_level
    # ,
    # ci_levels_for_stars = ci_levels_for_stars
  )

  # diagnostics: single or multiple (only when coda is available)
  diagnostics_single <- list(
    ess = rep(NA_real_, ncol(post_all)),
    geweke_z = rep(NA_real_, ncol(post_all)),
    ess_min = NA_real_,
    geweke_max_abs = NA_real_,
    converged = NA
  )
  diagnostics_multiple <- list(
    rhat = rep(NA_real_, ncol(best_chains[[1]]$post)),
    rhat_max = NA_real_,
    ess = rep(NA_real_, ncol(best_chains[[1]]$post)),
    ess_min = NA_real_
  )

  if (requireNamespace("coda", quietly = TRUE)) {
    if (n_chains >= 2) {
      mlist <- coda::mcmc.list(lapply(best_chains, function(x) coda::mcmc(x$post)))
      gd <- coda::gelman.diag(mlist, autoburnin = FALSE, multivariate = FALSE)$psrf
      diagnostics_multiple$rhat <- as.numeric(gd[, "Point est."])
      diagnostics_multiple$rhat_max <- suppressWarnings(max(diagnostics_multiple$rhat, na.rm = TRUE))
      ess_m <- coda::effectiveSize(mlist)
      diagnostics_multiple$ess <- as.numeric(ess_m)
      diagnostics_multiple$ess_min <- suppressWarnings(min(diagnostics_multiple$ess, na.rm = TRUE))
    } else {
      m <- coda::mcmc(post_all)
      ess <- as.numeric(coda::effectiveSize(m))
      gz <- as.numeric(coda::geweke.diag(m)$z)

      diagnostics_single$ess <- ess
      diagnostics_single$geweke_z <- gz
      diagnostics_single$ess_min <- suppressWarnings(min(ess, na.rm = TRUE))
      diagnostics_single$geweke_max_abs <- suppressWarnings(max(abs(gz), na.rm = TRUE))

      ess_ok <- is.finite(diagnostics_single$ess_min) && diagnostics_single$ess_min >= ess_threshold
      gz_ok <- is.finite(diagnostics_single$geweke_max_abs) && diagnostics_single$geweke_max_abs <= geweke_z_threshold
      diagnostics_single$converged <- isTRUE(ess_ok && gz_ok)
    }
  }

  draws_out <- NULL
  if (isTRUE(return_draws)) {
    if (n_chains == 1L) draws_out <- best_chains[[1]]$post
    else draws_out <- lapply(best_chains, `[[`, "post")
  }

  out <- list(
    ref_predictor = list(index = ref_idx_term, name = ref_name),
    severity = severity_df,
    grid_summary = grid_summary,
    best_settings = list(
      mu = best_run$mu_str,
      Sigma_diag = best_run$sigma_str,
      kappa = kappa_best,
      acceptance_rate = best_acceptance,
      prop_matched = best_prop_matched
    ),
    posterior_means = sum_final$posterior_means,
    standardized_coefs_back = sum_final$standardized_coefs_back,
    scaled_summary = sum_final$scaled_summary,
    burnin_step_trace_best = lapply(best_chains, `[[`, "burnin_step_trace"),
    step_size_final_best = vapply(best_chains, `[[`, numeric(1), "step_size_final")
  )

  if (n_chains >= 2) out$diagnostics_multiple <- diagnostics_multiple else out$diagnostics_single <- diagnostics_single
  if (isTRUE(return_draws)) out$draws <- draws_out

  out
}
