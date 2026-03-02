#' Severity-Adaptive MEP for Mixture (multi-predictor) Logistic
#'
#' Fits a multi-predictor logistic model with a Multivariate Exponential Power (MEP)
#' prior using a Random-Walk Metropolis-Hastings (RW-MH) sampler, where slope prior scales
#' are anchored by univariate DISCO severities. A small grid over intercept prior mean,
#' a global multiplier on slope scales, and the EP shape kappa is explored; one run is
#' selected via acceptance window, posterior predictive agreement, and (when available)
#' closeness to GLM coefficient ratios relative to a reference predictor.
#'
#' Design encoding. Predictors in X are expanded with model.matrix(~ ., data = X)
#' (intercept dropped for slopes). Numeric columns remain one column each. Factor columns
#' are expanded to treatment-contrast dummies (baseline is the first level). The RW-MH is
#' run on the standardized encoded design (z-scored columns). Optionally, coefficients can
#' be back-transformed to the original (unscaled) encoded columns via
#' transform_back = "logit"|"SAS"|"Long".
#'
#' @param y Numeric binary vector (0/1; logical or 2-level factor/character accepted; coerced to 0/1).
#' @param X Matrix or data.frame of predictors (no intercept). May include factors.
#'          Rows must align with y.
#' @param missing One of "complete" or "impute"; applied once to both modeling and severity diagnostics.
#'          Default "complete".
#' @param impute_args Optional list controlling simple external imputation when missing = "impute".
#'          Supported: numeric_method = "median"|"mean" (default "median"),
#'          factor_method = "mode" (default "mode").
#'
#' @param n_iter Integer; MH iterations per grid point (including burn-in). Default 10000.
#' @param burn_in Integer; burn-in iterations per grid point. Default 1000.
#' @param init_beta Initial value(s) for the MH chain. Default 0.01.
#' @param step_size Proposal standard deviation for RW-MH. Default 0.40.
#'
#' @param mu_intercept_offsets Numeric vector of offsets added to logit(mean(y)) for intercept prior mean grid.
#'          Default seq(-1, 1, by = 0.2).
#' @param sigma0_intercept Prior sd for the intercept (logit scale). Default 10.
#' @param sigma_global_multipliers Numeric vector of global multipliers applied to all slope prior scales.
#'          Default c(0.1, 0.5, 1, 2, 5, 10).
#' @param sigma_hi Slope prior sd under mild separation (s=0). Default 5.
#' @param sigma_lo Slope prior sd under severe separation (s=1). Default 0.15.
#' @param kappa_min,kappa_max EP shape at s=0 and s=1. Defaults 1 and 2.5.
#' @param kappa_delta Offsets around anchor-average kappa to form the grid.
#'          Default seq(-0.5, 0.5, by = 0.2), truncated to \eqn{[0.5, 3]}.
#'
#' @param accept_window Numeric length-2 vector; acceptable MH acceptance interval. Default c(0.30, 0.40).
#' @param accept_target Scalar acceptance target used if no grid point is inside accept_window. Default 0.35.
#'
#' @param ref Predictor name or index (in original X) to serve as ratio denominator. If NULL, uses max severity.
#' @param transform_back One of "none","logit","SAS","Long". Default "none".
#' @param ci_level Credible interval level in (0,1). Default 0.95.
#' @param seed Optional integer RNG seed.
#' @param return_draws Logical; if TRUE return post-burn draws for selected run.
#' @param ess_threshold Minimum ESS required (requires coda). Default 150.
#' @param geweke_z_threshold Maximum abs Geweke z allowed (requires coda). Default 2.
#' @param n_chains_best Integer; number of MH chains to rerun for the selected best grid point. Default 1.
#' @param chain_seeds_best Optional integer vector of length n_chains_best for best-point reruns.
#' @param combine_chains One of "stack" or "none". Default "stack".
#'
#' @param tune_threshold_hi,tune_threshold_lo Burn-in acceptance thresholds for multiplicative tuning
#'        (increase if > hi; decrease if < lo). Defaults 0.45 and 0.20.
#' @param tune_interval Iterations between tuning checks during burn-in. Default 1000.
#'
#' @return See function body for returned list components.
#'
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
    combine_chains = c("stack","none"),
    tune_threshold_hi = 0.45,
    tune_threshold_lo = 0.20,
    tune_interval = 1000
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
          if (!identical(fac_method, "mode")) {
            warning("Only factor_method = 'mode' is supported here; using mode.")
          }
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
        v2 <- suppressWarnings(as.numeric(v))
        if (anyNA(v2)) v2[is.na(v2)] <- 0
        v <- v2
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

  # normalize y to {0,1}
  if (!is.numeric(y)) {
    if (is.logical(y)) y <- as.integer(y)
    else if (is.factor(y) || is.character(y)) y <- as.integer(factor(y)) - 1L
    else stop("`y` must be numeric 0/1, logical, or 2-level factor/character.", call. = FALSE)
  }
  y <- as.numeric(y)

  # data prep
  X_df <- as.data.frame(X, stringsAsFactors = TRUE)
  if (nrow(X_df) != length(y)) stop("Rows of X must match length of y.", call. = FALSE)
  if (ncol(X_df) < 1L) stop("X must have at least one predictor.", call. = FALSE)

  rows_used <- NULL
  if (missing == "complete") {
    idx <- stats::complete.cases(data.frame(y = y, X_df, check.names = FALSE))
    if (!any(idx)) stop("No complete rows for y and X.", call. = FALSE)
    rows_used <- which(idx)
    y <- y[idx]
    X_df <- X_df[idx, , drop = FALSE]
  } else {
    keep <- !is.na(y)
    if (!any(keep)) stop("No rows with observed y.", call. = FALSE)
    rows_used <- which(keep)
    y <- y[keep]
    X_df <- X_df[keep, , drop = FALSE]
    X_df <- impute_frame(X_df, impute_args)
  }

  if (length(y) < 2L || length(unique(y)) < 2L) {
    stop("Outcome must contain both 0 and 1 after missing handling.", call. = FALSE)
  }

  # encoded design
  mm <- stats::model.matrix(~ . , data = X_df)
  assign_vec <- attr(mm, "assign")
  term_labs  <- attr(stats::terms(~ . , data = X_df), "term.labels")
  X_mm <- mm[, -1, drop = FALSE]
  assign_mm <- assign_vec[-1]
  p_enc <- ncol(X_mm)
  p_terms <- length(term_labs)

  # severities per original term
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
  ref_pos_enc <- 1 + ref_enc_cols[1]  # position in full beta (Intercept + slopes)

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
    tune_threshold_hi = 0.45, tune_threshold_lo = 0.20, tune_interval = 1000,
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
    colnames(chain) <- c("Intercept", colnames(X_enc))

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

    # record burn-in tuning checkpoints
    tune_pts <- if (burn_in >= tune_interval && tune_interval > 0) floor(burn_in / tune_interval) else 0L
    ss_trace <- if (tune_pts > 0L) numeric(tune_pts) else numeric(0)
    it_trace <- if (tune_pts > 0L) integer(tune_pts) else integer(0)
    ar_trace <- if (tune_pts > 0L) numeric(tune_pts) else numeric(0)
    idx_t <- 0L

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

      if (t <= burn_in && tune_interval > 0 && (t %% tune_interval) == 0) {
        ar <- acc / t
        if (ar > tune_threshold_hi) cur_step <- cur_step * 1.10
        if (ar < tune_threshold_lo) cur_step <- cur_step * 0.90

        idx_t <- idx_t + 1L
        if (idx_t <= length(ss_trace)) {
          it_trace[idx_t] <- t
          ar_trace[idx_t] <- ar
          ss_trace[idx_t] <- cur_step
        }
      }
    }

    post <- chain[(burn_in + 1):n_iter, , drop = FALSE]

    # posterior predictive check
    n_rep <- nrow(post)
    n_obs <- nrow(X)
    y_rep <- matrix(0L, nrow = n_rep, ncol = n_obs)
    for (i in seq_len(n_rep)) {
      pr_i <- 1 / (1 + exp(-(X %*% post[i, ])))
      y_rep[i, ] <- stats::rbinom(n_obs, 1, pr_i)
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
      diagnostics = sum_obj$diagnostics,
      step_size_final = cur_step,
      burnin_step_trace = data.frame(
        iter = it_trace,
        acceptance = ar_trace,
        step_size = ss_trace,
        row.names = NULL
      )
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
    step_size_final = numeric(),
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
          tune_threshold_hi = tune_threshold_hi,
          tune_threshold_lo = tune_threshold_lo,
          tune_interval = tune_interval,
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
            step_size_final = res$step_size_final,
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

  # GLM ratios on standardized encoded design
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
      tune_threshold_hi = tune_threshold_hi,
      tune_threshold_lo = tune_threshold_lo,
      tune_interval = tune_interval,
      chain_seed = chain_seeds_best[ch]
    )
  }

  prop_matched_best <- mean(vapply(best_chains, function(x) x$prop_matched, numeric(1)), na.rm = TRUE)
  acc_best <- mean(vapply(best_chains, function(x) x$acceptance_rate, numeric(1)), na.rm = TRUE)

  # multi-chain diagnostics
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
    rows_used = rows_used,
    missing_info = list(policy = missing, imputed = identical(missing, "impute")),
    burnin_step_trace_best = lapply(best_chains, `[[`, "burnin_step_trace"),
    step_size_final_best = vapply(best_chains, `[[`, numeric(1), "step_size_final"),
    draws = draws_out
  )
}
