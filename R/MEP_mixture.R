#' Severity-Adaptive MEP for Mixture (multi-predictor) Logistic
#'
#' Fits a multi-predictor logistic model with a \strong{Multivariate Exponential Power (MEP)}
#' prior using a simple Random-Walk Metropolis–Hastings (RW–MH) sampler, where the \emph{slope
#' prior scales} are \strong{anchored by univariate DISCO severities}. A small grid over the
#' intercept prior mean, a global multiplier on slope scales, and the EP shape \eqn{\kappa}
#' is explored; one run is selected via acceptance window, posterior predictive agreement,
#' and—when available—closeness to GLM coefficient ratios relative to a reference predictor.
#'
#' \strong{Design encoding.} Predictors in \code{X} are expanded with
#' \code{model.matrix(~ ., data = X)} (intercept dropped for slopes). Numeric columns remain
#' one column each. Factor columns are expanded to treatment-contrast dummies (baseline is
#' the first level). The RW–MH is run on the \emph{standardized} encoded design (z-scored
#' columns), but standardized back-transforms (A, SAS, Long) are reported for each encoded
#' column using the \emph{unscaled} encoded design.
#'
#' @section What this function does:
#' \itemize{
#'   \item For each \emph{original} predictor in \code{X}, compute a univariate \strong{severity}
#'         using \code{DISCO::uni_separation()} on the modeling rows.
#'   \item Map severity \eqn{s \in [0,1]} to an anchor slope prior scale \eqn{\sigma_j(s)}
#'         and an anchor \eqn{\kappa_j(s)}; the intercept uses a wide prior.
#'   \item Build a small grid: intercept mean offsets, global multipliers on the \eqn{\sigma_j},
#'         and \eqn{\kappa} around the anchor average; run RW–MH for each grid point.
#'   \item Select one run using a score favoring acceptance in a target window, high posterior
#'         predictive agreement, and closeness of standardized
#'         slope ratios to GLM ratios (with a single \emph{reference} denominator).
#' }
#'
#' @section Factor handling & column names:
#' \itemize{
#'   \item \strong{Numeric predictors} appear as a single column with their original name
#'         (e.g., \code{X3}). No suffixes are added.
#'   \item \strong{Factor predictors} are expanded by \code{model.matrix()} using treatment
#'         contrasts with the \emph{first level as the baseline}. For a two-level factor
#'         \code{X3} with levels \code{A} and \code{B} (baseline = \code{A}), the encoded
#'         design includes a single dummy column \code{X3B}, which equals 1 when
#'         \code{X3 == "B"} and 0 when \code{X3 == "A"}. The coefficient for \code{X3B}
#'         is the log-odds difference \emph{B vs A} (controlling for other predictors).
#'   \item The output \code{posterior$effects$Predictor} uses these encoded names: numeric
#'         predictors as \code{"Xk"}; factor dummies as \code{"FactorLevel"} (e.g., \code{X3B}).
#'         Back-transforms (\code{b_A_original}, \code{b_SAS_original}, \code{b_Long_original})
#'         are computed per encoded column using that column’s standard deviation
#'         (for a 0/1 dummy with prevalence \eqn{p}, \eqn{\mathrm{sd}=\sqrt{p(1-p)}}).
#'   \item \strong{Change the baseline} beforehand to alter dummy labels:
#' \preformatted{X$X3 <- stats::relevel(X$X3, ref = "B")  # baseline becomes B; dummy shows as X3A}
#'   \item If you truly intend to \emph{treat a factor as numeric}, convert it yourself:
#' \preformatted{X$X3 <- as.numeric(X$X3)  # no dummy expansion; column remains 'X3'}
#' }
#'
#' @param y Numeric binary vector (0/1; logical or 2-level factor/character accepted; coerced to 0/1).
#' @param X Matrix or data.frame of predictors (no intercept). May include factors.
#'          Rows must align with \code{y}.
#' @param missing How to treat missing values in \emph{modeling}: only \code{"complete"} is
#'          supported; rows with any NA in \code{y} or \code{X} are dropped. Default \code{"complete"}.
#' @param severity_missing How to treat missing values when computing per-predictor DISCO
#'          severities: passed to \code{DISCO::uni_separation()} (\code{"complete"} or \code{"impute"}).
#'          Default \code{"complete"}.
#' @param impute_args Optional list of imputation settings used only when
#'          \code{severity_missing = "impute"} (see \code{DISCO::uni_separation}).
#'
#' @param n_iter_grid Integer; MH iterations per grid point (including burn-in). Default \code{10000}.
#' @param burn_in_grid Integer; burn-in iterations per grid point. Default \code{1000}.
#' @param step_size Proposal standard deviation for RW–MH. Default \code{0.40}.
#'
#' @param mu_intercept_offsets Numeric vector of offsets added to \eqn{\mathrm{logit}(\bar{y})}
#'          for the intercept prior mean grid. Default \code{seq(-1, 1, by = 0.2)}.
#' @param sigma0_intercept Prior \eqn{\sigma_0} (sd) for the intercept (logit scale). Default \code{10}.
#' @param sigma_global_multipliers Numeric vector of global multipliers applied to all
#'          slope prior scales (after severity anchoring). Default \code{c(0.1, 0.5, 1, 2, 5, 10)}.
#' @param sigma_hi Slope prior sd under mild separation (\eqn{s=0}). Default \code{5}.
#' @param sigma_lo Slope prior sd under severe separation (\eqn{s=1}). Default \code{0.15}.
#' @param kappa_min,kappa_max EP shape at \eqn{s=0} and \eqn{s=1}. Defaults \code{1}, \code{2.5}.
#' @param kappa_delta Offsets added around the anchor-average \eqn{\kappa} to form the grid.
#'          Default \code{seq(-0.5, 0.5, by = 0.2)} (truncated to the interval \eqn{[0.5, 3]}).
#'
#' @param accept_window Numeric length-2 vector; acceptable MH acceptance interval.
#'          Default \code{c(0.30, 0.40)}.
#' @param accept_target Scalar acceptance target used as a fallback (closest is preferred)
#'          when no grid point falls inside \code{accept_window}. Default \code{0.35}.
#'
#' @param ref Predictor \emph{name} or \emph{index} (in the original \code{X}) to serve as
#'          ratio denominator. If \code{NULL} (default), the predictor with the highest
#'          univariate severity is used. If the reference is a factor, the denominator
#'          column is the \emph{first} dummy generated for that factor.
#' @param compare Logical; if \code{TRUE} (default), fit GLM and Firth (\code{logistf})
#'          on the standardized encoded design as comparators.
#' @param seed For reproducibility, default 2025.
#' @param return_draws Logical; if \code{TRUE}, return the post-burn MH chain for the
#'          selected grid run.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{ref_predictor}: list with \code{index} (1-based in original \code{X}) and \code{name}.
#'   \item \code{severity}: data.frame with per-\emph{original} predictor severities used to anchor prior scales.
#'   \item \code{grid_summary}: data.frame summarizing all grid runs
#'         (\code{grid_id}, \code{mu}, \code{sigma_diag}, \code{kappa}, \code{acceptance_rate},
#'          \code{prop_matched}, \code{posterior_ratio_std}).
#'   \item \code{best}: list describing the selected grid point
#'         (\code{grid_id}, \code{mu}, \code{sigma_diag}, \code{kappa},
#'          \code{acceptance_rate}, \code{prop_matched}).
#'   \item \code{posterior}: list with
#'         \itemize{
#'           \item \code{means_std}: posterior means on the standardized encoded design
#'                 (length = intercept + encoded slopes),
#'           \item \code{effects}: data.frame with per encoded column:
#'                 \code{Predictor}, \code{Scaled} (standardized slope),
#'                 \code{b_A_original}, \code{b_SAS_original}, \code{b_Long_original}.
#'         }
#'   \item \code{ratios}: list of vectors on the encoded design:
#'         \code{GLM_beta}, \code{GLM_ratio}, \code{Firth_beta}, \code{Firth_ratio},
#'         \code{MEP_beta_std}, \code{MEP_ratio_std},
#'         and back-transform families \code{MEP_beta_b_A_original},
#'         \code{MEP_beta_b_SAS_original}, \code{MEP_beta_b_Long_original},
#'         plus their ratios \code{MEP_ratio_b_A_original},
#'         \code{MEP_ratio_b_SAS_original}, \code{MEP_ratio_b_Long_original}.
#'         Ratios divide by the encoded column corresponding to the reference predictor
#'         (first dummy if factor).
#'   \item \code{comparators}: list with raw \code{glm} and \code{firth} fits (if available).
#'   \item \code{draws}: matrix of MH draws after burn-in for the selected run
#'         (returned only when \code{return_draws = TRUE}).
#' }
#'
#' @details
#' \strong{Standardization & back-transforms.} The sampler runs on z-scored encoded columns.
#' Reported \code{Scaled} slopes are posterior means on this working scale. Back-transforms
#' use the encoded column SDs: \code{A} is per-SD effect on log-odds;
#' \code{SAS} is \code{A} multiplied by \eqn{\pi/\sqrt{3}};
#' \code{Long} is \code{A} multiplied by \eqn{\pi/\sqrt{3} + 1}.
#'
#'
#' \strong{Reference predictor & ratios.} The ratio denominator is chosen from the original
#' \code{X} columns (highest DISCO severity by default). If that predictor is a factor,
#' the denominator is the \emph{first} dummy column generated for that factor by
#' \code{model.matrix()}. GLM/Firth ratios are computed on the standardized encoded design.
#'
#' \strong{Missing data.} This version supports only \code{missing = "complete"} for modeling
#' (drop rows with NA in \code{y} or any column of \code{X}). For severity anchoring you may use
#' \code{severity_missing = "impute"} with \code{impute_args} which are passed to
#' \code{DISCO::uni_separation()}.
#'
#' @seealso \code{\link[DISCO]{uni_separation}}, \code{\link[logistf]{logistf}}
#'
#' @examples
#' \donttest{
#'
#' ## Example 2: toy with uni- and latent separation
#' ## X3 alone perfectly separates y (A for 0s, B for 1s);
#' ## X1 and X2 do not separate alone, but z = X1 - 5*X2 does (latent separation).
#' set.seed(1)
#' y <- c(0,0,0,0, 1,1,1,1)
#' X_toy <- data.frame(
#'   X1 = c(-1.86, -0.81,  1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
#'   X2 = c( 0.52,  -0.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09),
#'   X3 = factor(c(rep("A", 4), rep("B", 4)))
#' )
#' # Univariate diagnostics
#' d3 <- DISCO::uni_separation(data.frame(y=y, X3=X_toy$X3), "X3", "y", "complete")
#' d1 <- DISCO::uni_separation(data.frame(y=y, X1=X_toy$X1), "X1", "y", "complete")
#' d2 <- DISCO::uni_separation(data.frame(y=y, X2=X_toy$X2), "X2", "y", "complete")
#' d3$separation_type; d1$separation_type; d2$separation_type
#'
#' # Latent (manual): z = X1 - 5*X2 separates perfectly
#' z <- X_toy$X1 - 5 * X_toy$X2
#' dz <- DISCO::uni_separation(data.frame(y = y, z = z), "z", "y", "complete")
#' dz$separation_type  # "Perfect separation"
#'
#' # Latent (automatic): inclusion-minimal separating subsets of {X1, X2}
#' lat <- DISCO::latent_separation(
#'   y = y,
#'   X = X_toy[, c("X1", "X2")],
#'   find_minimal  = TRUE,
#'   mode          = "either",
#'   missing       = "complete",
#'   scale_X       = FALSE
#' )
#' names(lat$minimal_subsets)              # e.g., "X1_X2"
#' lat$minimal_subsets[[1]]$type           # "Perfect separation"
#' lat$minimal_subsets[[1]]$vars           # c("X1","X2")
#'
#' # Fit MEP_mixture (factor X3 will appear as a dummy 'X3B' = B vs A)
#' fit_toy <- MEP_mixture(y, X_toy, n_iter_grid = 4000, burn_in_grid = 1000, seed = 9)
#' fit_toy$ref_predictor
#' fit_toy$posterior$effects     # includes "X3B" (dummy for B vs A)
#'
#' # Change baseline level to rename the dummy (now "X3A" = A vs B)
#' X_toy2 <- X_toy
#' X_toy2$X3 <- stats::relevel(X_toy2$X3, ref = "B")
#' fit_toy2 <- MEP_mixture(y, X_toy2, n_iter_grid = 2000, burn_in_grid = 500, seed = 9)
#' fit_toy2$posterior$effects
#' }
#'
#' @importFrom DISCO uni_separation
#' @importFrom stats glm binomial coef plogis qlogis sd model.matrix terms
#' @importFrom logistf logistf
#' @export
MEP_mixture <- function(
    y, X,
    missing = c("complete"),
    severity_missing = c("complete","impute"),
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
    compare = TRUE,
    seed = 2025,
    return_draws = FALSE
) {
  missing <- match.arg(missing)
  severity_missing <- match.arg(severity_missing)
  if (!is.null(seed)) set.seed(as.integer(seed))

  # normalize y
  if (!is.numeric(y)) {
    if (is.logical(y)) y <- as.integer(y)
    else if (is.factor(y) || is.character(y)) y <- as.integer(factor(y)) - 1L
    else stop("`y` must be numeric 0/1, logical, or 2-level factor/character.", call. = FALSE)
  }

  # data prep
  X_df <- as.data.frame(X, stringsAsFactors = TRUE)
  stopifnot(nrow(X_df) == length(y))

  if (missing == "complete") {
    idx <- stats::complete.cases(data.frame(y = y, X_df, check.names = FALSE))
    y <- y[idx]; X_df <- X_df[idx, , drop = FALSE]
  } else {
    stop("Only 'complete' modeling missingness is supported in this version.", call. = FALSE)
  }
  if (ncol(X_df) < 1L) stop("X must have at least one predictor.", call. = FALSE)

  # encoded design (no intercept), and mapping to original terms
  mm <- stats::model.matrix(~ . , data = X_df)  # includes intercept
  assign_vec <- attr(mm, "assign")
  term_labs  <- attr(stats::terms(~ . , data = X_df), "term.labels")  # <-- fixed here
  X_mm <- mm[, -1, drop = FALSE]
  assign_mm <- assign_vec[-1]
  p_enc <- ncol(X_mm)
  p_terms <- length(term_labs)

  # univariate severities per original predictor
  sev_vec <- rep(0, p_terms)
  for (j in seq_len(p_terms)) {
    pred <- X_df[[term_labs[j]]]
    out <- try(DISCO::uni_separation(
      data = data.frame(y = y, x = pred),
      predictor = "x", outcome = "y",
      missing = severity_missing, impute_args = impute_args
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
  run_MH_sampler <- function(n_iter, init_beta, step_size, X_enc, y, mu, Sigma, kappa, burn_in = 1000) {
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
    se   <- apply(post, 2, stats::sd)

    # standardized back-transforms (A, SAS, Long)
    s_x <- apply(X_enc, 2, stats::sd)
    s_x[!is.finite(s_x) | s_x == 0] <- 1
    logit_sd <- pi / sqrt(3)
    b_scaled <- pm[-1]
    b_A      <- b_scaled / s_x
    b_SAS    <- b_A * logit_sd
    b_Long   <- b_A * (logit_sd + 1)

    effects <- data.frame(
      Predictor = colnames(X_enc),
      Scaled = b_scaled,
      b_A_original = b_A,
      b_SAS_original = b_SAS,
      b_Long_original = b_Long,
      row.names = NULL, check.names = FALSE
    )

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
      se = se,
      effects = effects,
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
          burn_in = burn_in_grid
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

  # GLM/Firth comparators
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
  safe_div <- function(num, den) ifelse(is.finite(den) & abs(den) > 1e-12, num / den, NA_real_)

  glm_beta <- glm_coef[-1]
  glm_ratio <- safe_div(glm_beta, glm_beta[ref_enc_cols[1]])

  firth_beta <- rep(NA_real_, length(glm_beta))
  firth_ratio <- rep(NA_real_, length(glm_beta))
  fit_firth <- NULL
  if (isTRUE(compare)) {
    fit_firth <- try(logistf::logistf(y ~ ., data = data.frame(y = y, X_scaled)), silent = TRUE)
    if (!inherits(fit_firth, "try-error")) {
      cf <- stats::coef(fit_firth)
      firth_beta <- cf[-1]
      firth_ratio <- safe_div(firth_beta, firth_beta[ref_enc_cols[1]])
    }
  }

  mep_beta_std <- pm[-1]
  mep_ratio_std <- safe_div(mep_beta_std, mep_beta_std[ref_enc_cols[1]])

  denom_row <- effects$Predictor == colnames(X_mm)[ref_enc_cols[1]]
  bA  <- effects$b_A_original
  bS  <- effects$b_SAS_original
  bL  <- effects$b_Long_original
  ratio_bA <- safe_div(bA, bA[denom_row])
  ratio_bS <- safe_div(bS, bS[denom_row])
  ratio_bL <- safe_div(bL, bL[denom_row])

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
    ratios = list(
      GLM_beta = glm_beta,
      GLM_ratio = glm_ratio,
      Firth_beta = firth_beta,
      Firth_ratio = firth_ratio,
      MEP_beta_std = mep_beta_std,
      MEP_ratio_std = mep_ratio_std,
      MEP_beta_b_A_original = bA,
      MEP_beta_b_SAS_original = bS,
      MEP_beta_b_Long_original = bL,
      MEP_ratio_b_A_original = ratio_bA,
      MEP_ratio_b_SAS_original = ratio_bS,
      MEP_ratio_b_Long_original = ratio_bL
    ),
    comparators = list(glm = glm_fit, firth = fit_firth),
    draws = if (isTRUE(return_draws)) best_run$chain else NULL
  )
}
