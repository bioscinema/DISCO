# =============================================================================
# mep.R  —  MEP-regularized logistic regression via MH, with credible intervals
# =============================================================================

#' Metropolis–Hastings sampler for MEP-regularized logistic regression
#'
#' Run a simple random-walk Metropolis–Hastings sampler for a logistic model
#' with a Multivariate Exponential Power (MEP) prior on the coefficient vector.
#' Predictors can be optionally z-scored before fitting. The function returns
#' posterior draws (after burn-in), summary stats in the working (logit/scaled)
#' space, back-transformed coefficients to the original predictor scale, and
#' Bayesian credible intervals for all reported effects.
#'
#' @param n_iter Integer. Total MCMC iterations (including burn-in).
#' @param init_beta Numeric vector of length \eqn{p} (intercept + predictors).
#'   Starting values for the chain.
#' @param step_size Numeric. Proposal step-size (RW–MH standard deviation).
#' @param X_orig Matrix or data.frame of predictors (no intercept), \eqn{n \times (p-1)}.
#' @param y Binary outcome of length \eqn{n} (0/1, logical, or 2-level factor/character).
#' @param mu Numeric vector of length \eqn{p}. Prior location for the MEP prior (includes intercept).
#' @param Sigma Positive-definite \eqn{p \times p} matrix. Prior scale for the MEP prior.
#' @param kappa Positive numeric. Shape of the MEP prior; \eqn{\kappa = 1} resembles Laplace-like
#'   tails, larger values approach Gaussian.
#' @param burn_in Integer. Number of initial draws discarded (default \code{1000}).
#' @param scale_X Logical. If \code{TRUE} (default), columns of \code{X_orig} are z-scored
#'   before fitting; credible intervals are still reported on both the scaled and original scales.
#' @param ci_level Numeric in \eqn{(0,1)}. Credible interval level (default \code{0.95}).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{posterior_chain}}{Matrix of draws after burn-in (\code{iter} \eqn{\times} \eqn{p}).}
#'   \item{\code{posterior_means}}{Posterior means (length \eqn{p}).}
#'   \item{\code{se_estimates}}{Posterior standard deviations (length \eqn{p}).}
#'   \item{\code{acceptance_rate}}{Overall MH acceptance rate.}
#'   \item{\code{prop_matched}}{Posterior predictive check: proportion of observations
#'     whose replicated outcomes match the observed more than 80\% of the time.}
#'   \item{\code{scaled_summary}}{Data frame (including Intercept) with columns
#'     \code{Mean}, \code{SD}, \code{CI_low}, \code{CI_high} in the working (logit/possibly scaled) space.}
#'   \item{\code{standardized_coefs_back}}{Data frame (per predictor) with means and CIs for:
#'     \code{Scaled} (working space), and back-transformed effects on
#'     \code{b_A_original}, \code{b_SAS_original}, \code{b_Long_original}.}
#'   \item{\code{scaling_info}}{List with \code{scale_X} flag and centering/scaling vectors (if used).}
#'   \item{\code{ci_level}}{Credible interval level used.}
#' }
#'
#' @details
#' The MEP prior is proportional to \eqn{\exp\{-\tfrac12 ( ( \beta-\mu )^\top \Sigma^{-1} (\beta-\mu) )^\kappa\}}.
#' When \code{scale_X=TRUE}, the logistic design used in the sampler is z-scored; all
#' reported “original-scale” effects are computed by back-transforming the draws.
#'
#' Back-transformed effect scales:
#' \itemize{
#'   \item \strong{A}: per–SD-on-\code{X_orig} effect on the log-odds scale.
#'   \item \strong{SAS}: \code{A} multiplied by \eqn{\pi/\sqrt{3}}.
#'   \item \strong{Long}: \code{A} multiplied by \eqn{\pi/\sqrt{3} + 1}.
#'   }
#'
#' @seealso \code{\link{mep_grid_search}}, \code{\link{latent_separation}}, \code{\link{uni_separation}}
#' @importFrom stats quantile rbinom glm binomial coef sd plogis logLik
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- c(0,0,0,0,1,1,1,1)
#' X <- cbind(
#'   X1 = c(-1.86, -0.81,  1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
#'   X2 = c( 0.52,  1.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09)
#' )
#' p <- ncol(X) + 1
#' fit <- MEP_latent(
#'   n_iter = 10000, burn_in = 1000, init_beta = rep(0.01, p), step_size = 0.4,
#'   X_orig = X, y = y, mu = rep(0, p), Sigma = diag(p), kappa = 1,
#'   scale_X = TRUE, ci_level = 0.95
#' )
#' head(fit$scaled_summary)
#' head(fit$standardized_coefs_back)
#' }
#' @export
MEP_latent <- function(
    n_iter, init_beta, step_size, X_orig, y, mu, Sigma, kappa,
    burn_in = 1000, scale_X = TRUE, ci_level = 0.95
) {
  stopifnot(is.matrix(X_orig) || is.data.frame(X_orig))
  X_mat <- as.matrix(X_orig)
  if (is.null(colnames(X_mat))) colnames(X_mat) <- paste0("V", seq_len(ncol(X_mat)))

  # 0) Prepare design (optional z-score) + intercept
  X_work <- if (isTRUE(scale_X)) scale(X_mat) else X_mat
  X <- cbind(Intercept = 1, X_work)

  # 1) Log-posterior
  log_posterior <- function(beta, X, y, mu, Sigma, kappa) {
    Xb   <- as.vector(X %*% beta)
    pi_x <- stats::plogis(Xb)
    loglik <- sum(y * log(pi_x + 1e-12) + (1 - y) * log(1 - pi_x + 1e-12))

    diff <- beta - mu
    QuadForm <- as.numeric(t(diff) %*% solve(Sigma) %*% diff)
    logprior <- -0.5 * (QuadForm^kappa)

    p <- length(beta)
    logC <- log(kappa) + lgamma(p/2) - (p/2)*log(pi) - lgamma(p/(2*kappa)) - 0.5*log(det(Sigma))
    loglik + logprior + logC
  }

  # 2) Metropolis–Hastings
  MH_sampler <- function(n_iter, init_beta, step_size, X, y, mu, Sigma, kappa) {
    p <- length(init_beta)
    chain <- matrix(0, nrow = n_iter, ncol = p)
    chain[1, ] <- init_beta

    current_lp <- log_posterior(init_beta, X, y, mu, Sigma, kappa)
    accept <- 0L

    for (iter in 2:n_iter) {
      proposal <- chain[iter - 1, ] + step_size * rnorm(p)
      proposal_lp <- log_posterior(proposal, X, y, mu, Sigma, kappa)
      ratio <- proposal_lp - current_lp

      if (log(runif(1)) < ratio) {
        chain[iter, ] <- proposal
        current_lp <- proposal_lp
        accept <- accept + 1L
      } else {
        chain[iter, ] <- chain[iter - 1, ]
      }

      if (iter <= burn_in && iter %% 1000 == 0) {
        acc_rate <- accept / iter
        if (acc_rate > 0.6) step_size <- step_size * 1.1
        if (acc_rate < 0.2) step_size <- step_size * 0.9
      }
    }
    list(chain = chain, acceptance_rate = accept / (n_iter - 1))
  }

  # 3) Run sampler
  smp <- MH_sampler(n_iter, init_beta, step_size, X, y, mu, Sigma, kappa)
  chain <- smp$chain
  acc   <- smp$acceptance_rate
  chain_post <- chain[(burn_in + 1):n_iter, , drop = FALSE]

  # Basic summary in scaled/logit space
  post_means <- colMeans(chain_post)
  se_est     <- apply(chain_post, 2, stats::sd)

  qlo <- (1 - ci_level) / 2
  qhi <- 1 - qlo
  qfun <- function(M) t(apply(M, 2, stats::quantile, probs = c(qlo, qhi), na.rm = TRUE))

  ci_scaled <- qfun(chain_post)
  scaled_summary <- data.frame(
    Param   = c("Intercept", colnames(X_mat)),
    Mean    = post_means,
    SD      = se_est,
    CI_low  = ci_scaled[, 1],
    CI_high = ci_scaled[, 2],
    row.names = NULL, check.names = FALSE
  )

  # 4) Back-transform coefficients to original X scale (with CIs)
  back_transform <- function(chain_post, X_with_intercept, X_orig, y) {
    p <- ncol(X_with_intercept) - 1L
    predictors <- colnames(X_with_intercept)[-1]
    if (is.null(predictors)) predictors <- paste0("V", seq_len(p))

    # SDs on ORIGINAL X
    s_x <- apply(as.matrix(X_orig), 2, stats::sd)
    s_x[!is.finite(s_x) | s_x == 0] <- 1  # guard

    # samples in scaled/logit space (exclude intercept)
    samp_scaled <- chain_post[, -1, drop = FALSE]

    # A scale
    A_samp <- sweep(samp_scaled, 2, s_x, "/")

    # SAS & Long
    logit_sd <- pi / sqrt(3)
    SAS_samp  <- A_samp * logit_sd
    Long_samp <- A_samp * (logit_sd + 1)

    # M scale: need sample-wise R and s_lin
    df <- as.data.frame(X_with_intercept)
    df$y <- y
    ll_null <- as.numeric(stats::logLik(stats::glm(y ~ 1, data = df, family = stats::binomial())))

    n_rep <- nrow(chain_post)
    R2_vec   <- numeric(n_rep)
    s_linvec <- numeric(n_rep)

    for (i in seq_len(n_rep)) {
      lp  <- as.vector(X_with_intercept %*% chain_post[i, ])
      pr  <- stats::plogis(lp)
      ll  <- sum(y * log(pr + 1e-12) + (1 - y) * log(1 - pr + 1e-12))
      R2  <- 1 - (ll / ll_null)
      R2_vec[i] <- max(R2, 0)  # clamp at 0 to avoid sqrt(<0)
      s_linvec[i] <- stats::sd(lp)
    }
    R_vec <- sqrt(R2_vec)
    M_samp <- sweep(A_samp, 1, R_vec, "/")
    M_samp <- sweep(M_samp, 1, s_linvec, "*")

    summarise_mat <- function(M) {
      ci <- qfun(M)
      data.frame(
        Mean    = colMeans(M),
        CI_low  = ci[, 1],
        CI_high = ci[, 2],
        row.names = NULL
      )
    }

    out_scaled <- summarise_mat(samp_scaled)
    out_A      <- summarise_mat(A_samp)
    out_SAS    <- summarise_mat(SAS_samp)
    out_Long   <- summarise_mat(Long_samp)
    out_M      <- summarise_mat(M_samp)

    data.frame(
      Predictor = predictors,
      Scaled              = out_scaled$Mean,
      Scaled_CI_low       = out_scaled$CI_low,
      Scaled_CI_high      = out_scaled$CI_high,
      b_A_original        = out_A$Mean,
      b_A_CI_low          = out_A$CI_low,
      b_A_CI_high         = out_A$CI_high,
      b_SAS_original      = out_SAS$Mean,
      b_SAS_CI_low        = out_SAS$CI_low,
      b_SAS_CI_high       = out_SAS$CI_high,
      b_Long_original     = out_Long$Mean,
      b_Long_CI_low       = out_Long$CI_low,
      b_Long_CI_high      = out_Long$CI_high,
      row.names = NULL, check.names = FALSE
    )
  }

  std_back <- back_transform(
    chain_post = chain_post,
    X_with_intercept = X,
    X_orig = X_mat,
    y = y
  )

  # 5) Posterior predictive checks
  n_rep <- nrow(chain_post)
  n_obs <- nrow(X)
  y_rep <- matrix(0L, nrow = n_rep, ncol = n_obs)
  for (i in seq_len(n_rep)) {
    pr_i <- stats::plogis(X %*% chain_post[i, ])
    y_rep[i, ] <- rbinom(n_obs, size = 1, prob = pr_i)
  }
  p_match <- colMeans(sweep(y_rep, 2, y, FUN = "=="))
  prop_matched <- mean(p_match > 0.8)

  list(
    posterior_chain = chain_post,
    posterior_means = post_means,
    se_estimates    = se_est,
    acceptance_rate = acc,
    prop_matched    = prop_matched,
    scaled_summary  = scaled_summary,            # includes Intercept
    standardized_coefs_back = std_back,          # includes means + CIs for A/SAS/Long/M
    scaling_info = list(scale_X = scale_X,
                        center = if (isTRUE(scale_X)) attr(X_work, "scaled:center") else NULL,
                        scale  = if (isTRUE(scale_X)) attr(X_work, "scaled:scale")  else NULL),
    ci_level = ci_level
  )
}

#' @keywords internal
#' @noRd
.parse_ratio <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  as.numeric(strsplit(x, ",\\s*")[[1]])
}

#' Grid search over MEP priors with automatic selection and CIs
#'
#' Try a grid of prior settings \eqn{(\mu, \Sigma, \kappa)} for the MEP-regularized
#' logistic model, fit each via \code{\link{MEP_latent}}, and select a single
#' run using an acceptance-rate window and (if available) closeness to a GLM
#' coefficient-ratio reference. Credible intervals are carried through from the
#' selected run.
#'
#' @param y Binary outcome (0/1, logical, or 2-level factor/character).
#' @param X Matrix or data.frame of predictors (no intercept).
#' @param n_iter,burn_in,init_beta,step_size Passed to \code{\link{MEP_latent}}.
#' @param mu_vals Numeric vector; each value is repeated to length \eqn{p}
#'   (intercept + predictors) to form \eqn{\mu}.
#' @param Sigma_list List of \eqn{p \times p} positive-definite matrices for the prior scale.
#' @param kappa_vals Numeric vector of positive \eqn{\kappa} values.
#' @param accept_window Numeric length-2 vector giving the acceptable MH acceptance-rate
#'   window (default \code{c(0.3, 0.4)}).
#' @param accept_target Numeric. Target acceptance used to pick the closest run if
#'   the window yields no candidates (default \code{0.35}).
#' @param scale_X Logical. If \code{TRUE} (default), predictors are z-scored before fitting.
#' @param ci_level Numeric in \eqn{(0,1)}. Credible interval level (default \code{0.95}).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{best_settings}}{List with chosen \code{mu}, \code{Sigma}, \code{kappa}.}
#'   \item{\code{best_acceptance}}{Acceptance rate of the selected run.}
#'   \item{\code{best_prop_matched}}{Posterior predictive “match” statistic of the selected run.}
#'   \item{\code{posterior_means}}{Posterior means for the selected run.}
#'   \item{\code{standardized_coefs_back}}{Data frame of back-transformed effects with CIs
#'         (A / SAS / Long scales).}
#'   \item{\code{scaled_summary}}{Data frame of working-space summaries with CIs (includes intercept).}
#'   \item{\code{Ref_ratio}}{Character string of the GLM reference coefficient ratio (if available).}
#'   \item{\code{results_table}}{Data frame summarizing all grid runs (acceptance, ratios, etc.).}
#'   \item{\code{ci_level}}{Credible interval level used.}
#' }
#'
#' @seealso \code{\link{MEP_latent}}, \code{\link{latent_separation}}, \code{\link{uni_separation}}
#' @importFrom stats glm binomial coef
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- c(0,0,0,0,1,1,1,1)
#' X <- cbind(
#'   X1 = c(-1.86, -0.81,  1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
#'   X2 = c( 0.52,  1.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09)
#' )
#' gs <- mep_grid_search(y, X, n_iter = 10000, burn_in = 1000, scale_X = TRUE)
#' gs$best_settings
#' head(gs$standardized_coefs_back)
#' }
#' @export
mep_grid_search <- function(
    y, X,
    n_iter = 10000, burn_in = 1000, init_beta = NULL, step_size = 0.4,
    mu_vals = seq(-1, 1, by = 0.1),
    Sigma_list = list(diag(ncol(X)+1), diag(ncol(X)+1)*0.1, diag(ncol(X)+1)*0.5,
                      diag(ncol(X)+1)*2, diag(ncol(X)+1)*5),
    kappa_vals = c(0.5, 1, 2),
    accept_window = c(0.3, 0.4),
    accept_target = 0.35,
    scale_X = TRUE, ci_level = 0.95
) {
  X <- as.matrix(X)
  p <- ncol(X) + 1L
  if (is.null(init_beta)) init_beta <- rep(0.01, p)

  # reference GLM ratio (respect scaling)
  X_ref <- if (isTRUE(scale_X)) scale(X) else X
  glm_fit <- stats::glm(y ~ X_ref, family = stats::binomial())
  glm_coefs <- stats::coef(glm_fit)
  Ref_ratio <- if (length(glm_coefs) >= 3) {
    paste(round(glm_coefs[-c(1, 2)] / glm_coefs[2], 3), collapse = ", ")
  } else NA
  ref_ratio_vec <- .parse_ratio(Ref_ratio)

  mu_grid <- lapply(mu_vals, function(m) rep(m, p))

  results_df <- data.frame(
    mu = character(), sigma = character(), kappa = numeric(),
    acceptance_rate = numeric(), posterior_means = character(),
    prop_matched = numeric(), posterior_b_A_original_ratio = character(),
    stringsAsFactors = FALSE
  )
  grid_results_list <- list()
  grid_idx <- 1L

  for (mu in mu_grid) {
    mu_str <- paste(mu, collapse = ", ")
    for (Sigma in Sigma_list) {
      sigma_str <- paste(Sigma, collapse = ", ")
      for (kappa in kappa_vals) {

        run <- MEP_latent(
          n_iter = n_iter, init_beta = init_beta, step_size = step_size,
          X_orig = X, y = y, mu = mu, Sigma = Sigma, kappa = kappa,
          burn_in = burn_in, scale_X = scale_X, ci_level = ci_level
        )

        pm_str <- paste(round(run$posterior_means, 5), collapse = ", ")
        vals   <- run$posterior_means
        ratio_str <- if (length(vals) >= 3) {
          paste(round(vals[-c(1, 2)] / vals[2], 3), collapse = ", ")
        } else NA

        results_df <- rbind(
          results_df,
          data.frame(
            mu = mu_str, sigma = sigma_str, kappa = kappa,
            acceptance_rate = run$acceptance_rate,
            posterior_means = pm_str,
            prop_matched = run$prop_matched,
            posterior_b_A_original_ratio = ratio_str,
            stringsAsFactors = FALSE
          )
        )

        grid_results_list[[grid_idx]] <- list(
          mu_str = mu_str, sigma_str = sigma_str, kappa = kappa,
          acceptance_rate = run$acceptance_rate,
          posterior_means = run$posterior_means,
          prop_matched = run$prop_matched,
          standardized_coefs_back = run$standardized_coefs_back,  # has CIs
          scaled_summary = run$scaled_summary                      # has CIs incl. intercept
        )
        grid_idx <- grid_idx + 1L
      }
    }
  }

  # ---- selection criteria ----
  lo <- accept_window[1]; hi <- accept_window[2]
  filt <- results_df[results_df$acceptance_rate >= lo & results_df$acceptance_rate <= hi, ]
  if (nrow(filt) == 0) {
    closest_idx <- which.min(abs(results_df$acceptance_rate - accept_target))
    filt <- results_df[closest_idx, , drop = FALSE]
  }

  ref_vec <- ref_ratio_vec
  filt$Ratio_mean_difference <- sapply(filt$posterior_b_A_original_ratio, function(r) {
    post_vec <- .parse_ratio(r)
    if (all(is.na(c(post_vec, ref_vec)))) Inf else
      mean(abs(post_vec - ref_vec), na.rm = TRUE)
  })

  if (min(filt$prop_matched, na.rm = TRUE) < 0.90) {
    best_row <- filt[order(-filt$prop_matched, filt$Ratio_mean_difference), ][1, , drop = FALSE]
  } else {
    best_row <- filt[order(filt$Ratio_mean_difference), ][1, , drop = FALSE]
  }

  best <- NULL
  for (elem in grid_results_list) {
    if (elem$mu_str == best_row$mu && elem$sigma_str == best_row$sigma && elem$kappa == best_row$kappa) {
      best <- elem; break
    }
  }
  if (is.null(best)) stop("Best run not found (mep_grid_search).")

  list(
    best_settings = list(mu = best_row$mu, Sigma = best_row$sigma, kappa = best_row$kappa),
    best_acceptance = best$acceptance_rate,
    best_prop_matched = best$prop_matched,
    posterior_means = best$posterior_means,
    standardized_coefs_back = best$standardized_coefs_back,  # includes CIs
    scaled_summary = best$scaled_summary,                     # includes CIs incl. intercept
    Ref_ratio = Ref_ratio,
    results_table = results_df,
    ci_level = ci_level
  )
}
