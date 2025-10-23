# DISCO
**DI**agnosis of **S**eparation and **C**orrection of **O**dds-ratio inflation in logistic regression

> Fast diagnostics for **perfect** and **quasi-complete** separation in binary outcomes, with clear severity scoring, per-subset missing-data handling, and pretty `gt` tables — plus severity-adaptive and MEP-regularized Bayesian estimators.

---

## Why separation matters
When predictors **perfectly** (or almost perfectly) split the outcome, logistic regression can produce **infinite** (or severely inflated) odds ratios and unstable inference. **DISCO** helps you **detect** these cases early and (optionally) **stabilize** estimates with Bayesian priors.

---

## Overview

### Separation Diagnosis
- `uni_separation()` — univariate detector using:
  - Hungarian matching to align clusters,
  - a **vectorized Rand index**,
  - a **non-negative boundary threshold** to guard against boundary ties,
  - a **continuous severity score** in \[0,1\].
- `latent_separation()` — multivariate detector using **LP-based** feasibility with options to:
  - search for **inclusion-minimal** separating subsets (pruning),
  - choose **perfect**, **quasi**, or **either** as “hit” criteria,
  - handle missingness **per subset** (complete-case or imputation).

### Estimation (Bayesian)
- `EP_univariable()` — **severity-adaptive univariate** logistic regression with an **Exponential Power (EP)** prior (RW–MH sampler).  
  **Default output:** standardized coefficients `beta0` and `beta1` (z-scored X working scale).  
  **Optional back-transform:** set `transform_beta1 = "logit" | "SAS" | "Long"` to add `beta0_orig` and a slope per 1 unit of the **original** predictor.
- `mep_function()` — **core multi-predictor** logistic sampler with a **Multivariate Exponential Power (MEP)** prior for user-specified \((\mu,\Sigma,\kappa)\). Returns working-scale summaries **with CIs** and back-transformed effects (**A / SAS / Long**, with CIs).
- `MEP_latent()` — **grid search** over \((\mu,\Sigma,\kappa)\) for multi-predictor models using `mep_function()`; selects one run using an **acceptance-rate window**, **posterior predictive agreement**, and (when available) **GLM coefficient-ratio** closeness. Carries through **CIs** from the selected run.
- `MEP_mixture()` — **severity-adaptive multi-predictor** logistic with **numeric + factor** predictors. Encodes factors via `model.matrix(~ ., data = X)`, anchors slope prior scales by **per-predictor DISCO severities**, runs a small grid (intercept mean offsets, global slope multipliers, κ), and selects one run via acceptance window + checks. Reports posterior **means** on the standardized scale and back-transformed **A / SAS / Long** **means** per encoded column; also provides GLM/Firth/MEP **ratios**.

Planned: odds-ratio deflation after separation detection.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("bioscinema/DISCO")
library(DISCO)
```

## Separation Diagnosis
### Quick Start

```r
library(DISCO)

#### I. Univariate: quick diagnostic #############
df <- data.frame(Y = c(0,1,0,1), X = c(-2, 2, -1, 1))
res <- uni_separation(df, predictor = "X", outcome = "Y", missing = "complete")
res$separation_type     # "Perfect separation" | "Quasi-complete separation" | "No separation problem"
res$severity_score      # in [0,1]
res$boundary_threshold  # data-driven non-negative threshold
res$missing_info        # method, params, rows_used, n_used

#### II. Latent (multivariate): LP-based check ####
y <- c(0,0,0,0,1,1,1,1)
X <- cbind(
  X1 = c(-1.86, -0.81, 1.32, -0.40, 0.91, 2.49, 0.34, 0.25),
  X2 = c( 0.52,  1.07, 0.60,  0.67,-1.39, 0.16,-1.40,-0.09)
)
lat <- latent_separation(y, X, missing = "complete")  # scale_X = FALSE by default
lat$type           # "Perfect separation" | "Quasi-complete separation" | "No separation problem"
lat$message        # human-readable summary
lat$removed        # for quasi: variables whose removal yields perfect (else NULL)
lat$missing_info   # method, params, rows_used, n_used

# Optional: enable z-score scaling of the encoded predictors
lat_scaled <- latent_separation(y, X, missing = "complete", scale_X = TRUE)

# Search for inclusion-minimal separating subsets (pruned search)
lat_min <- latent_separation(
  y, X,
  find_minimal   = TRUE,
  min_vars       = 2,              # smallest subset size to consider
  max_vars       = ncol(X),        # largest subset size (default is all)
  stop_at_first  = FALSE
)
names(lat_min$minimal_subsets)     # e.g., "X1", "X2", "X1_X2"
one <- lat_min$minimal_subsets[[1]]
one$type          # "Perfect separation" or "Quasi-complete separation"
one$vars          # variables in the subset
one$removed       # for quasi-only; variables whose removal yields perfect
one$missing_info  # rows used, etc.
```

### Missing-Data Handling

- `missing = "complete"`: drop rows with **NA** in the tested variables; the **outcome is never imputed**.
- `missing = "impute"`: **impute predictors only**; outcome NAs are **always dropped**.
- Every result returns `missing_info` with:
  - `method` (`"complete"` or `"impute"`)
  - `params` (imputation settings or `custom_fn`)
  - `n_used` (number of rows analyzed after the missing-data rule)
  - `rows_used` (exact **original** row indices included in the fit)
- In latent (subset) searches, **each subset** has its own `missing_info` (its own `rows_used`/`n_used`) since missingness can differ by variable combination.

#### Imputation Params (`missing_info$params`)
A compact summary of per-type rules when `method = "impute"` (e.g., `Numeric=Mean; Categorical=Missing; Logical=Mode`):

- **Numeric** (`numeric_method`)
  - `Mean` — replace `NA` with the column mean  
  - `Median` — replace `NA` with the column median  
  - `Constant=<c>` — replace `NA` with a supplied constant (e.g., `Constant=0`)  
  - `Custom Imputer` — a user-provided function was used

- **Categorical** (`categorical_method`)
  - `Missing` — turn `NA` into an explicit `"Missing"` level (adds a level if needed)  
  - `Mode` — replace `NA` with the most frequent level  
  - `Constant=<level>` — replace `NA` with a supplied level  
  - `Custom Imputer` — a user-provided function was used

- **Logical** (`logical_method`)
  - `Mode` — replace `NA` with the most frequent value (`TRUE`/`FALSE`)  
  - `Constant=TRUE/FALSE` — replace `NA` with the given logical  
  - `Custom Imputer` — a user-provided function was used

**Mode definition.** *Mode* is the value/level occurring **most often** among non-missing entries. If there’s a tie, it is broken **deterministically** (e.g., first level / first in `table()` order).

Example with missing values:

```r
set.seed(2025)
df_miss <- data.frame(
  Y    = c(0,0,0,0,1,1,1,1, 0,1, NA, 1),
  X1   = c(0.5, 1.0, NA, 2.0, 5.0, 6.0, 7.0, NA, 1.5, 8.0, 9.0, 6.5),
  X2   = c(10, 9, 8, NA, 6, 5, NA, 3, 2, 1, 0, NA),
  Race = factor(c("A","A","B", NA, "C","C","B","A","B", NA, "C","A")),
  L1   = c(TRUE, NA, FALSE, TRUE, TRUE, NA, FALSE, TRUE, FALSE, TRUE, TRUE, NA)
)

# Complete-case
res_cc <- uni_separation(df_miss, "X1", "Y", missing = "complete")
res_cc$missing_info

# Impute (defaults = numeric median, categorical mode, logical mode)
res_imp <- uni_separation(df_miss, "X1", "Y", missing = "impute")
res_imp$missing_info

# Treat NA as a level for categorical predictors
res_cat <- uni_separation(
  df_miss, "Race", "Y", missing = "impute",
  impute_args = list(categorical_method = "missing")
)
res_cat$missing_info$params
```

### Pretty tables (gt) — optional

```r
# install.packages("gt")
library(gt)

# One predictor
res_uni_cc <- uni_separation(df_miss, "X1", "Y", missing = "complete")
gt_uni_separation(res_uni_cc, title = "Univariate (X1 vs Y) — Complete-case")

# All predictors vs outcome (one-shot summary)
gt_uni_separation_all(df_miss, outcome = "Y", missing = "complete")
```

### Example output

![Univariate DISCO table](man/figures/readme-uni-gt-all.png)

This table summarizes univariate screen results of each predictor against the outcome (`Y`) using **complete-case** data. It flags whether any single predictor causes separation in a logistic model.

- **Predictor / Outcome**  
  The variable tested and the binary outcome.

- **Indices (block of columns)**
  - **Separation Index** — scaled 0–1; **1.000 = Perfect Separation**, values near 0 indicate no separation concern.
  - **Severity** — scaled 0–1 and shown with color; higher = more severe (red), lower = minimal/none (green).
  - **Boundary Threshold** — the univariate cut that separates classes at the boundary (reported for reference).
  - **Single-Tie Boundary** — **Yes** if separation hinges on a single boundary case; **No** otherwise.
  - **Tie Count** — number of tied boundary rows involved in the separating cut.

- **Missing-Data Handling (block of columns)**
  - **Missing Method** — *Complete* = complete-case analysis.
  - **Imputation Params** — parameters shown when imputation is used (here they’re placeholders because complete-case was used).
  - **N Used** — number of subjects actually analyzed after the missing-data rule.

- **Separation**  
  Text label summarizing the result:
  - **Perfect Separation** — predictor alone perfectly separates the outcome (MLE will diverge).
  - **No Problem** — no separation detected.

- **Rows Used (Original Indices)**  
  Original dataset row numbers included in the analysis for that predictor (not re-indexed after filtering).

- Example

  - **X1** has **Separation Index = 1.000**, **Severity = 1.000**, and is labeled **Perfect Separation** with **N Used = 9**.  
    → X1 alone perfectly separates `Y` on the complete-case subset (rows `1, 2, 4, 5, 6, 7, 9, 10, 12`). Standard logistic regression MLE will not be finite; consider remedies (e.g., Firth penalization, Bayesian priors, or data/model modifications).

  - **Race, L1, X2** show **Separation Index ≈ 0.57–0.61**, **Severity = 0.000**, and **No Problem**.  
    → These predictors do **not** induce separation in univariate fits on their respective complete-case samples (N Used = 8–9).

  - All rows used are reported as **original indices**. Since the method is *Complete*, the full list is shown for each predictor.

---
```r
# Latent: minimal subsets, complete-case per subset
res_lat_cc <- latent_separation(
  y = df$Y,
  X = df[, c("X1","X2","Race","L1")],
  find_minimal = TRUE,
  missing = "complete"
)
gt_latent_separation(res_lat_cc, title = "Latent Minimal Subsets — Complete-Case")
```
![Latent DISCO table](man/figures/readme-latent-gt-complete.png)

This table summarizes subsets of predictors that yield separation in a **complete-case analysis** (rows with any missing values are excluded). Each row corresponds to one predictor subset.

- **Subset / Variables / # Of Predictors**  
  The predictor combination evaluated and its size.

- **Missing-Data Handling / Imputation Params / N Used**  
  The approach to missing data (*Complete* means complete-case) and how many subjects were analyzed for that subset. Imputation parameters are shown when applicable.

- **Separation**  
  The separation status observed with the listed subjects as used:
  - **Perfect Separation** — the subset already perfectly separates the outcome.
  - **Quasi-Complete Separation** — separation is nearly perfect (one side has boundary cases).

- **Removed And Rest Reach Perfect**  
  A **minimal** set of original row indices that, if dropped, would make the **remaining used** rows achieve **Perfect Separation**.
  - Empty/blank here means no removal is needed (the subset already has Perfect Separation).
  - When multiple indices appear, **all** must be removed; removing a strict subset may not suffice.

- **Rows Used (Original Indices)**  
  The original dataset row numbers included in the complete-case fit for that subset. We always report **original indices** (not re-indexed after filtering).

- Example
  - **X1_X2, X2_L1, Race_L1, X1_Race, X2_Race**  
    - **Separation:** Perfect Separation  
    - **Removed And Rest Reach Perfect:** *(blank)* → already Perfect; no rows need removal.  
    - **N Used:** 6–7, with original indices listed in the last column.

  - **X1_L1**  
    - **Separation:** Quasi-Complete Separation  
    - **Removed And Rest Reach Perfect:** `5`  
    - **Interpretation:** The complete-case fit used rows `1, 4, 5, 7, 9, 10`. If you drop row **5**, the remaining rows (`1, 4, 7, 9, 10`) would yield **Perfect Separation** for subset `{X1, L1}`.
---

```r
# Latent minimal subsets (imputed)
res_lat_imp <- latent_separation(
  y = df_miss$Y,
  X = df_miss[, c("X1","X2","Race","L1")],
  find_minimal = TRUE,
  missing = "impute",
  impute_args = list(numeric_method = "mean", categorical_method = "missing", logical_method = "mode")
)
gt_latent_separation(res_lat_imp, title = "Latent Minimal Subsets — Imputed")
```
![Latent DISCO table](man/figures/readme-latent-gt.png)

This table summarizes subsets of predictors that yield separation in an **imputation analysis**. Each row corresponds to one predictor subset.

- **Subset / Variables / # Of Predictors**  
  The subset name, which variables it includes, and its size.

- **Missing-Data Handling / Imputation Params / N Used**  
  - **Missing Method** — *Impute* indicates missing values in predictors were filled according to the listed recipe.  
  - **Imputation Params** — e.g., `Numeric=Mean; Categorical=Missing; Logical=Mode` means numeric features use mean imputation, categorical features add a literal “Missing” level, and logicals use the mode.  
  - **N Used** — number of subjects used after applying the imputation rule (can be less than total if the **outcome** is missing for some rows, since outcomes are not imputed).

- **Separation**  
  The separation status **after imputation** using the rows counted in *N Used*:
  - **Perfect Separation** — the subset perfectly separates the outcome on the imputed dataset.
  - **Quasi-Complete Separation** — nearly perfect; boundary cases exist (not shown in this example).

- **Removed And Rest Reach Perfect**  
  A **minimal** set of original row indices that, *if removed*, would make the **remaining used** rows achieve **Perfect Separation**.  
  - Blank = already Perfect; no removal needed.  
  - If multiple indices appear, **all** are required (removing a strict subset may not suffice).

- **Rows Used (Original Indices)**  
  Original dataset row numbers included in the imputed fit for that subset.  
  - In **imputed** mode this is hidden by default (shown as “—”) to avoid long lists.  
  - Set `show_rows_used = TRUE` to preview indices. Indices are always reported in the **original** numbering (not re-indexed after filtering/imputation).

- Example
  - All listed pairs (e.g., **X1_X2**, **X1_Race**, **X2_Race**, **Race_L1**) achieve **Perfect Separation** under the specified imputation scheme with **N Used = 11**.
  - **Removed And Rest Reach Perfect** is blank for each row → no removals are needed because separation is already perfect.
  - **Practical note:** Separation that appears **after imputation** can reflect either genuine structure or an artifact of the imputation rule (e.g., adding a “Missing” level). Consider sensitivity checks (complete-case vs. imputed analyses), penalized likelihood (Firth), or Bayesian priors when fitting logistic models in the presence of separation.

---
---

## Severity-adaptive Univariate Bayes (EP_univariable)

`EP_univariable()` fits **intercept + one predictor** logistic regression with a **severity-adaptive EP prior** informed by `uni_separation()`. It uses a random-walk MH sampler with light auto-tuning.

**Defaults (match SLURM script):**
- `n_iter = 20000`, `burn_in = 5000`
- Proposal s.d. blend: `step = 0.30*(1-severity) + 0.12*severity`
- Prior scales: `sigma0 = 10`, `sigma1_hi = 5`, `sigma1_lo = 0.15`
- Shape blend: `kappa = 1 + severity*(2.5 - 1)`
- Missing-data handling propagated from `uni_separation()`

**Output (by default):** standardized coefficients `beta0`, `beta1`  
*(z-scored X working scale; good for across-predictor comparability)*

**Optional back-transform:** choose one:
- `transform_beta1 = "logit"` → add `beta0_orig`, `beta1_logit` (per 1 unit of original X)  
- `transform_beta1 = "SAS"`   → add `beta0_orig`, `beta1_SAS = beta1_logit * π/√3`  
- `transform_beta1 = "Long"`  → add `beta0_orig`, `beta1_Long = beta1_logit * (π/√3 + 1)`

### Quick start

```r
set.seed(1)
n  <- 60
x  <- rnorm(n, mean = 0.3, sd = 1)
eta <- -0.2 + 1.0 * x
y  <- rbinom(n, size = 1, prob = plogis(eta))
df <- data.frame(y = y, x = x)

# (1) Default: STANDARDIZED beta0, beta1
fit_std <- EP_univariable(df, predictor = "x", outcome = "y",
                          n_iter = 6000, burn_in = 2000, seed = 42)
fit_std$posterior

# (2) Back-transform slope to ORIGINAL-x units on LOGIT scale
fit_logit <- EP_univariable(df, "x", "y",
                            transform_beta1 = "logit",
                            n_iter = 6000, burn_in = 2000, seed = 42)
fit_logit$posterior

# (3) Alternative effect scales on ORIGINAL-x units
fit_sas  <- EP_univariable(df, "x", "y", transform_beta1 = "SAS",
                           n_iter = 6000, burn_in = 2000, seed = 42)
fit_long <- EP_univariable(df, "x", "y", transform_beta1 = "Long",
                           n_iter = 6000, burn_in = 2000, seed = 42)
```

**Interpretation**
- `beta1` (default) — change in log-odds per **1 SD** increase in X (standardized scale).
- `beta1_logit` — change in log-odds per **1 unit** increase in original X.
- `beta1_SAS = beta1_logit * (π/√3)`; `beta1_Long = beta1_logit * (π/√3 + 1)`.
- `beta0_orig` aligns the intercept with original X-units.

**Dependencies:** `DISCO` (for severity); `logistf` optional (Firth comparator).

---

## Bayesian Estimation with MEP (multivariate; posterior means & CIs)

**When to use:** After your screen flags separation risk in multi-predictor settings, fit a logistic model with an MEP prior to stabilize estimation.

### Core sampler (fixed prior): `mep_function()`

```r
set.seed(1)
y <- c(0,0,0,0,1,1,1,1)
X <- cbind(
  X1 = c(-1.86, -0.81,  1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
  X2 = c( 0.52,  1.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09)
)
p <- ncol(X) + 1
fit <- mep_function(
  n_iter = 10000, burn_in = 1000, init_beta = rep(0.01, p), step_size = 0.4,
  X_orig = X, y = y, mu = rep(0, p), Sigma = diag(p), kappa = 1,
  scale_X = TRUE, ci_level = 0.95
)

# Working-space summary (includes Intercept) with CIs
head(fit$scaled_summary)

# Back-transformed effects with 95% CIs (per predictor; no intercept):
# Scaled (working/logit), A, SAS, Long
head(fit$standardized_coefs_back)
```

> Note: The output does **not** include an “M” scale; only **Scaled / A / SAS / Long** are returned.

### Prior grid search (auto-selected run): `MEP_latent()`

```r
set.seed(1)
y <- c(0,0,0,0,1,1,1,1)
X <- cbind(
  X1 = c(-1.86,-0.81, 1.32,-0.40, 0.91, 2.49, 0.34, 0.25),
  X2 = c( 0.52,-0.07, 0.60, 0.67,-1.39, 0.16,-1.40,-0.09)
)
gs <- MEP_latent(
  y, X,
  n_iter = 10000, burn_in = 1000, init_beta = NULL, step_size = 0.4,
  mu_vals = seq(-1, 1, by = 0.1),
  Sigma_list = list(diag(ncol(X)+1), diag(ncol(X)+1)*0.1, diag(ncol(X)+1)*0.5,
                    diag(ncol(X)+1)*2,  diag(ncol(X)+1)*5),
  kappa_vals = c(0.5, 1, 2),
  accept_window = c(0.3, 0.4), accept_target = 0.35,
  scale_X = TRUE, ci_level = 0.95
)

gs$best_settings        # chosen (mu, Sigma, kappa)
head(gs$scaled_summary) # working/logit space with CIs
head(gs$standardized_coefs_back) # A / SAS / Long with CIs
```

**Selection logic (summary):**
1. Filter runs by MH acceptance in `accept_window` (fallback to closest to `accept_target`).
2. Prefer higher `prop_matched` if < 0.90 is observed; else minimize mean absolute difference to a GLM coefficient-ratio reference (on the same scaling).

---

## Severity-Adaptive MEP for Mixture (numeric + factor predictors): `MEP_mixture()`

**What it does**
- Encodes `X` via `model.matrix(~ ., data = X)` (intercept dropped for slopes).
- For each **original** predictor (before encoding), computes a **DISCO severity** (on modeling rows).
- Maps severity \( s \in [0,1] \) to an **anchor slope prior sd** and **κ** (intercept uses a wide prior).
- Runs a small grid: **intercept mean offsets**, **global multiplier** on slope prior sds, and **κ** around the anchor average; chooses one run using an **acceptance window**, **posterior predictive agreement**, and (when available) **GLM ratio** closeness relative to a single **reference** predictor.
- Returns posterior **means** (not CIs) for standardized slopes and back-transforms (**A / SAS / Long**) per **encoded column**; also returns GLM/Firth/MEP **betas & ratios**.

**Factor handling & labels**
- Treatment contrasts with the **first level as baseline**.
- Two-level factor `X3` with levels `A` (baseline) and `B` yields dummy **`X3B`** (=1 for B, 0 for A). Its coefficient is **B vs A** (controlling for other predictors).
- Change baseline before calling to alter dummy labels:
  ```r
  X$X3 <- stats::relevel(X$X3, ref = "B")
  ```

**Example**
```r
set.seed(1)
y <- c(0,0,0,0, 1,1,1,1)
X <- data.frame(
  X1 = c(-1.86,-0.81, 1.32,-0.40, 0.91, 2.49, 0.34, 0.25),
  X2 = c( 0.52,-0.07, 0.60, 0.67,-1.39, 0.16,-1.40,-0.09),
  X3 = factor(c(rep("A",4), rep("B",4)))
)

fit <- MEP_mixture(
  y, X,
  n_iter_grid = 4000, burn_in_grid = 1000, step_size = 0.40,
  sigma_hi = 5, sigma_lo = 0.15,
  kappa_min = 1, kappa_max = 2.5,
  sigma_global_multipliers = c(0.1, 0.5, 1, 2, 5, 10),
  accept_window = c(0.30, 0.40),
  accept_target = 0.35,
  seed = 9, compare = TRUE
)

fit$ref_predictor       # chosen reference (original predictor index/name)
fit$posterior$effects   # means for Scaled / b_A_original / b_SAS_original / b_Long_original
fit$ratios$MEP_ratio_std
```

**Notes**
- Standardization/back-transforms use **unscaled encoded** column SDs; for a 0/1 dummy with prevalence \(p\), SD is \(\sqrt{p(1-p)}\).
- The **reference predictor** for ratios comes from original `X` (default = highest severity). If a factor, the denominator is its **first dummy**.
- Outputs for `effects` are **means without CIs** in this version.

---

## API cheatsheet
```r
# Univariate Detection
uni_separation(
  data, predictor, outcome = "Y",
  missing = c("complete","impute"),
  impute_args = list(
    numeric_method     = c("median","mean"),
    categorical_method = c("mode","missing"),
    logical_method     = c("mode","missing"),
    custom_fn          = NULL
  )
)

# Latent (multivariate) Detection
latent_separation(
  y, X,
  find_minimal = FALSE,
  test_combinations = FALSE,
  min_vars = 2, max_vars = NULL,
  mode = c("either","perfect","quasi"),
  stop_at_first = FALSE,
  missing = c("complete","impute"),
  impute_args = list(...),
  scale_X = FALSE     # set TRUE to z-score encoded predictors
)

# Severity-adaptive Univariate Bayes
EP_univariable(
  data, predictor, outcome = "y",
  missing = c("complete","impute"), impute_args = list(),
  n_iter = 20000, burn_in = 5000,
  step_hi = 0.30, step_lo = 0.12,
  ci_level = 0.95, compare = TRUE,
  return_draws = FALSE, seed = NULL,
  transform_beta1 = c("none","logit","SAS","Long"),
  sigma0 = 10, sigma1_hi = 5, sigma1_lo = 0.15,
  kappa_min = 1, kappa_max = 2.5,
  tune_threshold_hi = 0.45, tune_threshold_lo = 0.20, tune_interval = 500
)

# Core MEP logistic (fixed prior; returns CIs)
mep_function(
  n_iter, init_beta, step_size,
  X_orig, y, mu, Sigma, kappa,
  burn_in = 1000,
  scale_X = TRUE,
  ci_level = 0.95
)

# Grid search over (mu, Sigma, kappa); auto-select; returns CIs from the best run
MEP_latent(
  y, X,
  n_iter = 10000, burn_in = 1000, init_beta = NULL, step_size = 0.4,
  mu_vals = seq(-1, 1, by = 0.1),
  Sigma_list = list(...),
  kappa_vals = c(0.5, 1, 2),
  accept_window = c(0.3, 0.4), accept_target = 0.35,
  scale_X = TRUE, ci_level = 0.95
)

# Severity-anchored MEP with factors; grid + selection; returns means (no CIs)
MEP_mixture(
  y, X,
  missing = "complete",
  severity_missing = c("complete","impute"),
  impute_args = list(),
  n_iter_grid = 10000, burn_in_grid = 1000,
  step_size = 0.40,
  mu_intercept_offsets = seq(-1, 1, by = 0.2),
  sigma0_intercept = 10,
  sigma_global_multipliers = c(0.1, 0.5, 1, 2, 5, 10),
  sigma_hi = 5, sigma_lo = 0.15,
  kappa_min = 1, kappa_max = 2.5,
  kappa_delta = seq(-0.5, 0.5, by = 0.2),
  accept_window = c(0.30, 0.40), accept_target = 0.35,
  ref = NULL, compare = TRUE, seed = 2025,
  return_draws = FALSE
)
```

### Notes & assumptions
- Outcome is binary and will be normalized to `{0,1}` (supports logical or 2-level factor/character).
- Categorical predictors are handled directly (univariate) or via dummy encoding (latent / mixture).
- **MEP_mixture** outputs are per **encoded** column (e.g., `FactorLevel` dummies). Change baselines with `stats::relevel()` to alter dummy interpretation.

### Testing

The package includes unit tests covering:
- Univariate: no separation / perfect / quasi (numeric & categorical).
- Latent: no separation / perfect / quasi; minimal-subset search.
- Missingness: complete-case vs imputation (and custom imputer); rows_used/n_used always reported.

Run tests:
```r
devtools::test()
# or
testthat::test_dir("tests/testthat")
```
