library(testthat)
library(DiSCoB)

# -------------------------------------------------------------------
# Test on uni-separation
# -------------------------------------------------------------------

## No separation when mixed values ----------------------------------
test_that("No separation when mixed values", {
  df <- data.frame(Y = c(0,1,0,1), X = c(1,2,2,1))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "No separation problem")
})

test_that("No separation when mixed values", {
  df <- data.frame(Y = c("no", "yes", "no", "yes"), X = c(1,2,2,1))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "No separation problem")
})

test_that("No separation when mixed values", {
  df <- data.frame(Y = c("no", "yes", "no", "yes"), X = c("Male","Female","Female","Male"))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "No separation problem")
})


## Perfect separation when mixed values ------------------------------
test_that("Perfect separation when mixed values", {
  df <- data.frame(Y = c(0,1,0,1), X = c(-1,2,-2,1))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "Perfect separation")
})

test_that("Perfect separation when mixed values", {
  df <- data.frame(Y = c("no", "yes", "no", "yes"), X = c(-1,2,-2,1))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "Perfect separation")
})

test_that("Perfect separation when mixed values", {
  df <- data.frame(Y = c("no", "yes", "no", "yes"), X = factor(c("-1","2","-2","1")))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "Perfect separation")
})


test_that("Perfect separation when mixed values", {
  df <- data.frame(Y = c("no", "yes", "no", "yes"), X = c("High","Graduate","Middle","College"))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "Perfect separation")
})


## Quasi separation when mixed values ----------------------------------

test_that("Quasi separation when mixed values", {
  df <- data.frame(Y = c(0,1,0,1), X = c(-3,-2,-2,1))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "Quasi-complete separation")
})

test_that("Quasi separation when mixed values", {
  df <- data.frame(Y = factor(c("no", "yes", "no", "yes")), X = c(-3,-2,-2,1))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "Quasi-complete separation")
})

test_that("Quasi separation when mixed values", {
  df <- data.frame(Y = factor(c("no", "yes", "no", "yes")), X = factor(c("-3","-2","-2","1")))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "Quasi-complete separation")
})

test_that("Quasi separation when mixed values", {
  df <- data.frame(Y = factor(c("no", "yes", "no", "yes")), X = c("Low","Middle","Middle","High"))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "Quasi-complete separation")
})



# -------------------------------------------------------------------
# Test on latent-separation
# -------------------------------------------------------------------

## No latent separation ---------------------------------------------

test_that("Detect no latent separation", {
  set.seed(2025)
  y <- sample(0:1, 20, TRUE)
  X <- matrix(rnorm(20*3), 20, 3)
  res <- latent_separation(y, X)
  expect_equal(res$type, "no separation problem")
})


## Latent perfect separation ---------------------------------------------

test_that("Detect latent perfect separation", {
  set.seed(2025)
  y <- c(0,0,0,0,1,1,1,1)
  X <- cbind(
    X1 = c(-1.86, -0.81, 1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
    X2 = c( 0.52,  1.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09)
  )
  res <- latent_separation(y, X)
  expect_equal(res$type, "perfect separation")
})

test_that("Detect latent perfect separation", {
  set.seed(2025)
  y <- c("no","no","no","no","yes","yes","yes","yes")
  X <- cbind(
    X1 = c(-1.86, -0.81, 1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
    X2 = c( 0.52,  1.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09)
  )
  res <- latent_separation(y, X)
  expect_equal(res$type, "perfect separation")
})


test_that("Detect latent perfect separation", {
  set.seed(2025)
  y <- c(0,0,0,0,1,1,1,1)
  X <- cbind(
    X1 = c(-1.86, -0.81, 1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
    X2 = c( 0.52,  1.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09),
    X3 = c("Male", "Female", "Female", "Male","Male", "Female", "Female", "Male")
  )
  res <- latent_separation(y, X)
  expect_equal(res$type, "perfect separation")
})

## Latent quasi-complete separation ---------------------------------------------

test_that("Detect latent quasi-complete separation", {
  y <- c(0,0,0,0,1,1,1,1)
  X <- cbind(
    X1 = c(-1.86, -0.81, 1.32, 0.40,  0.91,  2.49,  0.34,  0.25),
    X2 = c( 0.52,  1.07, 0.60, -0.67, -1.39,  0.16, -1.40, -0.09)
  )
  res <- latent_separation(y, X)
  expect_equal(res$type, "quasi-complete separation")
})

test_that("Detect latent quasi-complete separation", {
  y <- c(0,0,0,0,1,1,1,1)
  X <- cbind(
    X1 = c(-1.86, -0.81, 1.32, 0.40,  0.91,  2.49,  0.34,  0.25),
    X2 = c( 0.52,  1.07, 0.60, -0.67, -1.39,  0.16, -1.40, -0.09),
    X3 = c("Male", "Female", "Female", "Male","Male", "Female", "Female", "Male")
  )
  res <- latent_separation(y, X)
  expect_equal(res$type, "quasi-complete separation")
})


test_that("Detect latent quasi-complete separation", {
  y <-  c("no","no","no","no","yes","yes","yes","yes")
  X <- cbind(
    X1 = c(-1.86, -0.81, 1.32, 0.40,  0.91,  2.49,  0.34,  0.25),
    X2 = c( 0.52,  1.07, 0.60, -0.67, -1.39,  0.16, -1.40, -0.09),
    X3 = c("Male", "Female", "Female", "Male","Male", "Female", "Female", "Male")
  )
  res <- latent_separation(y, X)
  expect_equal(res$type, "quasi-complete separation")
})



# -------------------------------------------------------------------
# Latent Separation:
# Subset of different combination of variables
# -------------------------------------------------------------------

## All Numeric -----------------------------------------------------

y <- c(0,0,0,0,1,1,1,1)
X <- cbind(
  X1 = c(-1.86, -0.81, 1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
  X2 = c( 0.52,  1.07, 0.60, 0.67, -1.39,  0.16, -1.40, -0.09),
  X3 = c(0, 0, 1, 1, 1, 1, 0, 0)
  )

# All 3 variables
res <- latent_separation(y, X, test_combinations = FALSE)
res

# Find all inclusion-minimal separating subsets (perfect or quasi)
res <- latent_separation(y, X, find_minimal = TRUE, min_vars = 2, mode = "either")
res$minimal_subsets


# Only perfect minimal subsets, stop at first:
res <- latent_separation(y, X, find_minimal = TRUE, min_vars = 2,
                         mode = "perfect", stop_at_first = TRUE)
res

## With Factors -----------------------------------------------------

y <- c(0,0,0,0,1,1,1,1)
X <- cbind(
  X1 = c(-1.86, -0.81, 1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
  X2 = c( 0.52,  1.07, 0.60, 0.67, -1.39,  0.16, -1.40, -0.09),
  X3 = c("Male", "Male", "Female", "Female", "Female", "Female", "Male", "Male")
)

# All 3 variables
res <- latent_separation(y, X, test_combinations = FALSE)
res

# Find all inclusion-minimal separating subsets (perfect or quasi)
res <- latent_separation(y, X, find_minimal = TRUE, min_vars = 2, mode = "either")
res$minimal_subsets

# Only perfect minimal subsets:
res <- latent_separation(y, X, find_minimal = TRUE, min_vars = 2,
                         mode = "perfect", stop_at_first = FALSE)
res

# Only perfect minimal subsets, stop at first:
res <- latent_separation(y, X, find_minimal = TRUE, min_vars = 2,
                         mode = "perfect", stop_at_first = TRUE)
res


## Set highest varuables ------------------------------------------------
y <- c(0,0,0,0,1,1,1,1)
X <- cbind(
  X1 = c(-1.86, -0.81, 1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
  X2 = c( 0.52,  1.07, 0.60, 0.67, -1.39,  0.16, -1.40, -0.09),
  X3 = c("Male", "Male", "Female", "Female", "Female", "Female", "Male", "Male"),
  X4 = c( 0.52,  0, 0.60, 0.67, -1.39,  0, -1.40, -0.09)
)

# All 4 variables
res <- latent_separation(y, X, test_combinations = FALSE)
res

# consider min as 2 variables and maximum as 3 vairables
res <- latent_separation(y, X, test_combinations = TRUE, find_minimal = TRUE, min_vars = 2, max_vars = 3)
res



# -------------------------------------------------------------------
# Test on data with missing values
# -------------------------------------------------------------------

set.seed(2025)
df <- data.frame(
  Y    = c(0,0,0,0,1,1,1,1, 0,1, NA, 1),
  X1   = c(0.5, 1.0, NA, 2.0, 5.0, 6.0, 7.0, NA, 1.5, 8.0, 9.0, 6.5),     # numeric w/ NAs
  X2   = c(10, 9, 8, NA, 6, 5, NA, 3, 2, 1, 0, NA),                        # numeric w/ NAs
  Race = factor(c("A","A","B", NA, "C","C","B","A", "B", NA, "C", "A")),   # multi-level w/ NAs
  L1   = c(TRUE, NA, FALSE, TRUE, TRUE, NA, FALSE, TRUE, FALSE, TRUE, TRUE, NA) # logical w/ NAs
)
df

# uni_separation() ---------------------------------------------------
# uni_separation() with complete-case handling
res_cc <- uni_separation(
  data = df,
  predictor = "X1",
  outcome   = "Y",
  missing   = "complete"   # <- drop rows with NA in Y or X1
)
res_cc$separation_type
res_cc$severity_score
res_cc$missing_info$n_used
res_cc$missing_info$rows_used

# uni_separation() with imputation (numeric median, categorical mode)
res_imp <- uni_separation(
  data = df,
  predictor = "X1",
  outcome   = "Y",
  missing   = "impute",
  impute_args = list(         # all optional; these are defaults
    numeric_method     = "median",
    categorical_method = "mode",
    logical_method     = "mode"
  )
)
res_imp$separation_type
res_imp$severity_score
res_cc$missing_info$n_used
res_cc$missing_info$rows_used

# Categorical predictor: treat missing as its own level
res_cat_missing_level <- uni_separation(
  data = df,
  predictor = "Race",
  outcome   = "Y",
  missing   = "impute",
  impute_args = list(categorical_method = "missing")  # add "Missing" level
)
res_cat_missing_level$separation_type
res_cat_missing_level$missing_info$params


# latent_separation() ---------------------------------------------------

# Minimal subsets with complete-case per subset
res_latent_cc <- latent_separation(
  y = df$Y,
  X = df[, c("X1","X2","Race","L1")],
  find_minimal  = TRUE,    # prune supersets once a minimal hit is found
  mode          = "either",
  missing       = "complete"
)

names(res_latent_cc$minimal_subsets)  # subset names that hit

sub_name <- names(res_latent_cc$minimal_subsets)[1]
res_latent_cc$minimal_subsets[[sub_name]]$type
res_latent_cc$minimal_subsets[[sub_name]]$missing_info$n_used
res_latent_cc$minimal_subsets[[sub_name]]$missing_info$rows_used


# Minimal subsets with imputation (numeric mean; missing category kept)

res_latent_imp <- latent_separation(
  y = df$Y,
  X = df[, c("X1","X2","Race","L1")],
  find_minimal  = TRUE,
  mode          = "either",
  missing       = "impute",
  impute_args   = list(
    numeric_method     = "mean",
    categorical_method = "missing",
    logical_method     = "mode"
  )
)

names(res_latent_imp$minimal_subsets)
sub_name <- names(res_latent_imp$minimal_subsets)[1]
res_latent_imp$minimal_subsets[[sub_name]]$type
res_latent_imp$minimal_subsets[[sub_name]]$missing_info$params
res_latent_imp$minimal_subsets[[sub_name]]$missing_info$n_used

# Using a custom imputer (plug any function you like)
my_imputer <- function(X) {
  X <- as.data.frame(X)
  for (j in names(X)) {
    if (is.numeric(X[[j]])) {
      fill <- stats::median(X[[j]], na.rm = TRUE)
      if (is.na(fill)) fill <- 0
      X[[j]][is.na(X[[j]])] <- fill
    } else if (is.logical(X[[j]])) {
      X[[j]][is.na(X[[j]])] <- FALSE
    } else { # treat as categorical
      f <- if (is.factor(X[[j]])) X[[j]] else factor(X[[j]], exclude = NULL)
      if (!("Missing" %in% levels(f))) levels(f) <- c(levels(f), "Missing")
      f[is.na(f)] <- "Missing"
      X[[j]] <- droplevels(f)
    }
  }
  X
}

res_latent_custom <- latent_separation(
  y = df$Y,
  X = df[, c("X1","X2","Race","L1")],
  find_minimal = TRUE,
  missing     = "impute",
  impute_args = list(custom_fn = my_imputer)
)

names(res_latent_custom$minimal_subsets)
res_latent_custom$minimal_subsets[[1]]$missing_info$params$custom_fn
res_latent_custom$minimal_subsets[[1]]$missing_info$n_used


# -------------------------------------------------------------------
# User Friendly Output
# -------------------------------------------------------------------

# install.packages(c("gt", "recipes", "VIM", "mice", "missForest"))

# install.packages("gt")
library(gt)

set.seed(2025)
df <- data.frame(
  Y    = c(0,0,0,0,1,1,1,1, 0,1, NA, 1),
  X1   = c(0.5, 1.0, NA, 2.0, 5.0, 6.0, 7.0, NA, 1.5, 8.0, 9.0, 6.5),
  X2   = c(10, 9, 8, NA, 6, 5, NA, 3, 2, 1, 0, NA),
  Race = factor(c("A","A","B", NA, "C","C","B","A", "B", NA, "C", "A")),
  L1   = c(TRUE, NA, FALSE, TRUE, TRUE, NA, FALSE, TRUE, FALSE, TRUE, TRUE, NA)
)

# Univariate -----------------------------------------------------------------

# 1. Univariate: complete-case
res_uni_cc <- uni_separation(df, "X1", "Y", missing = "complete")
gt_uni_separation(res_uni_cc, title = "Univariate (X1 vs Y) — Complete-case")

# All predictors vs Y, complete-case
tab_uni_all <- gt_uni_separation_all(df, outcome = "Y", missing = "complete")
gt::gtsave(tab_uni_all, "man/figures/readme-uni-gt-all.png")

# 2. Univariate: impute with explicit rules
res_uni_imp <- uni_separation(
  df, "X1", "Y", missing = "impute",
  impute_args = list(numeric_method = "median", categorical_method = "missing", logical_method = "mode")
)
gt_uni_separation(res_uni_imp, title = "Univariate (X1 vs Y) — Imputed", subtitle = "Median/mode; Missing=level")

# 3. Impute with custom choices, keep only Perfect/Quasi
gt_uni_separation_all(
  df, outcome = "Y", missing = "impute",
  impute_args = list(numeric_method = "mean", categorical_method = "missing"),
  only_hits = TRUE
)


# 4. Keep only Perfect/Quasi rows in the all-predictors summary
gt_uni_separation_all(
  df, outcome = "Y", missing = "impute",
  impute_args = list(backend = "recipes", backend_args = list(method = "median_mode")),
  only_hits = TRUE
)

# Latent ----------------------------------------------------------------------

# 1. Latent: minimal subsets, complete-case per subset
res_lat_cc <- latent_separation(
  y = df$Y,
  X = df[, c("X1","X2","Race","L1")],
  find_minimal = TRUE,
  missing = "complete"
)
tab_lat_com <- gt_latent_separation(res_lat_cc, title = "Latent Minimal Subsets — Complete-case", show_rows_used = TRUE)
gt::gtsave(tab_lat_com, "man/figures/readme-latent-gt-complete.png")

# 2. Latent: minimal subsets, imputed per subset
res_lat_imp <- latent_separation(
  y = df$Y,
  X = df[, c("X1","X2","Race","L1")],
  find_minimal = TRUE,
  missing = "impute",
  impute_args = list(numeric_method = "mean", categorical_method = "missing", logical_method = "mode")
)
tab_lat <- gt_latent_separation(res_lat_imp, title = "Latent Minimal Subsets — Imputed", subtitle = "Mean; Missing=level")
gt::gtsave(tab_lat, "man/figures/readme-latent-gt.png")

# 3. Impute via Mice
res_lat_mice <- latent_separation(
  y = df$Y,
  X = df[, c("X1","X2","Race","L1")],
  find_minimal = TRUE,
  mode = "either",
  missing = "impute",
  impute_args = list(
    backend = "mice",
    backend_args = list(m = 1, maxit = 5, printFlag = FALSE, seed = 2025)
  )
)
gt_latent_separation(res_lat_mice, title = "Latent Minimal Subsets — Imputed (Mice or fallback)")



