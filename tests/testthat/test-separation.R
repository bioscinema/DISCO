library(testthat)
library(DiSCoB)

# -------------------------------------------------------------------
# Test on uni-separation
# -------------------------------------------------------------------

test_that("No separation when mixed values", {
  df <- data.frame(Y = c(0,1,0,1), X = c(1,2,2,1))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "No separation problem")
})

test_that("Perfect separation when mixed values", {
  df <- data.frame(Y = c(0,1,0,1), X = c(-1,2,-2,1))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "Perfect separation")
})

test_that("Quasi separation when mixed values", {
  df <- data.frame(Y = c(0,1,0,1), X = c(-3,-2,-2,1))
  res <- uni_separation(df, "X")
  expect_equal(res$separation_type, "Quasi-complete separation")
})


# -------------------------------------------------------------------
# Test on latent-separation
# -------------------------------------------------------------------
test_that("Detect no latent separation", {
  set.seed(2025)
  y <- sample(0:1, 20, TRUE)
  X <- matrix(rnorm(20*3), 20, 3)
  res <- latent_separation(y, X)
  expect_equal(res$type, "no separation problem")
})

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

test_that("Detect latent quasi-complete separation", {
  y <- c(0,0,0,0,1,1,1,1)
  X <- cbind(
    X1 = c(-1.86, -0.81, 1.32, 0.40,  0.91,  2.49,  0.34,  0.25),
    X2 = c( 0.52,  1.07, 0.60, -0.67, -1.39,  0.16, -1.40, -0.09)
  )
  res <- latent_separation(y, X)
  expect_equal(res$type, "quasi-complete separation")
})


# -------------------------------------------------------------------
# Example: subset of different combination of variables
# -------------------------------------------------------------------

y <- c(0,0,0,0,1,1,1,1)
X <- cbind(
  X1 = c(-1.86, -0.81, 1.32, 0.40,  0.91,  2.49,  0.34,  0.25),
  X2 = c( 0.52,  1.07, 0.60, -0.67, -1.39,  0.16, -1.40, -0.09),
  X3 = c(0, 0, 1, 1, 1, 1, 0, 0)
  )

# All subset diagnostics (quasi‐complete and perfect)
res <- latent_separation(y, X, test_combinations = TRUE, min_vars = 2)
res

# Only the subsets exhibiting perfect separation on all subjects
res <- latent_separation(y, X, test_combinations = TRUE, min_vars = 2, only_perfect = TRUE)
res
