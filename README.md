# DISCO
**DI**agnosis of **S**eparation and **C**orrection of **O**dds-ratio inflation in logistic regression

## Overview
**DISCO** provides tools to:
### Separation Diagnosis
- `uni_separation()` — univariate detector using Hungarian matching + a **vectorized Rand index**, **non-negative boundary threshold**, and a **continuous severity score** in [0,1].
- `latent_separation()` — multivariate, **LP-based** feasibility checks for **perfect** or **quasi-complete** separation; can optionally test **all variable subsets**.

### Estimation Correction

## Installation
Install the development version directly from GitHub:
```r
# install.packages("devtools")
devtools::install_github("bioscinema/DISCO")
library(DISCO)
```

## Separation Diagnosis Examples

### Univariate Separation
```r

# No separation
df <- data.frame(Y = c(0,1,0,1), X = c(1,2,2,1))
uni_separation(df, "X")
#> separation_type: "No separation problem"

# Perfect separation
df <- data.frame(Y = c(0,1,0,1), X = c(-1,2,-2,1))
uni_separation(df, "X")
#> separation_type: "Perfect separation"

# Quasi-complete separation
df <- data.frame(Y = c(0,1,0,1), X = c(-3,-2,-2,1))
uni_separation(df, "X")
#> separation_type: "Quasi-complete separation"
```

### Latent (multivariate) Separation

```r
set.seed(2025)
y <- sample(0:1, 20, TRUE)
X <- matrix(rnorm(20*3), 20, 3)

latent_separation(y, X)
#> type: "no separation problem"


# Example with true separation
y <- c(0,0,0,0,1,1,1,1)
X <- cbind(
  X1 = c(-1.86, -0.81, 1.32, -0.40,  0.91,  2.49,  0.34,  0.25),
  X2 = c( 0.52,  1.07,  0.60,  0.67, -1.39,  0.16, -1.40, -0.09)
)

latent_separation(y, X)
#> type: "perfect separation"


# Quasi-complete separation example
y <- c(0,0,0,0,1,1,1,1)
X <- cbind(
  X1 = c(-1.86, -0.81, 1.32, 0.40,  0.91,  2.49,  0.34,  0.25),
  X2 = c( 0.52,  1.07, 0.60, -0.67, -1.39,  0.16, -1.40, -0.09)
)

latent_separation(y, X)
#> type: "quasi-complete separation"


# ------------------------------------------------------------
# Subset search
# ------------------------------------------------------------

y <- c(0,0,0,0,1,1,1,1)
X <- cbind(
  X1 = c(-1.86, -0.81, 1.32, 0.40,  0.91,  2.49,  0.34,  0.25),
  X2 = c( 0.52,  1.07, 0.60, -0.67, -1.39,  0.16, -1.40, -0.09),
  X3 = c(0, 0, 1, 1, 1, 1, 0, 0)
)

# All subset diagnostics (quasi and perfect)
res_all <- latent_separation(y, X, test_combinations = TRUE, min_vars = 2)
res_all

# Only subsets with perfect separation
res_perfect <- latent_separation(y, X, test_combinations = TRUE, min_vars = 2, only_perfect = TRUE)
res_perfect
```
