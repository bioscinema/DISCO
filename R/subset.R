# -- internal: single-subset separation check -----------------------------------
# mode: "perfect" => perfect only
#       "quasi"   => quasi only (via leave-one-out)
#       "either"  => perfect first, else quasi
.check_sep <- function(y, X, mode = c("either","perfect","quasi")) {
  mode <- match.arg(mode)
  X <- as.matrix(X)

  # perfect?
  full_res <- feasibility_check_LP(y, X)
  perfect_hit <- (full_res$ge$status == 0 || full_res$le$status == 0)
  if (mode %in% c("either","perfect") && perfect_hit) {
    return(list(hit = TRUE, type = "perfect separation", removed = NULL))
  }
  if (mode == "perfect") {
    return(list(hit = FALSE, type = "no separation problem", removed = NULL))
  }

  # quasi via leave-one-out perfect check
  quasi_inds <- which(vapply(seq_along(y), function(i) {
    sub <- feasibility_check_LP(y[-i], X[-i,, drop = FALSE])
    sub$ge$status == 0 || sub$le$status == 0
  }, logical(1)))

  if (length(quasi_inds) > 0) {
    return(list(hit = TRUE, type = "quasi-complete separation", removed = quasi_inds))
  }
  list(hit = FALSE, type = "no separation problem", removed = NULL)
}
