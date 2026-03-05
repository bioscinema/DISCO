# ======= Drop all "Missing-Data Handling" columns from ALL gt outputs =======
# Changes:
# - Remove missing_method / impute_params / n_used from tidy_* outputs
# - Remove "Missing-Data Handling" spanners and labels
# - Simplify rows_used display (no longer depends on missing method)

# Build a compact "Rows Used" preview string
.rows_used_string <- function(v, show_rows_used = FALSE) {
  if (length(v) == 0) return("—")

  # Detect "all subjects" as exactly 1:N
  is_all <- length(v) == max(v) && identical(v, seq_len(max(v)))
  if (is_all) return(paste(v, collapse = ", "))

  if (isTRUE(show_rows_used)) {
    paste0(
      paste(head(v, 10), collapse = ", "),
      if (length(v) > 10) " …" else ""
    )
  } else {
    "—"
  }
}

# Always show full list if it's exactly 1..N; otherwise respect preview toggle
.rows_used_display <- function(v, show_rows_used = FALSE) {
  .rows_used_string(v, show_rows_used = show_rows_used)
}

.to_titlecase <- function(x) {
  if (is.null(x)) return(x)
  tools::toTitleCase(as.character(x))
}

.type_icon <- function(x) {
  if (isTRUE(grepl("^Perfect",  x))) return("Perfect Separation")
  if (isTRUE(grepl("^Quasi",    x))) return("Quasi-Complete Separation")
  if (isTRUE(grepl("^No ",      x))) return("No Problem")
  if (isTRUE(grepl("^Constant", x))) return("All Constant")
  .to_titlecase(x)
}

.palette01_fn <- function(x) {
  pal <- grDevices::colorRampPalette(c("#e8f5e9", "#fff59d", "#ef9a9a"))(101)
  if (length(x) == 0) return(character(0))

  out <- rep("#f2f2f2", length(x))  # grey for NA
  x_clamp <- pmin(pmax(x, 0), 1)
  idx <- ifelse(is.na(x_clamp), NA_integer_, as.integer(round(x_clamp * 100)) + 1L)

  ok <- !is.na(idx)
  out[ok] <- pal[idx[ok]]
  out
}

# ---------- Univariate result (single predictor) ----------

#' Tidy a uni_separation() result (missing-handling columns removed)
#' @keywords internal
#' @noRd
tidy_uni_separation <- function(res) {
  mi <- res$missing_info %||% list()
  tibble::tibble(
    predictor          = res$predictor %||% NA_character_,
    outcome            = res$outcome   %||% NA_character_,
    type               = res$separation_type,
    separation_index   = as.numeric(res$separation_index),
    severity_score     = as.numeric(res$severity_score),
    boundary_threshold = as.numeric(res$boundary_threshold),
    single_tie         = isTRUE(res$single_tie_boundary),
    tie_count          = as.integer(res$tie_rows_boundary %||% 0L),
    rows_used          = list(mi$rows_used %||% integer())
  )
}

#' Create a gt table for a univariate result (missing-handling columns removed)
#' @param res Result of uni_separation()
#' @param title,subtitle Title/subtitle strings
#' @param digits Numeric formatting decimals
#' @param show_rows_used If TRUE, show a preview of row indices (unless all rows are used)
#' @export
gt_uni_separation <- function(
    res,
    title = "Univariate Separation",
    subtitle = NULL,
    digits = 3,
    show_rows_used = FALSE
) {
  stopifnot(requireNamespace("gt", quietly = TRUE))
  df <- tidy_uni_separation(res)
  df$Type <- vapply(df$type, .type_icon, character(1))

  df$rows_used_str <- vapply(
    df$rows_used,
    .rows_used_display,
    character(1),
    show_rows_used = show_rows_used
  )

  # Title-case character body cells
  char_cols <- vapply(df, is.character, logical(1))
  df[names(df)[char_cols]] <- lapply(df[names(df)[char_cols]], .to_titlecase)

  gt::gt(df) |>
    gt::tab_header(title = .to_titlecase(title), subtitle = .to_titlecase(subtitle)) |>
    gt::cols_label(
      predictor          = "Predictor",
      outcome            = "Outcome",
      Type               = "Separation",
      separation_index   = "Separation Index",
      severity_score     = "Severity",
      boundary_threshold = "Boundary Threshold",
      single_tie         = "Single-Tie Boundary",
      tie_count          = "Tie Count",
      rows_used_str      = "Rows Used (Original Indices)"
    ) |>
    gt::cols_hide(columns = c(type, rows_used)) |>
    gt::fmt_number(
      columns  = c(separation_index, severity_score, boundary_threshold),
      decimals = digits
    ) |>
    gt::fmt_markdown(columns = c(Type)) |>
    gt::data_color(columns = c(severity_score), fn = .palette01_fn) |>
    gt::tab_spanner(
      label   = "Indices",
      columns = c(separation_index, severity_score, boundary_threshold)
    ) |>
    gt::fmt(
      columns = c(single_tie),
      fns = function(x) ifelse(is.na(x), "—", ifelse(x, "Yes", "No"))
    ) |>
    gt::tab_options(table.font.size = gt::px(14), data_row.padding = gt::px(4))
}

# ---------- Latent (minimal subsets / exhaustive) ----------

#' Tidy latent_separation() results (missing-handling columns removed)
#' @keywords internal
#' @noRd
tidy_latent_separation <- function(res) {
  as_row <- function(name, item) {
    mi <- item$missing_info %||% list()
    di <- item$diagnostics %||% list()

    n_diag <- di$n %||% mi$n_used %||% NA_integer_
    K_rel  <- di$K_relax %||% NA_real_

    K_rel_norm <- if (!is.null(K_rel) && !is.null(n_diag) &&
                      is.finite(K_rel) && n_diag > 0L) {
      as.numeric(K_rel) / as.numeric(n_diag)
    } else {
      NA_real_
    }

    score_val <- if (!is.na(K_rel_norm)) 1 - K_rel_norm else NA_real_

    tibble::tibble(
      subset_name   = name,
      vars          = if (is.null(item$vars)) "All predictors" else paste(item$vars, collapse = ", "),
      k             = if (is.null(item$vars)) NA_integer_ else length(item$vars),
      type          = item$type %||% res$type %||% NA_character_,
      removed       = paste(item$removed %||% character(), collapse = ", "),
      K_relax       = as.numeric(K_rel %||% NA_real_),
      K_relax_per_n = round(as.numeric(K_rel_norm), 3),
      score         = score_val,
      n             = as.integer(n_diag),
      rows_used     = list(mi$rows_used %||% integer())
    )
  }

  if (!is.null(res$minimal_subsets)) {
    lst <- res$minimal_subsets
    if (!length(lst)) {
      return(
        tibble::tibble(
          subset_name   = character(),
          vars          = character(),
          k             = integer(),
          type          = character(),
          removed       = character(),
          K_relax       = numeric(),
          K_relax_per_n = numeric(),
          score         = numeric(),
          n             = integer(),
          rows_used     = list()
        )
      )
    }
    dplyr::bind_rows(Map(as_row, names(lst), unname(lst)))
  } else if (!is.null(res$type)) {
    as_row("(all predictors)", res)
  } else if (is.list(res) && length(res) && !is.null(res[[1]]$type)) {
    dplyr::bind_rows(Map(as_row, names(res), unname(res)))
  } else {
    stop("Unrecognized structure for latent_separation() result.")
  }
}

#' Create a gt table for latent results (missing-handling columns removed)
#' @param res Result of latent_separation()
#' @param title,subtitle Title/subtitle strings
#' @param sort_by Columns to sort by (in priority order)
#' @param show_rows_used If TRUE, show preview of row indices per subset
#' @export
gt_latent_separation <- function(
    res,
    title = "Latent Separation: Minimal Subsets",
    subtitle = NULL,
    sort_by = c("type_order", "K_relax_per_n", "k", "n"),
    show_rows_used = FALSE,
    digits = 3
) {
  stopifnot(requireNamespace("gt", quietly = TRUE))
  df <- tidy_latent_separation(res)

  if (nrow(df) == 0L) {
    return(
      gt::gt(tibble::tibble(Note = "No Separating Subset Found.")) |>
        gt::tab_header(title = .to_titlecase(title), subtitle = .to_titlecase(subtitle))
    )
  }

  df$Type <- vapply(df$type, .type_icon, character(1))
  df$type_order <- dplyr::case_when(
    grepl("^Perfect", df$type) ~ 1L,
    grepl("^Quasi",   df$type) ~ 2L,
    grepl("^No ",     df$type) ~ 3L,
    TRUE                        ~ 9L
  )

  for (key in rev(sort_by)) {
    if (key %in% names(df)) {
      df <- df[order(df[[key]]), , drop = FALSE]
    }
  }

  df$rows_used_str <- vapply(
    df$rows_used,
    .rows_used_display,
    character(1),
    show_rows_used = show_rows_used
  )

  char_cols <- vapply(df, is.character, logical(1))
  df[names(df)[char_cols]] <- lapply(df[names(df)[char_cols]], .to_titlecase)

  gt::gt(df) |>
    gt::tab_header(title = .to_titlecase(title), subtitle = .to_titlecase(subtitle)) |>
    gt::cols_label(
      subset_name   = "Subset",
      vars          = "Variables",
      k             = "# Of Predictors",
      Type          = "Separation",
      K_relax       = "K_relax",
      score         = "Score",
      n             = "n In LP",
      rows_used_str = "Rows Used (Original Indices)"
    ) |>
    gt::cols_hide(columns = c(type, type_order, rows_used, removed)) |>
    gt::fmt_number(columns = c(K_relax), decimals = 0) |>
    gt::fmt_number(columns = c(score), decimals = digits) |>
    gt::fmt_markdown(columns = c(Type)) |>
    gt::sub_missing(columns = c(score), missing_text = "No separation") |>
    gt::data_color(columns = c(score), fn = .palette01_fn) |>
    gt::tab_spanner(label = "Subset",   columns = c(subset_name, vars, k)) |>
    gt::tab_spanner(label = "Severity", columns = c(K_relax, score, n)) |>
    gt::tab_options(table.font.size = gt::px(14), data_row.padding = gt::px(4))
}

# ---------- Univariate across all predictors ----------

#' GT table: univariate separation across all predictors (missing-handling columns removed)
#'
#' Runs uni_separation() for each predictor and renders a single gt table.
#'
#' @param data A data.frame/tibble.
#' @param outcome Outcome column name (default "Y").
#' @param predictors Optional character vector of predictor names. If NULL,
#'   uses all columns except `outcome`.
#' @param missing "complete" or "impute" — passed to uni_separation().
#' @param impute_args List of imputation args — passed to uni_separation().
#' @param include_constant If FALSE (default), drop rows with "Constant outcome/predictor".
#' @param only_hits If TRUE, keep only "Perfect" or "Quasi" rows.
#' @param digits Number formatting for numeric columns (default 3).
#' @param show_rows_used If TRUE, preview row indices (unless all rows are used)
#' @param title,subtitle Title/subtitle for the table.
#' @export
gt_uni_separation_all <- function(
    data,
    outcome = "Y",
    predictors = NULL,
    missing = c("complete","impute"),
    impute_args = list(),
    include_constant = FALSE,
    only_hits = FALSE,
    digits = 3,
    show_rows_used = FALSE,
    title = "Univariate Separation — All Predictors",
    subtitle = NULL
) {
  stopifnot(requireNamespace("gt", quietly = TRUE))
  missing <- match.arg(missing)

  if (is.null(predictors)) {
    predictors <- setdiff(names(data), outcome)
  }

  res_list <- lapply(predictors, function(p) {
    out <- try(
      uni_separation(
        data        = data,
        predictor   = p,
        outcome     = outcome,
        missing     = missing,
        impute_args = impute_args
      ),
      silent = TRUE
    )
    if (inherits(out, "try-error")) return(NULL)
    tidy_uni_separation(out)
  })
  res_list <- Filter(Negate(is.null), res_list)

  if (!length(res_list)) {
    return(
      gt::gt(tibble::tibble(Note = "No Predictors Or No Valid Results.")) |>
        gt::tab_header(title = .to_titlecase(title), subtitle = .to_titlecase(subtitle))
    )
  }

  df <- dplyr::bind_rows(res_list)

  if (!include_constant) df <- df[!grepl("^Constant", df$type), , drop = FALSE]
  if (isTRUE(only_hits)) df <- df[grepl("^(Perfect|Quasi)", df$type), , drop = FALSE]

  if (nrow(df) == 0L) {
    return(
      gt::gt(tibble::tibble(Note = "No Rows After Filtering.")) |>
        gt::tab_header(title = .to_titlecase(title), subtitle = .to_titlecase(subtitle))
    )
  }

  df$Type <- vapply(df$type, .type_icon, character(1))
  df$type_order <- dplyr::case_when(
    grepl("^Perfect", df$type, ignore.case = TRUE) ~ 1L,
    grepl("^Quasi",   df$type, ignore.case = TRUE) ~ 2L,
    grepl("^No ",     df$type, ignore.case = TRUE) ~ 3L,
    TRUE ~ 9L
  )

  df$rows_used_str <- vapply(
    df$rows_used,
    .rows_used_display,
    character(1),
    show_rows_used = show_rows_used
  )

  df <- df[order(df$type_order, -df$severity_score, -df$separation_index, df$predictor), , drop = FALSE]

  char_cols <- vapply(df, is.character, logical(1))
  df[names(df)[char_cols]] <- lapply(df[names(df)[char_cols]], .to_titlecase)

  gt::gt(df) |>
    gt::tab_header(title = .to_titlecase(title), subtitle = .to_titlecase(subtitle)) |>
    gt::cols_label(
      predictor          = "Predictor",
      outcome            = "Outcome",
      Type               = "Separation",
      separation_index   = "Separation Index",
      severity_score     = "Severity",
      boundary_threshold = "Boundary Threshold",
      single_tie         = "Single-Tie Boundary",
      tie_count          = "Tie Count",
      rows_used_str      = "Rows Used (Original Indices)"
    ) |>
    gt::cols_hide(columns = c(type, type_order, rows_used)) |>
    gt::fmt_number(columns = c(separation_index, severity_score, boundary_threshold), decimals = digits) |>
    gt::fmt_markdown(columns = c(Type)) |>
    gt::data_color(columns = c(severity_score), fn = .palette01_fn) |>
    gt::tab_spanner(
      label   = "Indices",
      columns = c(separation_index, severity_score, boundary_threshold)
    ) |>
    gt::fmt(
      columns = c(single_tie),
      fns = function(x) ifelse(is.na(x), "—", ifelse(x, "Yes", "No"))
    ) |>
    gt::tab_options(table.font.size = gt::px(14), data_row.padding = gt::px(4))
}
