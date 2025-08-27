#' GT table: univariate separation across all predictors
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
#' @param show_rows_used If TRUE, preview row indices (unless all rows are used,
#'   in which case the full 1..N is shown automatically).
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

  # Run uni_separation() per predictor, tidy, bind
  res_list <- lapply(predictors, function(p) {
    out <- try(
      uni_separation(
        data       = data,
        predictor  = p,
        outcome    = outcome,
        missing    = missing,
        impute_args= impute_args
      ),
      silent = TRUE
    )
    if (inherits(out, "try-error")) return(NULL)
    tidy_uni_separation(out)
  })
  res_list <- Filter(Negate(is.null), res_list)
  if (!length(res_list)) {
    return(
      gt::gt(tibble::tibble(Note = "No predictors or no valid results.")) |>
        gt::tab_header(title = title, subtitle = subtitle)
    )
  }
  df <- dplyr::bind_rows(res_list)

  # Filtering options
  if (!include_constant) df <- df[!grepl("^Constant", df$type), , drop = FALSE]
  if (isTRUE(only_hits)) df <- df[grepl("^(Perfect|Quasi)", df$type), , drop = FALSE]

  if (nrow(df) == 0L) {
    return(
      gt::gt(tibble::tibble(Note = "No rows after filtering.")) |>
        gt::tab_header(title = title, subtitle = subtitle)
    )
  }

  # Decorate
  df$Type <- vapply(df$type, .type_icon, character(1))
  df$type_order <- dplyr::case_when(
    grepl("^Perfect", df$type) ~ 1L,
    grepl("^Quasi",   df$type) ~ 2L,
    grepl("^No ",     df$type) ~ 3L,
    TRUE ~ 9L
  )
  df$rows_used_str <- vapply(df$rows_used, .rows_used_string, character(1),
                             show_rows_used = show_rows_used)

  # Sort: best to worst by type, then severity desc, then index desc
  df <- df[order(df$type_order, -df$severity_score, -df$separation_index,
                 df$predictor), , drop = FALSE]

  gt::gt(df) |>
    gt::tab_header(title = title, subtitle = subtitle) |>
    gt::cols_label(
      predictor          = "Predictor",
      outcome            = "Outcome",
      Type               = "Separation",
      separation_index   = "Separation Index",
      severity_score     = "Severity",
      boundary_threshold = "Boundary Threshold",
      single_tie         = "Single-Tie Boundary",
      tie_count          = "Tie Count",
      missing_method     = "Missing Method",
      impute_params      = "Imputation Params",
      n_used             = "N Used",
      rows_used_str      = "Rows Used"
    ) |>
    gt::cols_hide(columns = c(type, type_order, rows_used)) |>
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
    gt::tab_spanner(
      label   = "Missing-Data Handling",
      columns = c(missing_method, impute_params, n_used)
    ) |>
    gt::fmt(
      columns = c(single_tie),
      fns = function(x) ifelse(is.na(x), "—", ifelse(x, "Yes", "No"))
    ) |>
    gt::tab_options(table.font.size = gt::px(14), data_row.padding = gt::px(4))
}
