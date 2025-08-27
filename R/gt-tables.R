#' GT tables for separation results (no overlap column)
#' @keywords internal
#' @noRd

.format_impute_params <- function(params) {
  if (is.null(params)) return(NA_character_)
  if (!is.null(params$custom_fn)) return("custom imputer")
  kv <- c(
    paste0("numeric=", params$numeric_method  %||% NA_character_),
    # tolerate a mis-typed field if it exists in older results
    paste0("categorical=", params$categorical_method %||% params$categororical_method %||% NA_character_),
    paste0("logical=", params$logical_method %||% NA_character_)
  )
  paste(kv, collapse = "; ")
}

.type_icon <- function(x) {
  if (isTRUE(grepl("^Perfect", x))) return("perfect separation")
  if (isTRUE(grepl("^Quasi", x)))   return("quasi-complete separation")
  if (isTRUE(grepl("^No ", x)))     return("no problem")
  if (isTRUE(grepl("^Constant", x)))return("all constant")
  x
}

# New: palette function for gt::data_color(fn = ...)
.palette01_fn <- function(x) {
  pal <- grDevices::colorRampPalette(c("#e8f5e9", "#fff59d", "#ef9a9a"))(101)
  if (length(x) == 0) return(character(0))
  # clamp to [0,1], map to 1..101, keep NAs
  x_clamp <- pmin(pmax(x, 0), 1)
  idx <- ifelse(is.na(x_clamp), NA_integer_, as.integer(round(x_clamp * 100)) + 1L)
  out <- rep(NA_character_, length(x))
  ok <- !is.na(idx)
  out[ok] <- pal[idx[ok]]
  out
}

# Centralized string builder for Rows Used
.rows_used_string <- function(v, show_rows_used = FALSE) {
  if (length(v) == 0) return("—")
  # Detect "all subjects" as exactly 1:N
  is_all <- length(v) == max(v) && identical(v, seq_len(max(v)))
  if (is_all) {
    # Always show the full list when all subjects are used
    return(paste(v, collapse = ", "))
  }
  # Otherwise, follow the preview toggle
  if (isTRUE(show_rows_used)) {
    paste0(paste(head(v, 10), collapse = ", "),
           if (length(v) > 10) " …" else "")
  } else {
    "—"
  }
}

#' Tidy a uni_separation() result (no overlap column)
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
    # overlap removed
    missing_method     = mi$method %||% NA_character_,
    impute_params      = .format_impute_params(mi$params %||% NULL),
    n_used             = as.integer(mi$n_used %||% NA_integer_),
    rows_used          = list(mi$rows_used %||% integer())
  )
}

#' Create a gt table for a univariate result (no overlap column)
#' @param res Result of uni_separation()
#' @param title,subtitle Title/subtitle strings
#' @param digits Numeric formatting decimals
#' @param show_rows_used If TRUE, show a preview of row indices
#' @export
gt_uni_separation <- function(res, title = "Univariate Separation",
                              subtitle = NULL, digits = 3,
                              show_rows_used = FALSE) {
  stopifnot(requireNamespace("gt", quietly = TRUE))
  df <- tidy_uni_separation(res)
  df$Type <- vapply(df$type, .type_icon, character(1))
  df$rows_used_str <- vapply(df$rows_used, .rows_used_string, character(1),
                             show_rows_used = show_rows_used)

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
      # overlap removed
      missing_method     = "Missing Method",
      impute_params      = "Imputation Params",
      n_used             = "N Used",
      rows_used_str      = "Rows Used"
    ) |>
    gt::cols_hide(columns = c(type, rows_used)) |>
    gt::fmt_number(
      columns  = c(separation_index, severity_score, boundary_threshold),
      decimals = digits
    ) |>
    gt::fmt_markdown(columns = c(Type)) |>
    # Updated: use fn= instead of colors=
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

#' Tidy latent_separation() results (minimal-subsets or exhaustive)
#' @keywords internal
#' @noRd
tidy_latent_separation <- function(res) {
  as_row <- function(name, item) {
    mi <- item$missing_info %||% list()
    tibble::tibble(
      subset_name    = name,
      vars           = paste(item$vars %||% NA_character_, collapse = ", "),
      k              = length(item$vars %||% character()),
      type           = item$type %||% res$type %||% NA_character_,
      removed        = paste(item$removed %||% character(), collapse = ", "),
      missing_method = mi$method %||% NA_character_,
      impute_params  = .format_impute_params(mi$params %||% NULL),
      n_used         = as.integer(mi$n_used %||% NA_integer_),
      rows_used      = list(mi$rows_used %||% integer())
    )
  }

  if (!is.null(res$minimal_subsets)) {
    lst <- res$minimal_subsets
    if (!length(lst)) {
      return(tibble::tibble(subset_name=character(), vars=character(), k=integer(),
                            type=character(), removed=character(),
                            missing_method=character(), impute_params=character(),
                            n_used=integer(), rows_used=list()))
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

#' Create a gt table for latent results
#' @param res Result of latent_separation()
#' @param title,subtitle Title/subtitle strings
#' @param sort_by Columns to sort by (in priority order)
#' @param show_rows_used If TRUE, show preview of row indices per subset
#' @export
gt_latent_separation <- function(res, title = "Latent Separation: Minimal Subsets",
                                 subtitle = NULL, sort_by = c("type","k","n_used"),
                                 show_rows_used = FALSE) {
  stopifnot(requireNamespace("gt", quietly = TRUE))
  df <- tidy_latent_separation(res)
  if (nrow(df) == 0L) {
    return(
      gt::gt(tibble::tibble(Note = "No separating subset found.")) |>
        gt::tab_header(title = title, subtitle = subtitle)
    )
  }

  df$Type <- vapply(df$type, .type_icon, character(1))
  df$type_order <- dplyr::case_when(
    grepl("^Perfect", df$type) ~ 1L,
    grepl("^Quasi",   df$type) ~ 2L,
    grepl("^No ",     df$type) ~ 3L,
    TRUE ~ 9L
  )

  # sort by user preference (stable)
  for (key in rev(sort_by)) if (key %in% names(df)) df <- df[order(df[[key]]), , drop = FALSE]

  df$rows_used_str <- vapply(df$rows_used, .rows_used_string, character(1),
                             show_rows_used = show_rows_used)

  gt::gt(df) |>
    gt::tab_header(title = title, subtitle = subtitle) |>
    gt::cols_label(
      subset_name   = "Subset",
      vars          = "Variables",
      k             = "# of predictors",
      Type          = "Separation",
      removed       = "Removed to Perfect (quasi)",
      missing_method= "Missing Method",
      impute_params = "Imputation Params",
      n_used        = "N Used",
      rows_used_str = "Rows Used"
    ) |>
    gt::cols_hide(columns = c(type, type_order, rows_used)) |>
    gt::fmt_markdown(columns = c(Type)) |>
    gt::tab_spanner(label = "Subset", columns = c(subset_name, vars, k)) |>
    gt::tab_spanner(label = "Missing-Data Handling", columns = c(missing_method, impute_params, n_used)) |>
    gt::tab_options(table.font.size = gt::px(14), data_row.padding = gt::px(4))
}
