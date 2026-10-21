#' Compare a table with a reference table
#'
#' Used to show that refactored code reproduces a reference analysis. Rows are matched on key columns;
#' character columns must be identical; numeric columns must agree within a relative tolerance (absolute near
#' zero). When the reference was written rounded, `round_ref` declares it, and the new values are rounded the
#' same way before comparison; such matches are reported as `different_explained` with the rounding rule, never
#' as identical. Tolerances are inputs fixed before the comparison, not tuned after it.
#'
#' @param ref,new `data.frame`s with the same column names (map them before calling).
#' @param keys `character`. Key columns that identify a row.
#' @param rel_tol,abs_tol `numeric(1)`. Relative and absolute tolerance for numeric columns.
#' @param round_ref Named `list` column -> rounding rule: `list(digits = d)` (`round`) or
#'   `list(signif = d)` (significant digits).
#' @param columns `character` or `NULL`. Columns to compare (default: all non-key columns of `ref`).
#' @return `list(status, max_rel_diff, n_ref, n_new, n_matched, n_mismatch, mismatched_columns, note)`;
#'   status is `"identical"`, `"equal_within_tolerance"`, `"different_explained"` or `"different"`.
#' @examples
#' a <- data.frame(id = c("x", "y"), v = c(1, 2))
#' compare_tables(a, transform(a, v = v * (1 + 1e-12)), keys = "id")$status
#' @export
compare_tables <- function(ref, new, keys, rel_tol = 1e-8, abs_tol = 1e-12, round_ref = list(), columns = NULL) {
  columns <- columns %||% setdiff(names(ref), keys)
  missing <- setdiff(c(keys, columns), names(new))
  if (length(missing) > 0) {
    return(list(
      status = "different", max_rel_diff = NA_real_, n_ref = nrow(ref), n_new = nrow(new), n_matched = 0L,
      n_mismatch = NA_integer_, mismatched_columns = missing,
      note = paste("missing columns in new:", paste(missing, collapse = ", "))
    ))
  }
  kr <- do.call(paste, c(ref[keys], sep = "\r"))
  kn <- do.call(paste, c(new[keys], sep = "\r"))
  if (anyDuplicated(kr) || anyDuplicated(kn)) {
    abort("compare_tables(): keys ", paste(keys, collapse = ","), " are not unique")
  }
  i <- match(kr, kn)
  n_matched <- sum(!is.na(i))
  bad_cols <- character()
  max_rel <- 0
  exact <- TRUE
  rounded_used <- FALSE
  for (col in columns) {
    a <- ref[[col]]
    b <- new[[col]][i]
    if (is.numeric(a) && is.numeric(b)) {
      rule <- round_ref[[col]]
      if (!is.null(rule)) {
        b2 <- if (!is.null(rule$digits)) round(b, rule$digits) else signif(b, rule$signif)
        ok <- (is.na(a) & is.na(b2)) | (!is.na(a) & !is.na(b2) & abs(a - b2) <= pmax(abs_tol, rel_tol * abs(a)))
        if (!all(ok)) bad_cols <- c(bad_cols, col)
        if (!isTRUE(all(a == b, na.rm = TRUE)) || !identical(is.na(a), is.na(b))) rounded_used <- TRUE
        next
      }
      both_na <- is.na(a) & is.na(b)
      d <- abs(a - b)
      rel <- ifelse(abs(a) > abs_tol, d / abs(a), d)
      ok <- both_na | (!is.na(d) & (d <= abs_tol | rel <= rel_tol))
      if (any(!both_na & !is.na(d))) max_rel <- max(max_rel, rel[!both_na & !is.na(d)])
      if (!isTRUE(all(a == b | both_na, na.rm = FALSE))) exact <- FALSE
      if (!all(ok)) bad_cols <- c(bad_cols, col)
    } else {
      ok <- identical(as.character(a), as.character(b))
      if (!ok) bad_cols <- c(bad_cols, col)
    }
  }
  n_mismatch <- if (length(bad_cols)) length(bad_cols) else 0L
  status <- if (n_matched < nrow(ref) || nrow(new) != nrow(ref) || length(bad_cols) > 0) {
    "different"
  } else if (rounded_used) {
    "different_explained"
  } else if (exact) {
    "identical"
  } else {
    "equal_within_tolerance"
  }
  note <- if (status == "different_explained") {
    rules <- vapply(round_ref, function(r) paste(names(r), unlist(r)), "")
    paste("reference rounded:", paste(names(round_ref), rules, collapse = "; "))
  } else if (n_matched < nrow(ref)) {
    paste0(nrow(ref) - n_matched, " reference row(s) without a match")
  } else if (nrow(new) != nrow(ref)) {
    paste0("row count differs: ref ", nrow(ref), ", new ", nrow(new))
  } else {
    ""
  }
  list(
    status = status, max_rel_diff = max_rel, n_ref = nrow(ref), n_new = nrow(new), n_matched = n_matched,
    n_mismatch = n_mismatch, mismatched_columns = bad_cols, note = note
  )
}
