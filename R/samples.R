#' Read a sample sheet
#'
#' Sample metadata come only from the sample sheet, never from file or column names, because provider
#' sample names are often inconsistent. `input_column` links each sample to its file (`{input}` placeholder)
#' or to a column of a matrix input.
#'
#' @param path `character(1)`. Tab-separated file with at least `sample_id` and `input_column`.
#' @param factors `list`. Named list factor -> ordered levels (reference first), from `samples.factors`.
#' @return A `data.frame` with one row per sample, `sample_id` as row names, factor columns as factors and a
#'   list column `inputs_available` (character vectors of data levels; `NA` means all).
#' @examples
#' f <- tempfile(fileext = ".tsv")
#' writeLines(c("sample_id\tinput_column\tcondition", "a\tA\tctrl", "b\tB\ttrt"), f)
#' read_sample_sheet(f, list(condition = c("ctrl", "trt")))
#' @export
read_sample_sheet <- function(path, factors = list()) {
  assert_file(path, "sample sheet")
  ss <- as.data.frame(data.table::fread(path, sep = "\t", colClasses = "character", na.strings = ""))
  assert_columns(ss, c("sample_id", "input_column"), paste0("sample sheet ", path))
  bad <- ss$sample_id[!grepl("^[A-Za-z0-9_.-]+$", ss$sample_id)]
  if (length(bad) > 0) abort("sample sheet: sample_id must match [A-Za-z0-9_.-]+: ", paste(bad, collapse = ", "))
  for (col in c("sample_id", "input_column")) {
    dup <- unique(ss[[col]][duplicated(ss[[col]])])
    if (length(dup) > 0) abort("sample sheet: duplicated ", col, ": ", paste(dup, collapse = ", "))
  }
  for (f in names(factors)) {
    if (!f %in% names(ss)) abort("samples.factors: column '", f, "' is not in the sample sheet")
    lev <- as.character(unlist(factors[[f]]))
    extra <- setdiff(unique(ss[[f]]), lev)
    if (length(extra) > 0) {
      abort("sample sheet column '", f, "': values not in samples.factors: ", paste(extra, collapse = ", "))
    }
    ss[[f]] <- factor(ss[[f]], levels = lev)
  }
  avail <- if ("inputs_available" %in% names(ss)) ss$inputs_available else rep(NA_character_, nrow(ss))
  ss$inputs_available <- lapply(avail, function(a) if (is.na(a)) NA_character_ else trimws(strsplit(a, ",")[[1]]))
  rownames(ss) <- ss$sample_id
  ss
}

# Samples that have a given data level (inputs_available = NA means every configured level).
samples_with_level <- function(samples, level) {
  keep <- vapply(samples$inputs_available, function(a) all(is.na(a)) || level %in% a, logical(1))
  samples[keep, , drop = FALSE]
}
