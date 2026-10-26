# Internal helpers shared by adapters, modules and writers.

`%||%` <- function(x, y) if (is.null(x)) y else x

# Errors name the object and the field at fault, so a user can fix the config or the file without reading code.
abort <- function(..., call = NULL) {
  stop(paste0(...), call. = FALSE)
}

assert_columns <- function(df, required, what) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    abort(
      what, ": missing column(s) ", paste(sQuote(missing, FALSE), collapse = ", "),
      "; available: ", paste(names(df), collapse = ", ")
    )
  }
  invisible(TRUE)
}

assert_file <- function(path, what) {
  if (!is.character(path) || length(path) != 1 || !file.exists(path)) {
    abort(what, ": file not found: ", path)
  }
  invisible(path)
}

# A path is kept as given when absolute; otherwise it is resolved against `base`, the folder of the config file.
resolve_path <- function(path, base) {
  if (is.null(path)) {
    return(NULL)
  }
  if (grepl("^(/|[A-Za-z]:[/\\\\]|~)", path)) {
    return(path.expand(path))
  }
  normalizePath(file.path(base, path), mustWork = FALSE)
}

md5_file <- function(path) unname(tools::md5sum(path))

is_whole_number <- function(x, tol = 0) {
  x <- x[!is.na(x)]
  all(abs(x - round(x)) <= tol)
}

package_file <- function(...) {
  path <- system.file(..., package = "bulkrnaseq")
  if (!nzchar(path)) abort("bulkrnaseq installation is incomplete: ", file.path(...), " not found")
  path
}
