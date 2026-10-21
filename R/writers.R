#' Plot theme of the package
#'
#' @param base_size `numeric(1)`. Base font size (points).
#' @return A `ggplot2` theme.
#' @examples
#' theme_bulkrnaseq()
#' @export
theme_bulkrnaseq <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e6e5e1", linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

#' Write a table with its metadata sidecar
#'
#' Writes `<dir>/<name>.tsv` (or `.tsv.gz` when `gzip = TRUE`) and `<dir>/<name>.meta.json`. Numbers are
#' written with 15 significant digits; `NA` is written as an empty field (data.table convention). Gzip output
#' carries no timestamp, so two runs produce identical files.
#'
#' @param x `data.frame`.
#' @param dir `character(1)`. Output folder (created if needed).
#' @param name `character(1)`. Table name without extension (snake_case).
#' @param meta `list`. Metadata (description, units, filters, seed, ...); the writer adds `columns`, `n_rows`,
#'   `function` and `package_version`.
#' @param gzip `logical(1)`. Compress the TSV.
#' @return Path of the TSV, invisibly.
#' @examples
#' d <- tempfile()
#' write_table(data.frame(a = 1:2, b = c("x", "y")), d, "demo", list(description = "example"))
#' list.files(d)
#' @export
write_table <- function(x, dir, name, meta = list(), gzip = FALSE) {
  if (!grepl("^[a-z0-9_]+$", name)) abort("write_table(): table name '", name, "' must be snake_case")
  bad <- names(x)[!grepl("^[a-z0-9_]+$", names(x)) & !names(x) %in% (meta$data_columns %||% character())]
  if (length(bad) > 0) {
    abort("write_table(): column(s) of '", name, "' not in snake_case: ", paste(bad, collapse = ", "))
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, if (gzip) ".tsv.gz" else ".tsv"))
  if (gzip) {
    tmp <- tempfile(fileext = ".tsv")
    data.table::fwrite(x, tmp, sep = "\t")
    gzip_deterministic(tmp, path)
    unlink(tmp)
  } else {
    data.table::fwrite(x, path, sep = "\t")
  }
  meta$table <- name
  meta$n_rows <- nrow(x)
  meta$columns <- names(x)
  meta$package_version <- as.character(utils::packageVersion("bulkrnaseq"))
  write_json(meta, file.path(dir, paste0(name, ".meta.json")))
  invisible(path)
}

# gzip without file name and mtime in the header, so the bytes depend only on the content.
gzip_deterministic <- function(src, dest) {
  status <- system2("gzip", c("-n", "-c", shQuote(src)), stdout = dest)
  if (!identical(as.integer(status), 0L)) abort("gzip failed for ", src)
  invisible(dest)
}

write_json <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA), path)
}

#' Write a figure deterministically
#'
#' PNG with fixed size, resolution and white background; label repulsion in package plots is seeded, so two
#' runs give identical files.
#'
#' @param plot A `ggplot` object.
#' @param path `character(1)`. Output `.png` path.
#' @param width,height `numeric(1)`. Size in inches.
#' @param dpi `numeric(1)`. Resolution.
#' @return `path`, invisibly.
#' @examples
#' p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
#'   ggplot2::geom_point()
#' write_figure(p, tempfile(fileext = ".png"))
#' @export
write_figure <- function(plot, path, width = 7, height = 5, dpi = 150) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi, bg = "white", device = "png")
  invisible(path)
}

#' Write an Excel workbook deterministically
#'
#' One sheet per table plus a first sheet `Info` (parameter, value). Document dates and zip timestamps are set
#' to a fixed value so that two runs produce byte-identical files.
#'
#' @param tables Named `list` of `data.frame`s (sheet names are truncated to 31 characters).
#' @param info `data.frame` with columns `parameter` and `value`.
#' @param path `character(1)`. Output `.xlsx`.
#' @return `path`, invisibly.
#' @examples
#' write_xlsx_deterministic(
#'   list(a = data.frame(x = 1)), data.frame(parameter = "p", value = "v"),
#'   tempfile(fileext = ".xlsx")
#' )
#' @export
write_xlsx_deterministic <- function(tables, info, path) {
  wb <- openxlsx::createWorkbook(creator = "bulkrnaseq")
  sheets <- c(list(Info = info), tables)
  for (nm in names(sheets)) {
    sheet <- substr(nm, 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, sheets[[nm]])
  }
  tmp <- tempfile(fileext = ".xlsx")
  openxlsx::saveWorkbook(wb, tmp, overwrite = TRUE)
  normalise_zip(tmp, path)
  unlink(tmp)
  invisible(path)
}

# Rewrites an Office zip with fixed document dates and entry timestamps, entries in sorted order.
normalise_zip <- function(src, dest, fixed_time = "2000-01-01T00:00:00Z") {
  ex <- tempfile("xlsx_")
  dir.create(ex)
  on.exit(unlink(ex, recursive = TRUE), add = TRUE)
  utils::unzip(src, exdir = ex)
  core <- file.path(ex, "docProps", "core.xml")
  if (file.exists(core)) {
    txt <- readLines(core, warn = FALSE, encoding = "UTF-8")
    stamp_re <- "(<dcterms:(created|modified)[^>]*>)[^<]*(</dcterms:(created|modified)>)"
    txt <- gsub(stamp_re, paste0("\\1", fixed_time, "\\3"), txt)
    writeLines(txt, core, useBytes = TRUE)
  }
  files <- sort(list.files(ex, recursive = TRUE, all.files = TRUE))
  stamp <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
  for (f in files) Sys.setFileTime(file.path(ex, f), stamp)
  if (file.exists(dest)) unlink(dest)
  # Absolute destination: the zip is written from inside the extraction folder.
  dest_abs <- file.path(normalizePath(dirname(dest)), basename(dest))
  withr::with_dir(ex, zip::zip(dest_abs, files = files, mode = "mirror"))
  invisible(dest)
}

#' Write a module's outputs
#'
#' Writes every table (with meta), every figure, and optionally the module's Excel workbook. An
#' `inference_scope` column of a table is copied into the workbook's `Info` sheet.
#'
#' @param result `list(tables, figures, meta)` returned by a module.
#' @param dir `character(1)`. Module folder (`NN_<module>`).
#' @param params `list`. Module parameters (written to the `Info` sheet).
#' @param xlsx_path `character(1)` or `NULL`. Excel output, `NULL` to skip.
#' @param gzip_tables `character`. Tables written as `.tsv.gz`.
#' @param figure_sizes Named `list` of `c(width, height)` per figure.
#' @return Character vector of written files, invisibly.
#' @examples
#' res <- list(tables = list(t1 = data.frame(a = 1)), figures = list(), meta = list())
#' write_module_outputs(res, tempfile(), list(p = 1))
#' @export
write_module_outputs <- function(result, dir, params = list(), xlsx_path = NULL, gzip_tables = character(),
                                 figure_sizes = list()) {
  written <- character()
  for (nm in names(result$tables)) {
    written <- c(written, write_table(result$tables[[nm]], file.path(dir, "tables"), nm, result$meta[[nm]] %||% list(),
      gzip = nm %in% gzip_tables
    ))
  }
  for (nm in names(result$figures)) {
    size <- figure_sizes[[nm]] %||% c(7, 5)
    png_path <- file.path(dir, "figures", paste0(nm, ".png"))
    written <- c(written, write_figure(result$figures[[nm]], png_path, size[[1]], size[[2]]))
  }
  if (!is.null(xlsx_path)) {
    scope_of <- function(t) if ("inference_scope" %in% names(t)) unique(stats::na.omit(t$inference_scope))
    scopes <- unique(unlist(lapply(result$tables, scope_of)))
    flat <- unlist(params)
    info <- data.frame(
      parameter = c(names(flat), "inference_scope"),
      value = c(as.character(flat), if (length(scopes)) paste(scopes, collapse = ", ") else "no p-values"),
      stringsAsFactors = FALSE
    )
    dir.create(dirname(xlsx_path), recursive = TRUE, showWarnings = FALSE)
    written <- c(written, write_xlsx_deterministic(result$tables, info, xlsx_path))
  }
  invisible(written)
}
