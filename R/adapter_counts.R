#' Read STAR ReadsPerGene count files
#'
#' Handles both STAR's 4-column `ReadsPerGene.out.tab` (unstranded, forward, reverse; first rows `N_*`) and
#' the 2-column files some providers deliver after choosing the strand. `N_*` rows never enter the count
#' matrix: they are returned as QC metrics, because they are reads that were not assigned to any gene.
#'
#' @param files Named `character`. One file per sample; names are `sample_id`s.
#' @param strand_column `character(1)`. `"unstranded"`, `"forward"` or `"reverse"` for 4-column files;
#'   `"single"` for 2-column files.
#' @param on_duplicate `character(1)`. What to do with duplicated gene IDs: `"error"` (default), `"sum"`, `"first"`.
#' @return A list with `counts` (integer matrix, genes x samples, in file order), `qc` (QC long table) and
#'   `provenance` (`data.frame` of file, md5, adapter, params).
#' @examples
#' f <- tempfile()
#' writeLines(c(
#'   "N_unmapped\t5\t5\t5", "N_multimapping\t2\t2\t2", "N_noFeature\t1\t3\t4",
#'   "N_ambiguous\t0\t0\t0", "g1\t10\t1\t9", "g2\t4\t0\t4"
#' ), f)
#' read_star_counts(c(s1 = f), strand_column = "reverse")$counts
#' @export
read_star_counts <- function(files, strand_column = c("single", "unstranded", "forward", "reverse"),
                             on_duplicate = c("error", "sum", "first")) {
  strand_column <- match.arg(strand_column)
  on_duplicate <- match.arg(on_duplicate)
  if (is.null(names(files)) || anyDuplicated(names(files))) {
    abort("read_star_counts(): `files` must be named by unique sample_id")
  }
  col_index <- c(unstranded = 2L, forward = 3L, reverse = 4L)
  per_file <- lapply(names(files), function(s) {
    f <- assert_file(files[[s]], paste0("STAR counts for sample '", s, "'"))
    x <- data.table::fread(f, header = FALSE, sep = "\t", colClasses = list(character = 1))
    if (ncol(x) == 2) {
      if (strand_column != "single") {
        abort("STAR counts ", f, ": 2 columns (strand already selected); set inputs.counts.strand_column: single")
      }
      val <- x[[2]]
    } else if (ncol(x) == 4) {
      if (strand_column == "single") {
        abort("STAR counts ", f, ": 4 columns; set inputs.counts.strand_column to unstranded, forward or reverse")
      }
      val <- x[[col_index[[strand_column]]]]
    } else {
      abort("STAR counts ", f, ": expected 2 or 4 tab-separated columns, found ", ncol(x))
    }
    if (!is_whole_number(val)) abort("STAR counts ", f, ": non-integer counts")
    is_n <- startsWith(x[[1]], "N_")
    list(
      ids = x[[1]][!is_n], counts = as.integer(val[!is_n]),
      qc = data.frame(
        sample_id = rep(s, sum(is_n)), metric = tolower(x[[1]][is_n]), value = as.numeric(val[is_n]),
        source = rep("star_counts", sum(is_n)), stringsAsFactors = FALSE
      )
    )
  })
  ids <- per_file[[1]]$ids
  for (k in seq_along(per_file)) {
    if (!identical(per_file[[k]]$ids, ids)) {
      if (!setequal(per_file[[k]]$ids, ids)) {
        abort("STAR counts: gene IDs differ between ", files[[1]], " and ", files[[k]])
      }
      abort(
        "STAR counts: gene order differs between ", files[[1]], " and ", files[[k]],
        "; files must come from one annotation"
      )
    }
  }
  m <- vapply(per_file, `[[`, integer(length(ids)), "counts")
  m <- matrix(m, ncol = length(files), dimnames = list(ids, names(files)))
  m <- resolve_duplicates(m, on_duplicate, "STAR counts")
  storage.mode(m) <- "integer"
  list(
    counts = m, qc = validate_qc_table(do.call(rbind, lapply(per_file, `[[`, "qc"))),
    provenance = file_provenance(
      files, "star_reads_per_gene",
      list(strand_column = strand_column, on_duplicate = on_duplicate)
    )
  )
}

#' Read StringTie gene abundance tables
#'
#' StringTie can report the same `Gene ID` on several rows (one per locus); `on_duplicate` decides whether
#' this is an error, or whether TPM/FPKM are summed or the first row kept. Rows are aligned on the first
#' file's gene order; a missing gene name coded as `"-"`, `"."` or empty becomes `NA`.
#'
#' @inheritParams read_star_counts
#' @param files Named `character`. One `gene_abundance` file per sample; names are `sample_id`s.
#' @return A list with `tpm`, `fpkm` (numeric matrices, genes x samples), `features` (`data.frame` with
#'   `feature_id`, `gene_name`, `reference`, `strand`, `start`, `end` from the first file),
#'   `name_mismatches` (number of genes whose name differs between samples) and `provenance`.
#' @examples
#' f <- tempfile()
#' writeLines(c(
#'   "Gene ID\tGene Name\tReference\tStrand\tStart\tEnd\tCoverage\tFPKM\tTPM",
#'   "g1\tA\t1\t+\t1\t100\t1.0\t2.0\t3.0", "g2\t-\t1\t-\t5\t50\t0\t0\t0"
#' ), f)
#' read_stringtie(c(s1 = f))$tpm
#' @export
read_stringtie <- function(files, on_duplicate = c("error", "sum", "first")) {
  on_duplicate <- match.arg(on_duplicate)
  need <- c("Gene ID", "Gene Name", "Reference", "Strand", "Start", "End", "FPKM", "TPM")
  per_file <- lapply(names(files), function(s) {
    f <- assert_file(files[[s]], paste0("StringTie abundance for sample '", s, "'"))
    x <- as.data.frame(data.table::fread(f,
      sep = "\t", header = TRUE, check.names = FALSE,
      colClasses = list(character = c("Gene ID", "Gene Name", "Reference"))
    ))
    assert_columns(x, need, paste0("StringTie file ", f))
    x[["Gene Name"]][x[["Gene Name"]] %in% c("-", ".", "")] <- NA_character_
    dup <- duplicated(x[["Gene ID"]])
    if (any(dup)) {
      if (on_duplicate == "error") {
        abort(
          "StringTie file ", f, ": duplicated Gene ID(s) ",
          paste(utils::head(unique(x[["Gene ID"]][dup]), 5), collapse = ", "),
          "; set inputs.abundance.on_duplicate to sum or first"
        )
      }
      vals <- as.matrix(x[, c("TPM", "FPKM")])
      first <- x[!dup, , drop = FALSE]
      if (on_duplicate == "sum") {
        summed <- rowsum(vals, x[["Gene ID"]], reorder = FALSE)
        first$TPM <- summed[first[["Gene ID"]], "TPM"]
        first$FPKM <- summed[first[["Gene ID"]], "FPKM"]
      }
      x <- first
    }
    x
  })
  ids <- per_file[[1]][["Gene ID"]]
  aligned <- lapply(seq_along(per_file), function(k) {
    x <- per_file[[k]]
    if (!setequal(x[["Gene ID"]], ids)) abort("StringTie: gene IDs differ between ", files[[1]], " and ", files[[k]])
    x[match(ids, x[["Gene ID"]]), , drop = FALSE]
  })
  names_mat <- vapply(aligned, function(x) x[["Gene Name"]], character(length(ids)))
  mismatch <- sum(apply(matrix(names_mat, nrow = length(ids)), 1, function(v) length(unique(v)) > 1))
  mk <- function(col) {
    matrix(vapply(aligned, function(x) as.numeric(x[[col]]), numeric(length(ids))),
      ncol = length(files),
      dimnames = list(ids, names(files))
    )
  }
  first <- aligned[[1]]
  list(
    tpm = mk("TPM"), fpkm = mk("FPKM"),
    features = data.frame(
      feature_id = ids, gene_name = first[["Gene Name"]], reference = first[["Reference"]],
      strand = first[["Strand"]], start = as.integer(first[["Start"]]),
      end = as.integer(first[["End"]]), stringsAsFactors = FALSE
    ),
    name_mismatches = mismatch,
    provenance = file_provenance(files, "stringtie_gene_abundance", list(on_duplicate = on_duplicate))
  )
}

#' Read a gene x sample matrix from TSV, CSV or Excel
#'
#' Used for counts, abundances and provider-normalised values delivered as one table. Excel sheets are read
#' with `readxl`, which reads every row even when a workbook declares a wrong used range (`<dimension
#' ref="A1"/>`).
#'
#' @param file `character(1)`. `.tsv`, `.txt`, `.csv`, `.xlsx` or `.xls`.
#' @param id_column `character(1)`. Column with gene identifiers.
#' @param sheet `character(1)` or `NULL`. Excel sheet name (first sheet when `NULL`).
#' @param kind `character(1)`. `"counts"` (integers required), `"tpm"`, `"fpkm"` or `"normalized"`.
#' @param on_duplicate See [read_star_counts()].
#' @return A list with `values` (numeric or integer matrix, columns named as in the file), and `provenance`.
#' @examples
#' f <- tempfile(fileext = ".tsv")
#' writeLines(c("id\tA\tB", "g1\t1\t2", "g2\t3\t4"), f)
#' read_matrix_table(f, "id", kind = "counts")$values
#' @export
read_matrix_table <- function(file, id_column, sheet = NULL, kind = c("counts", "tpm", "fpkm", "normalized"),
                              on_duplicate = c("error", "sum", "first")) {
  kind <- match.arg(kind)
  on_duplicate <- match.arg(on_duplicate)
  x <- read_any_table(file, sheet)
  assert_columns(x, id_column, paste0("matrix table ", file))
  value_cols <- setdiff(names(x), id_column)
  non_num <- value_cols[!vapply(x[value_cols], is.numeric, logical(1))]
  if (length(non_num) > 0) abort("matrix table ", file, ": non-numeric column(s) ", paste(non_num, collapse = ", "))
  m <- as.matrix(x[value_cols])
  rownames(m) <- as.character(x[[id_column]])
  if (anyNA(rownames(m))) abort("matrix table ", file, ": missing values in id column '", id_column, "'")
  if (kind == "counts") {
    if (!is_whole_number(m)) {
      abort(
        "matrix table ", file,
        ": kind 'counts' needs integers (normalised values are never converted back to counts)"
      )
    }
    storage.mode(m) <- "integer"
  }
  m <- resolve_duplicates(m, on_duplicate, paste0("matrix table ", file))
  params <- list(id_column = id_column, sheet = sheet %||% NA, kind = kind)
  list(values = m, provenance = file_provenance(file, "matrix", params))
}

read_any_table <- function(file, sheet = NULL) {
  assert_file(file, "table")
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("xlsx", "xls")) {
    sheets <- readxl::excel_sheets(file)
    if (!is.null(sheet) && !sheet %in% sheets) {
      abort("Excel file ", file, ": sheet '", sheet, "' not found; sheets: ", paste(sheets, collapse = " | "))
    }
    return(as.data.frame(readxl::read_excel(file, sheet = sheet %||% sheets[[1]], guess_max = 1e6)))
  }
  sep <- if (ext == "csv") "," else "\t"
  as.data.frame(data.table::fread(file, sep = sep, header = TRUE, check.names = FALSE))
}

resolve_duplicates <- function(m, how, what) {
  ids <- rownames(m)
  dup <- duplicated(ids)
  if (!any(dup)) {
    return(m)
  }
  if (how == "error") {
    abort(
      what, ": duplicated gene ID(s) ", paste(utils::head(unique(ids[dup]), 5), collapse = ", "),
      " (", sum(dup), " rows); set on_duplicate to sum or first"
    )
  }
  if (how == "first") {
    return(m[!dup, , drop = FALSE])
  }
  out <- rowsum(m, ids, reorder = FALSE)
  storage.mode(out) <- storage.mode(m)
  out
}

file_provenance <- function(files, adapter, params) {
  data.frame(
    file = unname(files), md5 = unname(vapply(files, md5_file, "")), adapter = adapter,
    params = as.character(jsonlite::toJSON(params, auto_unbox = TRUE)), stringsAsFactors = FALSE
  )
}
