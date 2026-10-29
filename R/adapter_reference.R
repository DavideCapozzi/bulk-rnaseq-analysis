#' Read an MSigDB gene set export
#'
#' Reads the long table written from `msigdbr` (one row per set and gene): columns `gs_collection`,
#' `gs_subcollection`, `gs_name`, `ensembl_gene`, `gene_symbol`, `db_version`.
#'
#' @param path `character(1)`. TSV, optionally gzipped.
#' @param sets `character` or `NULL`. Set names to keep; every one must be present.
#' @return `data.frame` with the columns above (all character).
#' @examples
#' f <- tempfile(fileext = ".tsv")
#' writeLines(c(
#'   "gs_collection\tgs_subcollection\tgs_name\tensembl_gene\tgene_symbol\tdb_version",
#'   "C8\t\tSET_A\tENSG00000000003\tTSPAN6\tv1"
#' ), f)
#' read_gene_sets(f)
#' @export
read_gene_sets <- function(path, sets = NULL) {
  assert_file(path, "gene set export")
  cols <- c("gs_collection", "gs_subcollection", "gs_name", "ensembl_gene", "gene_symbol", "db_version")
  x <- as.data.frame(data.table::fread(path, sep = "\t", colClasses = "character", select = cols, na.strings = ""))
  assert_columns(x, cols, paste0("gene set export ", path))
  if (!is.null(sets)) {
    missing <- setdiff(sets, x$gs_name)
    if (length(missing) > 0) abort("gene set export ", basename(path), ": set(s) not found: ", toString(missing))
    x <- x[x$gs_name %in% sets, , drop = FALSE]
  }
  x
}

#' Read reference expression profiles of the Human Protein Atlas
#'
#' Reads an HPA download (zipped or plain TSV): `rna_immune_cell_monaco.tsv` (`kind = "hpa_immune_cell"`,
#' column `Immune cell`, values `nTPM`) or `rna_single_cell_type.tsv` (`kind = "hpa_single_cell"`, column
#' `Cell type`, values `nCPM`). Only the requested cell types are kept; a missing one is an error, so that a
#' renamed label in a new HPA release stops the run instead of producing empty columns.
#'
#' @param path `character(1)`. `.tsv.zip` or `.tsv`.
#' @param kind `character(1)`. `"hpa_immune_cell"` or `"hpa_single_cell"`.
#' @param cells `character`. Cell types to keep.
#' @return `data.frame` with `feature_id` (the HPA `Gene` column, an Ensembl gene ID), `cell`, `value`.
#' @examples
#' f <- tempfile(fileext = ".tsv")
#' writeLines(c("Gene\tGene name\tImmune cell\tTPM\tpTPM\tnTPM", "ENSG1\tA\tT-reg\t1\t1\t2.5"), f)
#' read_reference_profiles(f, "hpa_immune_cell", "T-reg")
#' @export
read_reference_profiles <- function(path, kind = c("hpa_immune_cell", "hpa_single_cell"), cells) {
  kind <- match.arg(kind)
  assert_file(path, "reference profiles")
  cell_col <- c(hpa_immune_cell = "Immune cell", hpa_single_cell = "Cell type")[[kind]]
  value_col <- c(hpa_immune_cell = "nTPM", hpa_single_cell = "nCPM")[[kind]]
  file <- path
  if (grepl("\\.zip$", path)) {
    dir <- withr::local_tempdir()
    inner <- utils::unzip(path, list = TRUE)$Name
    inner <- inner[grepl("\\.tsv$", inner)]
    if (length(inner) != 1) abort("reference profiles ", basename(path), ": expected one .tsv inside the zip")
    file <- utils::unzip(path, files = inner, exdir = dir)
  }
  x <- data.table::fread(
    file,
    sep = "\t", select = c("Gene", cell_col, value_col), colClasses = list(character = c("Gene", cell_col))
  )
  assert_columns(x, c("Gene", cell_col, value_col), paste0("reference profiles ", basename(path)))
  missing <- setdiff(cells, x[[cell_col]])
  if (length(missing) > 0) {
    abort(
      "reference profiles ", basename(path), ": cell type(s) not found: ", toString(missing),
      " (available: ", toString(sort(unique(x[[cell_col]]))), ")"
    )
  }
  x <- x[x[[cell_col]] %in% cells]
  out <- data.frame(feature_id = x$Gene, cell = x[[cell_col]], value = as.numeric(x[[value_col]]))
  if (anyDuplicated(out[c("feature_id", "cell")])) {
    abort("reference profiles ", basename(path), ": a gene appears more than once for the same cell type")
  }
  out
}
