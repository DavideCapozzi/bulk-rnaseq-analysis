#' Read gene records from a GTF file
#'
#' Streams the (gzipped) GTF through `awk`, so the full file is never held in memory, and keeps one row per
#' `gene` record. Attribute parsing takes the first `gene_id`, `gene_name` and `gene_biotype` of each record,
#' as Ensembl GTFs are written.
#'
#' @param path `character(1)`. GTF or GTF.gz from the resource cache.
#' @return `data.frame` with `ensembl_gene_id`, `gene_symbol` (`NA` when absent), `biotype`, `chromosome`,
#'   `start`, `end`, `strand`.
#' @examples
#' f <- tempfile(fileext = ".gtf")
#' writeLines(paste("1", "src", "gene", "11", "20", ".", "+", ".",
#'   'gene_id "G1"; gene_name "A"; gene_biotype "protein_coding";',
#'   sep = "\t"
#' ), f)
#' read_gtf_genes(f)
#' @export
read_gtf_genes <- function(path) {
  assert_file(path, "GTF")
  prog <- c(
    paste0(
      "BEGIN{OFS=\"\\t\"; print \"ensembl_gene_id\",\"gene_symbol\",\"biotype\",",
      "\"chromosome\",\"start\",\"end\",\"strand\"}"
    ),
    "$3==\"gene\" {",
    "  id=\"\"; nm=\"\"; bt=\"\";",
    "  if (match($9, /gene_id \"[^\"]+\"/))      id=substr($9, RSTART+9,  RLENGTH-10);",
    "  if (match($9, /gene_name \"[^\"]+\"/))    nm=substr($9, RSTART+11, RLENGTH-12);",
    "  if (match($9, /gene_biotype \"[^\"]+\"/)) bt=substr($9, RSTART+14, RLENGTH-15);",
    "  print id, (nm==\"\"?\"NA\":nm), bt, $1, $4, $5, $7 }"
  )
  x <- stream_awk(path, prog)
  x$chromosome <- as.character(x$chromosome)
  x
}

#' Read exon records from a GTF file
#'
#' @param path `character(1)`. GTF or GTF.gz.
#' @param gene_symbols `character` or `NULL`. Keep only exons of these gene names (`NULL` keeps all).
#' @return `data.frame` with `ensembl_gene_id`, `gene_symbol`, `chromosome`, `start`, `end`, `strand`,
#'   `transcript_id`, `transcript_biotype`.
#' @examples
#' f <- tempfile(fileext = ".gtf")
#' writeLines(paste("1", "src", "exon", "11", "20", ".", "+", ".",
#'   'gene_id "G1"; transcript_id "T1"; gene_name "A"; transcript_biotype "protein_coding";',
#'   sep = "\t"
#' ), f)
#' read_gtf_exons(f)
#' @export
read_gtf_exons <- function(path, gene_symbols = NULL) {
  assert_file(path, "GTF")
  prog <- c(
    paste0(
      "BEGIN{OFS=\"\\t\"; print \"ensembl_gene_id\",\"gene_symbol\",\"chromosome\",\"start\",",
      "\"end\",\"strand\",\"transcript_id\",\"transcript_biotype\"}"
    ),
    "$3==\"exon\" {",
    "  id=\"\"; nm=\"\"; tid=\"\"; tb=\"\";",
    "  if (match($9, /gene_id \"[^\"]+\"/))            id=substr($9, RSTART+9,  RLENGTH-10);",
    "  if (match($9, /gene_name \"[^\"]+\"/))          nm=substr($9, RSTART+11, RLENGTH-12);",
    "  if (match($9, /transcript_id \"[^\"]+\"/))      tid=substr($9, RSTART+15, RLENGTH-16);",
    "  if (match($9, /transcript_biotype \"[^\"]+\"/)) tb=substr($9, RSTART+20, RLENGTH-21);",
    "  print id, (nm==\"\"?\"NA\":nm), $1, $4, $5, $7, tid, tb }"
  )
  x <- stream_awk(path, prog)
  x$chromosome <- as.character(x$chromosome)
  if (!is.null(gene_symbols)) x <- x[x$gene_symbol %in% gene_symbols, , drop = FALSE]
  rownames(x) <- NULL
  x
}

# Runs an awk program over a (possibly gzipped) text file and reads the tab-separated result. awk is a
# documented system dependency; it keeps memory use independent of the file size.
stream_awk <- function(path, prog) {
  awk_file <- tempfile(fileext = ".awk")
  on.exit(unlink(awk_file))
  writeLines(prog, awk_file)
  reader <- if (grepl("\\.gz$", path)) "gzip -dc" else "cat"
  cmd <- sprintf("%s %s | awk -F'\\t' -f %s", reader, shQuote(path), shQuote(awk_file))
  as.data.frame(data.table::fread(
    cmd = cmd, sep = "\t", header = TRUE, na.strings = "NA",
    colClasses = list(character = c(1, 2))
  ))
}
