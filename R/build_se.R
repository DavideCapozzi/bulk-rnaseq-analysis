#' Assemble the sample-level SummarizedExperiment
#'
#' Aligns the available matrices on one feature order, renames columns from `input_column` to `sample_id`
#' through the sample sheet, and annotates features from the cached GTF. Provider-normalised values that cover
#' fewer genes than the counts are kept with `NA` for the missing genes; they are never converted to counts.
#'
#' @param samples Sample sheet from [read_sample_sheet()].
#' @param counts,tpm,fpkm,normalized Matrices with columns named by `input_column` (any may be `NULL`).
#' @param annotation `data.frame` from [read_gtf_genes()] or `NULL`.
#' @param id_type `character(1)`. `"ensembl_gene_id"` or `"gene_symbol"`: what the input IDs are.
#' @param abundance_features `data.frame` from [read_stringtie()]`$features`, used when
#'   `symbol_source = "abundance"`.
#' @param symbol_source `character(1)`. `"annotation"` or `"abundance"` (gene names from StringTie, falling
#'   back to the feature ID when missing).
#' @param provenance `data.frame` of input files (file, md5, adapter, params).
#' @return A list with `se` (a validated `SummarizedExperiment`) and `report` (named list of ingest facts).
#' @examples
#' ss <- data.frame(sample_id = c("s1", "s2"), input_column = c("A", "B"), row.names = c("s1", "s2"))
#' ss$inputs_available <- list(NA_character_, NA_character_)
#' m <- matrix(c(1L, 2L, 3L, 4L), 2, dimnames = list(c("g1", "g2"), c("A", "B")))
#' build_se(ss, counts = m, id_type = "gene_symbol")$se
#' @export
build_se <- function(samples, counts = NULL, tpm = NULL, fpkm = NULL, normalized = NULL, annotation = NULL,
                     id_type = c("ensembl_gene_id", "gene_symbol", "entrez_id"), abundance_features = NULL,
                     symbol_source = c("annotation", "abundance"), provenance = NULL) {
  id_type <- match.arg(id_type)
  symbol_source <- match.arg(symbol_source)
  mats <- Filter(Negate(is.null), list(counts = counts, tpm = tpm, fpkm = fpkm, normalized = normalized))
  if (length(mats) == 0) abort("build_se(): no sample-level matrix provided")
  primary <- names(mats)[[1]]
  report <- list()
  level_of <- c(counts = "counts", tpm = "abundance", fpkm = "abundance", normalized = "normalized")
  ss <- samples_with_level(samples, level_of[[primary]])
  feature_ids <- rownames(mats[[primary]])
  assays <- list()
  for (nm in names(mats)) {
    m <- mats[[nm]]
    missing_cols <- setdiff(ss$input_column, colnames(m))
    if (length(missing_cols) > 0) {
      abort(
        "assay '", nm, "': sample sheet input_column value(s) not found among the input columns: ",
        paste(missing_cols, collapse = ", "), " (available: ", paste(colnames(m), collapse = ", "), ")"
      )
    }
    ignored <- setdiff(colnames(m), ss$input_column)
    if (length(ignored) > 0) report[[paste0("ignored_columns_", nm)]] <- ignored
    m <- m[, ss$input_column, drop = FALSE]
    colnames(m) <- ss$sample_id
    if (nm == primary) {
      assays[[nm]] <- m
      next
    }
    if (nm == "normalized") {
      # Providers often normalise a filtered gene set; the remaining genes stay NA rather than being dropped.
      extra <- setdiff(rownames(m), feature_ids)
      if (length(extra) > 0) {
        abort("assay 'normalized': ", length(extra), " gene(s) absent from the primary assay, e.g. ", extra[[1]])
      }
      full <- matrix(NA_real_, length(feature_ids), ncol(m), dimnames = list(feature_ids, colnames(m)))
      full[rownames(m), ] <- m
      report$normalized_genes <- nrow(m)
      assays[[nm]] <- full
    } else {
      if (!setequal(rownames(m), feature_ids)) {
        abort("assay '", nm, "': gene IDs differ from the primary assay '", primary, "'")
      }
      assays[[nm]] <- m[feature_ids, , drop = FALSE]
    }
  }
  rd <- annotate_features(feature_ids, annotation, id_type, abundance_features, symbol_source)
  report <- c(report, attr(rd, "report"))
  cd <- S4Vectors::DataFrame(ss[, setdiff(names(ss), "inputs_available"), drop = FALSE], check.names = FALSE)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays, rowData = S4Vectors::DataFrame(rd, check.names = FALSE),
    colData = cd
  )
  rownames(se) <- feature_ids
  S4Vectors::metadata(se)$provenance <- provenance
  S4Vectors::metadata(se)$id_type <- id_type
  list(se = validate_se(se), report = report)
}

#' Annotate features from a GTF table
#'
#' For Ensembl IDs the annotation is joined by ID; for gene symbols, a symbol maps to an Ensembl ID only when
#' it is unique in the annotation, and ambiguous or missing symbols are counted in the report.
#' `gene_symbol_unique` appends the feature ID to duplicated symbols so that plots and tables keep one label
#' per row.
#'
#' @inheritParams build_se
#' @param feature_ids `character`. Input gene IDs in assay order.
#' @return `data.frame` (rowData) with attribute `report`.
#' @examples
#' ann <- data.frame(
#'   ensembl_gene_id = "E1", gene_symbol = "A", biotype = "protein_coding",
#'   chromosome = "1", start = 1L, end = 9L, strand = "+"
#' )
#' annotate_features(c("E1", "E2"), ann, "ensembl_gene_id")
#' @export
annotate_features <- function(feature_ids, annotation = NULL, id_type = "ensembl_gene_id", abundance_features = NULL,
                              symbol_source = "annotation") {
  rd <- data.frame(feature_id = feature_ids, stringsAsFactors = FALSE)
  report <- list(n_features = length(feature_ids))
  ann_cols <- c("biotype", "chromosome", "start", "end", "strand")
  if (id_type == "ensembl_gene_id") {
    rd$ensembl_gene_id <- feature_ids
    if (!is.null(annotation)) {
      i <- match(feature_ids, annotation$ensembl_gene_id)
      report$unmapped_to_annotation <- sum(is.na(i))
      ann_symbol <- annotation$gene_symbol[i]
      for (col in ann_cols) rd[[col]] <- annotation[[col]][i]
    } else {
      ann_symbol <- rep(NA_character_, length(feature_ids))
    }
    if (symbol_source == "abundance") {
      if (is.null(abundance_features)) abort("annotation.symbol_source 'abundance' needs an abundance input")
      st_name <- abundance_features$gene_name[match(feature_ids, abundance_features$feature_id)]
      rd$gene_symbol <- ifelse(is.na(st_name), feature_ids, st_name)
    } else {
      rd$gene_symbol <- ifelse(is.na(ann_symbol), feature_ids, ann_symbol)
    }
    rd$symbol_matches_annotation <- (rd$gene_symbol == ann_symbol) | (is.na(ann_symbol) & rd$gene_symbol == feature_ids)
  } else if (id_type == "gene_symbol") {
    rd$gene_symbol <- feature_ids
    if (!is.null(annotation)) {
      sym <- annotation$gene_symbol
      dup_sym <- unique(sym[duplicated(sym) & !is.na(sym)])
      i <- match(feature_ids, sym)
      ambiguous <- feature_ids %in% dup_sym
      i[ambiguous] <- NA_integer_
      rd$ensembl_gene_id <- annotation$ensembl_gene_id[i]
      for (col in ann_cols) rd[[col]] <- annotation[[col]][i]
      report$symbols_mapped <- sum(!is.na(i))
      report$symbols_unmapped <- sum(is.na(i) & !ambiguous)
      report$symbols_ambiguous <- sum(ambiguous)
    }
  } else {
    rd$entrez_id <- feature_ids
  }
  dup <- duplicated(rd$gene_symbol) | duplicated(rd$gene_symbol, fromLast = TRUE)
  rd$symbol_duplicated <- dup
  rd$gene_symbol_unique <- ifelse(dup, paste0(rd$gene_symbol, "|", feature_ids), rd$gene_symbol)
  report$duplicated_symbols <- length(unique(rd$gene_symbol[dup]))
  attr(rd, "report") <- report
  rd
}
