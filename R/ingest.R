#' Ingest all configured inputs
#'
#' Runs the adapters for every configured input, assembles the SummarizedExperiment and the imported contrast
#' tables, and collects the facts for the ingest report (mapping losses, duplicates, ignored columns, checks of
#' provider fold-change conventions). No file is written here; see [write_ingest_outputs()].
#'
#' @param cfg Config from [read_config()].
#' @return `list(se, contrast_stats, qc, samples, report, provenance, resources, annotation,
#'   normalization_check)`; absent levels are `NULL`. `normalization_check` (see [check_normalization()]) is
#'   filled when both counts and provider-normalised values are configured.
#' @examples
#' \donttest{
#' # see the quickstart in README.md for a complete example on the shipped fixtures
#' }
#' @export
run_ingest <- function(cfg) {
  inp <- cfg$inputs %||% list()
  samples <- if (!is.null(cfg$samples$sheet)) {
    read_sample_sheet(config_path(cfg, cfg$samples$sheet), cfg$samples$factors %||% list())
  }
  res_paths <- verify_resources(cfg)
  gtf_key <- cfg$annotation$gtf %||% "gtf"
  annotation <- if (gtf_key %in% names(res_paths)) read_gtf_genes(res_paths[[gtf_key]])
  id_type <- cfg$annotation$id_type %||% "ensembl_gene_id"
  report <- list()
  prov <- list()
  qc <- NULL
  counts <- tpm <- fpkm <- normalized <- abundance_features <- NULL
  if (!is.null(inp$counts)) {
    if (inp$counts$adapter == "star_reads_per_gene") {
      r <- read_star_counts(
        per_sample_files(cfg, samples, inp$counts$files, "counts"), inp$counts$strand_column %||% "single",
        inp$counts$on_duplicate %||% "error"
      )
      qc <- rbind(qc, map_qc_samples(r$qc, samples))
      counts <- r$counts
    } else {
      r <- read_matrix_table(
        input_path(cfg, inp$counts$file), inp$counts$id_column, inp$counts$sheet, "counts",
        inp$counts$on_duplicate %||% "error"
      )
      counts <- r$values
    }
    prov$counts <- r$provenance
  }
  if (!is.null(inp$abundance)) {
    if (inp$abundance$adapter == "stringtie_gene_abundance") {
      r <- read_stringtie(
        per_sample_files(cfg, samples, inp$abundance$files, "abundance"), inp$abundance$on_duplicate %||% "error"
      )
      tpm <- r$tpm
      fpkm <- r$fpkm
      abundance_features <- r$features
      report$abundance_name_mismatches <- r$name_mismatches
    } else {
      r <- read_matrix_table(
        input_path(cfg, inp$abundance$file), inp$abundance$id_column, inp$abundance$sheet,
        inp$abundance$kind, inp$abundance$on_duplicate %||% "error"
      )
      if (identical(inp$abundance$kind, "tpm")) tpm <- r$values else fpkm <- r$values
    }
    prov$abundance <- r$provenance
  }
  if (!is.null(inp$normalized)) {
    nz <- inp$normalized
    r <- read_matrix_table(input_path(cfg, nz$file), nz$id_column, nz$sheet, "normalized")
    normalized <- r$values
    prov$normalized <- r$provenance
  }
  se <- NULL
  provenance <- do.call(rbind, unname(prov))
  if (!is.null(counts) || !is.null(tpm) || !is.null(normalized)) {
    b <- build_se(samples,
      counts = counts, tpm = tpm, fpkm = fpkm, normalized = normalized, annotation = annotation,
      id_type = id_type, abundance_features = abundance_features,
      symbol_source = cfg$annotation$symbol_source %||% "annotation", provenance = provenance
    )
    se <- b$se
    report <- c(report, b$report)
  }
  norm_check <- NULL
  if (!is.null(counts) && !is.null(normalized)) {
    nc <- check_normalization(
      SummarizedExperiment::assay(se, "counts"), SummarizedExperiment::assay(se, "normalized")
    )
    norm_check <- nc$table
    report <- c(report, nc$facts)
  }
  ct <- NULL
  for (t in inp$contrast_tables %||% list()) {
    x <- read_de_table(input_path(cfg, t$file), t$id, t$adapter, t$sheet, t$columns, t$df_residual)
    report[[paste0("contrast_table_", t$id)]] <- attr(x, "checks")
    provenance <- rbind(provenance, attr(x, "provenance"))
    if (id_type == "gene_symbol") {
      x$gene_symbol <- x$feature_id
    } else if (!is.null(annotation)) {
      x$gene_symbol <- annotation$gene_symbol[match(x$feature_id, annotation$ensembl_gene_id)]
    }
    ct <- bind_contrast_tables(list(ct, x))
  }
  if (!is.null(ct) && !is.null(annotation) && id_type == "gene_symbol") {
    ids <- unique(ct$feature_id)
    report$contrast_symbols_in_annotation <- sum(ids %in% annotation$gene_symbol)
    report$contrast_symbols_not_in_annotation <- sum(!ids %in% annotation$gene_symbol)
  }
  list(
    se = se, contrast_stats = ct, qc = qc, samples = samples, report = report, provenance = provenance,
    resources = res_paths, annotation = annotation, normalization_check = norm_check
  )
}

# One file per sample with the given data level: `pattern` with {input} replaced by its input_column.
per_sample_files <- function(cfg, samples, pattern, level) {
  ss <- samples_with_level(samples, level)
  paths <- vapply(ss$input_column, function(i) input_path(cfg, gsub("{input}", i, pattern, fixed = TRUE)), "")
  stats::setNames(paths, ss$input_column)
}

# Every input file the config points to, so that a missing one is reported before the run starts.
configured_input_files <- function(cfg, samples) {
  inp <- cfg$inputs %||% list()
  files <- character()
  for (lv in c("counts", "abundance")) {
    x <- inp[[lv]]
    # [[ ]]: `x$file` would partially match `files`.
    if (!is.null(x[["files"]])) files <- c(files, per_sample_files(cfg, samples, x[["files"]], lv))
    if (!is.null(x[["file"]])) files <- c(files, input_path(cfg, x[["file"]]))
  }
  if (!is.null(inp$normalized)) files <- c(files, input_path(cfg, inp$normalized$file))
  for (t in inp$contrast_tables %||% list()) files <- c(files, input_path(cfg, t$file))
  unname(files)
}

map_qc_samples <- function(qc, samples) {
  # Adapters name QC rows by input_column; outputs use sample_id.
  qc$sample_id <- samples$sample_id[match(qc$sample_id, samples$input_column)]
  qc
}

#' Write ingest outputs
#'
#' @param ing Result of [run_ingest()].
#' @param out `character(1)`. Project output folder.
#' @return Character vector of written files, invisibly.
#' @examples
#' \donttest{
#' # called by the workflow wrapper after run_ingest()
#' }
#' @export
write_ingest_outputs <- function(ing, out) {
  dir <- file.path(out, "01_ingest")
  obj <- file.path(dir, "objects")
  dir.create(obj, recursive = TRUE, showWarnings = FALSE)
  written <- character()
  if (!is.null(ing$se)) {
    se <- ing$se
    for (a in SummarizedExperiment::assayNames(se)) {
      m <- SummarizedExperiment::assay(se, a)
      d <- data.frame(feature_id = rownames(m), m, check.names = FALSE)
      written <- c(written, write_table(d, dir, a, list(
        description = paste("Assay", a, "(genes x samples)"),
        units = list(values = a), data_columns = colnames(m)
      ), gzip = TRUE))
    }
    rd <- as.data.frame(SummarizedExperiment::rowData(se))
    meta <- list(description = "Feature annotation (rowData).")
    written <- c(written, write_table(rd, dir, "features", meta, gzip = TRUE))
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    cd[] <- lapply(cd, as.character)
    written <- c(written, write_table(cd, dir, "samples", list(
      description = "Sample metadata (colData).",
      data_columns = names(cd)
    )))
    saveRDS(se, file.path(obj, "se.rds"), version = 3)
  }
  if (!is.null(ing$contrast_stats)) {
    written <- c(written, write_table(ing$contrast_stats, dir, "contrasts",
      list(description = "Imported contrast tables (canonical form)."),
      gzip = TRUE
    ))
    saveRDS(ing$contrast_stats, file.path(obj, "contrast_stats.rds"), version = 3)
  }
  if (!is.null(ing$annotation)) {
    written <- c(written, write_table(ing$annotation, dir, "annotation_genes",
      list(description = "Gene records of the cached GTF."),
      gzip = TRUE
    ))
  }
  if (!is.null(ing$qc) && nrow(ing$qc) > 0) {
    written <- c(written, write_table(ing$qc, dir, "qc_metrics", list(description = "QC metrics (long).")))
    saveRDS(ing$qc, file.path(obj, "qc.rds"), version = 3)
  }
  if (!is.null(ing$normalization_check)) {
    written <- c(written, write_table(ing$normalization_check, dir, "normalization_check", list(
      description = "Provider normalisation vs raw counts: implied and recomputed size factors per sample.",
      units = list(implied_size_factor = "raw / normalised (median over genes)")
    )))
  }
  if (!is.null(ing$provenance)) {
    prov <- ing$provenance
    prov$file <- basename(prov$file)
    written <- c(written, write_table(prov, dir, "input_files", list(description = "Input files with md5.")))
  }
  writeLines(ingest_report_md(ing), file.path(dir, "ingest_report.md"))
  invisible(c(written, file.path(dir, "ingest_report.md")))
}

ingest_report_md <- function(ing) {
  lines <- c("# Ingest report", "")
  if (!is.null(ing$se)) {
    lines <- c(lines, sprintf(
      "- Samples: %d; features: %d; assays: %s.", ncol(ing$se), nrow(ing$se),
      paste(SummarizedExperiment::assayNames(ing$se), collapse = ", ")
    ))
  }
  if (!is.null(ing$contrast_stats)) {
    for (cid in unique(ing$contrast_stats$contrast_id)) {
      x <- ing$contrast_stats[ing$contrast_stats$contrast_id == cid, ]
      lines <- c(lines, sprintf(
        "- Contrast table `%s` (%s): %d features; significance metric: %s.", cid, x$source[[1]], nrow(x),
        x$significance_metric[[1]]
      ))
    }
  }
  lines <- c(lines, "", "## Facts", "")
  for (nm in names(ing$report)) {
    v <- ing$report[[nm]]
    txt <- if (is.list(v)) {
      paste(names(v), unlist(lapply(v, paste, collapse = ", ")), sep = ": ", collapse = "; ")
    } else {
      paste(v, collapse = ", ")
    }
    lines <- c(lines, paste0("- ", nm, ": ", txt))
  }
  lines
}
