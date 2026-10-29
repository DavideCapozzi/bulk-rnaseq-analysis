#' Read MultiQC reports
#'
#' Decodes the plot data embedded in MultiQC HTML reports with `inst/python/multiqc_decode.py` (Python standard
#' library: no R package decodes LZString) and derives per-sample QC metrics from the sections it finds:
#' STAR (`star_input_reads`, `pct_unique`, ..., `assigned_reverse`, `strand_ratio_star`), RSeQC
#' (`rseqc_*_pct`, `pct_cds`, ..., `genebody_3to5`) and, only in reports marked with `read`, FastQC on raw reads
#' (`mean_q_<read>`, `adapter_pct_last_pos_<read>`, ...; `raw_read_pairs` from the General Statistics of `r1`).
#'
#' A MultiQC sample name belongs to the sample whose `input_column` it contains, or to the one named in
#' `aliases`; a name matching no sample or several samples is an error.
#'
#' @param reports `list` of `list(id, file, read)`: report id, path of the HTML report, and `"r1"`/`"r2"` for a
#'   report of FastQC on raw reads (`NULL` otherwise: its FastQC sections are not read).
#' @param samples Sample sheet from [read_sample_sheet()].
#' @param aliases Named `list`: MultiQC sample name -> `input_column`.
#' @return `list(qc, files, adapter_positions, provenance)`: `qc` is a QC long table (`source` =
#'   `multiqc:<id>`), `files` the decoded files (written to `01_ingest/multiqc/` by [write_ingest_outputs()]),
#'   `adapter_positions` the read position at which adapter content was read, per report.
#' @examples
#' \donttest{
#' # needs python3 and a MultiQC HTML report
#' }
#' @export
read_multiqc <- function(reports, samples, aliases = list()) {
  python <- Sys.which("python3")
  if (!nzchar(python)) abort("MultiQC reports: python3 not found on PATH (it decodes the reports)")
  files <- vapply(reports, function(r) assert_file(r$file, paste0("MultiQC report '", r$id, "'")), "")
  ids <- vapply(reports, `[[`, "", "id")
  dir <- tempfile("multiqc_")
  script <- system.file("python", "multiqc_decode.py", package = "bulkrnaseq", mustWork = TRUE)
  out <- suppressWarnings(system2(python, c("-I", shQuote(script), shQuote(dir), shQuote(paste0(ids, "=", files))),
    stdout = TRUE, stderr = TRUE
  ))
  if (!is.null(attr(out, "status"))) abort("MultiQC decoding failed:\n", paste(out, collapse = "\n"))
  map <- function(nm) multiqc_sample_ids(nm, samples, aliases)
  metrics <- list()
  positions <- list()
  for (r in reports) {
    plots <- jsonlite::read_json(file.path(dir, paste0(r$id, "_plotdata.json")))
    src <- paste0("multiqc:", r$id)
    m <- c(multiqc_star(plots, map), multiqc_rseqc(plots, map))
    if (!is.null(r$read)) {
      fq <- multiqc_fastqc(plots, map, r$read)
      positions[[r$id]] <- attr(fq, "adapter_position")
      m <- c(m, fq)
      gs <- utils::read.delim(file.path(dir, paste0(r$id, "_general_stats.tsv")),
        colClasses = "character", check.names = FALSE
      )
      if (r$read == "r1" && "Reads" %in% names(gs)) {
        m <- c(list(raw_read_pairs = stats::setNames(trunc(as.numeric(gs$Reads)), map(gs[["Sample Name"]]))), m)
      }
    }
    metrics <- c(metrics, lapply(m, function(v) list(value = v, source = src)))
  }
  dup <- unique(names(metrics)[duplicated(names(metrics))])
  if (length(dup) > 0) abort("MultiQC reports: metric(s) found in more than one report: ", paste(dup, collapse = ", "))
  if (!is.null(metrics$raw_read_pairs) && !is.null(metrics$star_input_reads)) {
    s <- metrics$star_input_reads
    raw <- metrics$raw_read_pairs$value[names(s$value)]
    metrics <- append(metrics, list(pct_raw_to_star = list(value = 100 * s$value / raw, source = s$source)),
      after = match("star_input_reads", names(metrics))
    )
  }
  qc <- do.call(rbind, lapply(names(metrics), function(k) {
    v <- metrics[[k]]$value
    v <- v[order(match(names(v), samples$sample_id))]
    data.frame(sample_id = names(v), metric = k, value = unname(v), source = metrics[[k]]$source)
  }))
  rownames(qc) <- NULL
  list(
    qc = validate_qc_table(qc), adapter_positions = positions,
    files = file.path(dir, c(paste0(rep(ids, each = 2), c("_plotdata.json", "_general_stats.tsv")), "plot_index.tsv")),
    provenance = file_provenance(files, "multiqc", list(ids = ids))
  )
}

# sample_id of each MultiQC sample name: an alias, else the one input_column the name contains.
multiqc_sample_ids <- function(names, samples, aliases) {
  ic <- samples$input_column
  vapply(names, function(n) {
    hit <- if (!is.null(aliases[[n]])) aliases[[n]] else ic[vapply(ic, grepl, logical(1), x = n, fixed = TRUE)]
    if (length(hit) != 1 || !hit %in% ic) {
      abort(
        "MultiQC sample '", n, "' matches ", if (length(hit) == 0) "no" else paste(hit, collapse = ", "),
        " input_column; add it to inputs.qc.sample_aliases"
      )
    }
    samples$sample_id[match(hit, ic)]
  }, "", USE.NAMES = FALSE)
}

# Bar plot dataset -> list(category -> values named by sample_id).
multiqc_bar <- function(plot, map, ds = 1) {
  ids <- map(unlist(plot$samples[[ds]]))
  cats <- plot$datasets[[ds]]
  stats::setNames(
    lapply(cats, function(k) stats::setNames(vapply(k$data, as.numeric, 0), ids)),
    vapply(cats, `[[`, "", "name")
  )
}

# Line plot dataset -> list of series (x, y), named by sample_id (series name before " - ").
multiqc_xy <- function(plot, map, ds = 1) {
  series <- plot$datasets[[ds]]
  nm <- vapply(series, function(s) strsplit(s$name, " - ", fixed = TRUE)[[1]][[1]], "")
  stats::setNames(lapply(series, function(s) {
    list(
      x = vapply(s$data, function(p) as.numeric(p[[1]]), 0),
      y = vapply(s$data, function(p) as.numeric(if (length(p) > 1) p[[2]] else p[[1]]), 0)
    )
  }), map(nm))
}

need_categories <- function(b, cats, plot) {
  miss <- setdiff(cats, names(b))
  if (length(miss) > 0) abort("MultiQC ", plot, ": missing categories ", paste(miss, collapse = ", "))
}

multiqc_star <- function(plots, map) {
  out <- list()
  p <- plots$star_alignment_plot
  if (is.null(p)) {
    return(out)
  }
  b <- multiqc_bar(p, map)
  cats <- c(
    "Uniquely mapped" = "unique", "Mapped to multiple loci" = "multi", "Mapped to too many loci" = "toomany",
    "Unmapped: too short" = "too_short", "Unmapped: other" = "unmapped_other"
  )
  need_categories(b, names(cats), "star_alignment_plot")
  tot <- Reduce(`+`, b[names(cats)])
  out$star_input_reads <- tot
  for (k in names(cats)) out[[paste0("pct_", cats[[k]])]] <- 100 * b[[k]] / tot
  g <- plots$star_gene_counts
  if (!is.null(g)) {
    labels <- vapply(g$config$data_labels, function(l) if (is.list(l)) l$name %||% "" else l, "")
    strand <- match(c("Reverse Stranded", "Same Stranded"), labels)
    if (anyNA(strand)) abort("MultiQC star_gene_counts: no 'Reverse Stranded' / 'Same Stranded' dataset")
    rev <- multiqc_bar(g, map, strand[[1]])[["Overlapping Genes"]]
    same <- multiqc_bar(g, map, strand[[2]])[["Overlapping Genes"]][names(rev)]
    out$assigned_reverse <- rev
    out$assigned_same <- same
    out$pct_assigned_input <- 100 * rev / tot[names(rev)]
    out$pct_assigned_unique <- 100 * rev / b[["Uniquely mapped"]][names(rev)]
    out$strand_ratio_star <- same / rev
  }
  out
}

multiqc_rseqc <- function(plots, map) {
  out <- list()
  if (!is.null(p <- plots$rseqc_infer_experiment_plot)) {
    b <- multiqc_bar(p, map)
    need_categories(b, c("Antisense", "Sense", "Undetermined"), "rseqc_infer_experiment_plot")
    out$rseqc_antisense_pct <- b$Antisense
    out$rseqc_sense_pct <- b$Sense[names(b$Antisense)]
    out$rseqc_undetermined_pct <- b$Undetermined[names(b$Antisense)]
  }
  if (!is.null(p <- plots$rseqc_read_distribution_plot)) {
    b <- multiqc_bar(p, map)
    need_categories(b, c("CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns"), "rseqc_read_distribution_plot")
    tot <- Reduce(`+`, b)
    out$pct_cds <- 100 * b$CDS_Exons / tot
    out$pct_utr <- 100 * (b[["5'UTR_Exons"]] + b[["3'UTR_Exons"]]) / tot
    out$pct_intron <- 100 * b$Introns / tot
    flank <- names(b)[startsWith(names(b), "TSS") | startsWith(names(b), "TES") | startsWith(names(b), "Other")]
    out$pct_intergenic_or_flank <- 100 * Reduce(`+`, b[flank]) / tot
  }
  if (!is.null(p <- plots$rseqc_gene_body_coverage_plot)) {
    # Mean coverage of percentiles 71-90 over 11-30 of the gene body.
    out$genebody_3to5 <- vapply(multiqc_xy(p, map), function(s) {
      y <- s$y[match(c(71:90, 11:30), as.integer(s$x))]
      if (anyNA(y)) abort("MultiQC rseqc_gene_body_coverage_plot: percentiles 11-30 or 71-90 missing")
      mean(y[1:20]) / mean(y[21:40])
    }, 0)
  }
  out
}

multiqc_fastqc <- function(plots, map, read) {
  out <- list()
  pos <- NULL
  if (!is.null(p <- plots$fastqc_adapter_content_plot)) {
    # Adapter content at the last position of the plot (highest adapter when several are plotted).
    xy <- multiqc_xy(p, map)
    pos <- unique(vapply(xy, function(s) max(s$x), 0))
    v <- vapply(xy, function(s) s$y[which.max(s$x)], 0)
    out[[paste0("adapter_pct_last_pos_", read)]] <- tapply(v, factor(names(v), unique(names(v))), max)
  }
  if (!is.null(p <- plots$fastqc_per_sequence_quality_scores_plot)) {
    out[[paste0("mean_q_", read)]] <- vapply(multiqc_xy(p, map), function(s) sum(s$x * s$y) / sum(s$y), 0)
  }
  if (!is.null(p <- plots$fastqc_per_base_sequence_quality_plot)) {
    xy <- multiqc_xy(p, map)
    out[[paste0("first_pos_meanq_lt20_", read)]] <- vapply(xy, function(s) {
      if (any(s$y < 20)) trunc(s$x[s$y < 20][[1]]) else NA_real_
    }, 0)
    out[[paste0("meanq_pos30plus_", read)]] <- vapply(xy, function(s) mean(s$y[s$x >= 30]), 0)
  }
  if (!is.null(p <- plots$fastqc_sequence_duplication_levels_plot)) {
    i100 <- match(">100", unlist(p$config$categories))
    if (is.na(i100)) abort("MultiQC fastqc_sequence_duplication_levels_plot: no '>100' category")
    xy <- multiqc_xy(p, map)
    out[[paste0("dup_lvl1_pct_", read)]] <- vapply(xy, function(s) s$y[[1]], 0)
    out[[paste0("dup_ge100_pct_", read)]] <- vapply(xy, function(s) sum(s$y[i100:length(s$y)]), 0)
  }
  if (!is.null(p <- plots$fastqc_per_sequence_gc_content_plot)) {
    out[[paste0("mean_gc_", read)]] <- vapply(multiqc_xy(p, map), function(s) sum(s$x * s$y) / sum(s$y), 0)
  }
  out <- lapply(out, function(v) stats::setNames(as.numeric(v), names(v)))
  attr(out, "adapter_position") <- pos
  out
}
