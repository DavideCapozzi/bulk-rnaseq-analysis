#' Expression-level QC metrics per sample
#'
#' Computed on the full count matrix (no gene filter): library size (reads assigned to genes), numbers of
#' detected genes, and the share of the library taken by mitochondrial, rRNA, ribosomal-protein,
#' protein-coding, lncRNA and hemoglobin genes and by the most expressed genes.
#'
#' @param se `SummarizedExperiment` with a `counts` assay (and optionally `tpm`) and rowData `gene_symbol`,
#'   `biotype`, `chromosome`.
#' @param ribosomal_protein_regex `character(1)`. Regex on gene symbols; only protein-coding genes count.
#' @param hemoglobin_genes `character`. Hemoglobin gene symbols.
#' @param rrna_biotypes `character`. Biotypes counted as rRNA.
#' @param top_n_share `integer(1)`. Number of most expressed genes in `pct_counts_top<N>_genes`.
#' @param organism `character(1)`. Organism key (mitochondrial chromosome name).
#' @return QC long table (`sample_id`, `metric`, `value`, `source = "expression"`).
#' @examples
#' m <- matrix(c(90L, 10L, 0L, 50L, 40L, 10L), 3, dimnames = list(c("a", "b", "c"), c("s1", "s2")))
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   rowData = S4Vectors::DataFrame(
#'     feature_id = rownames(m), gene_symbol = c("A", "RPL3", "HBB"),
#'     biotype = "protein_coding", chromosome = c("MT", "1", "1")
#'   )
#' )
#' expression_qc(se, hemoglobin_genes = "HBB", top_n_share = 1)
#' @export
expression_qc <- function(se, ribosomal_protein_regex = "^RP[SL][0-9]+[A-Z]?[0-9]*$", hemoglobin_genes = character(),
                          rrna_biotypes = c("rRNA", "Mt_rRNA", "rRNA_pseudogene"), top_n_share = 50,
                          organism = "human") {
  org <- organism_info(organism)
  m <- SummarizedExperiment::assay(se, "counts")
  rd <- SummarizedExperiment::rowData(se)
  lib <- colSums(m)
  is_mt <- rd$chromosome %in% org$chromosome_mito
  is_pc <- rd$biotype %in% "protein_coding"
  is_rp <- grepl(ribosomal_protein_regex, rd$gene_symbol) & is_pc
  pct <- function(rows) 100 * colSums(m[rows, , drop = FALSE]) / lib
  top_n <- min(top_n_share, nrow(m))
  metrics <- list(
    library_size = lib,
    genes_count_ge1 = colSums(m >= 1),
    genes_count_ge10 = colSums(m >= 10),
    pc_genes_count_ge10 = colSums(m[is_pc, , drop = FALSE] >= 10),
    pct_counts_mito = pct(is_mt),
    pct_counts_rrna = pct(rd$biotype %in% rrna_biotypes),
    pct_counts_ribo_protein = pct(is_rp),
    pct_counts_protein_coding = pct(is_pc),
    pct_counts_lncrna = pct(rd$biotype %in% "lncRNA"),
    pct_counts_hemoglobin = pct(rd$gene_symbol %in% hemoglobin_genes)
  )
  metrics[[paste0("pct_counts_top", top_n, "_genes")]] <- apply(m, 2, function(x) {
    100 * sum(sort(x, decreasing = TRUE)[seq_len(top_n)]) / sum(x)
  })
  if ("tpm" %in% SummarizedExperiment::assayNames(se)) {
    tpm <- SummarizedExperiment::assay(se, "tpm")
    metrics$genes_tpm_ge1 <- colSums(tpm >= 1)
    metrics$pct_tpm_mito <- colSums(tpm[is_mt, , drop = FALSE]) / 1e4
  }
  out <- do.call(rbind, lapply(names(metrics), function(nm) {
    data.frame(
      sample_id = colnames(m), metric = nm, value = as.numeric(metrics[[nm]]), source = "expression",
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  validate_qc_table(out)
}

#' Most expressed genes per sample
#'
#' @param se `SummarizedExperiment` with a `counts` assay and rowData `gene_symbol`, `biotype`.
#' @param n `integer(1)`. Genes per sample.
#' @return `data.frame` with `sample_id`, `rank`, `feature_id`, `gene_symbol`, `biotype`, `pct_of_library`
#'   (percentage of the sample's assigned reads).
#' @examples
#' m <- matrix(c(90L, 10L, 50L, 40L), 2, dimnames = list(c("a", "b"), c("s1", "s2")))
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   rowData = S4Vectors::DataFrame(feature_id = c("a", "b"), gene_symbol = c("A", "B"), biotype = "lncRNA")
#' )
#' top_genes(se, 1)
#' @export
top_genes <- function(se, n = 15) {
  m <- SummarizedExperiment::assay(se, "counts")
  rd <- SummarizedExperiment::rowData(se)
  lib <- colSums(m)
  n <- min(n, nrow(m))
  out <- do.call(rbind, lapply(colnames(m), function(s) {
    o <- order(m[, s], decreasing = TRUE)[seq_len(n)]
    data.frame(
      sample_id = s, rank = seq_len(n), feature_id = rownames(m)[o], gene_symbol = rd$gene_symbol[o],
      biotype = rd$biotype[o] %||% NA_character_, pct_of_library = 100 * m[o, s] / lib[[s]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

#' QC metrics as one row per sample
#'
#' @param qc QC long table (see [validate_qc_table()]), possibly from several sources.
#' @param samples `character`. Sample order of the output.
#' @return `data.frame` with `sample_id` and one column per metric, in order of first appearance. A metric
#'   reported by two sources gets the source appended (`<metric>_<source>`).
#' @examples
#' qc <- data.frame(sample_id = c("s1", "s2"), metric = "n_unmapped", value = c(3, 4), source = "star_counts")
#' qc_wide(qc, c("s1", "s2"))
#' @export
qc_wide <- function(qc, samples) {
  key <- qc$metric
  clash <- unique(qc$metric[duplicated(paste(qc$sample_id, qc$metric, sep = "\r"))])
  key[qc$metric %in% clash] <- paste0(qc$metric, "_", qc$source)[qc$metric %in% clash]
  out <- data.frame(sample_id = samples, stringsAsFactors = FALSE)
  for (k in unique(key)) out[[k]] <- qc$value[key == k][match(samples, qc$sample_id[key == k])]
  out
}

# Plot label of each sample: the sample sheet column named by samples.label_column, else sample_id.
sample_labels <- function(se, cfg) {
  cd <- SummarizedExperiment::colData(se)
  col <- cfg$samples$label_column
  if (!is.null(col) && col %in% names(cd)) as.character(cd[[col]]) else colnames(se)
}

# Experimental factors declared in samples.factors that are present in colData, in config order.
sample_factors <- function(se, cfg) {
  intersect(names(cfg$samples$factors %||% list()), names(SummarizedExperiment::colData(se)))
}

#' Bar chart of QC metrics per sample
#'
#' @param wide `data.frame` from [qc_wide()].
#' @param metrics `character`. Metrics to show (absent ones are skipped).
#' @param labels `character`. Sample labels, in the order of `wide$sample_id`.
#' @param fill `character` or `NULL`. Sample grouping used as bar colour.
#' @return A `ggplot` object, one facet per metric with its own y scale.
#' @examples
#' w <- data.frame(sample_id = c("s1", "s2"), library_size = c(1e6, 2e6), pct_counts_mito = c(2, 3))
#' plot_qc_metrics(w, c("library_size", "pct_counts_mito"), c("S1", "S2"))
#' @export
plot_qc_metrics <- function(wide, metrics, labels, fill = NULL) {
  metrics <- intersect(metrics, names(wide))
  if (length(metrics) == 0) abort("plot_qc_metrics(): none of the requested metrics is available")
  d <- do.call(rbind, lapply(metrics, function(mt) {
    data.frame(
      label = factor(labels, levels = unique(labels)), metric = mt, value = wide[[mt]],
      group = if (is.null(fill)) "all" else as.character(fill), stringsAsFactors = FALSE
    )
  }))
  d$metric <- factor(d$metric, levels = metrics)
  d$text <- ifelse(abs(d$value) >= 100, formatC(d$value, format = "fg", big.mark = ","), sprintf("%.1f", d$value))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$label, y = .data$value, fill = .data$group)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = .data$text), vjust = -0.3, size = 2.6, colour = "#2b2b2b") +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 3) +
    ggplot2::scale_fill_manual(values = palette_bulkrnaseq(unique(d$group)), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(title = "QC metrics per sample", x = NULL, y = NULL) +
    theme_bulkrnaseq() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = if (is.null(fill)) "none" else "right"
    )
}

# Categorical colours for sample groups (validated categorical palette slots, in order).
palette_bulkrnaseq <- function(levels) {
  pal <- c("#2a78d6", "#eb6834", "#1f9e89", "#9b59b6", "#c9a227", "#5a5a5a")
  stats::setNames(rep_len(pal, length(levels)), levels)
}

#' QC module
#'
#' Expression QC from the counts (when available), joined with the QC metrics imported at ingest (STAR `N_*`
#' rows, MultiQC reports, extra tables) into one row per sample.
#'
#' @param inputs `list` with `se` (optional, with counts) and `qc` (optional QC long table from ingest).
#' @param params `list`. `modules.qc` from the config.
#' @param context `list` with `cfg`, `gene_lists`, `organism`.
#' @return `list(tables, figures, meta)`; tables `qc_summary`, `expression_qc`, `top_genes`, figure
#'   `qc_metrics`.
#' @examples
#' m <- matrix(c(90L, 10L, 50L, 40L), 2, dimnames = list(c("a", "b"), c("s1", "s2")))
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   rowData = S4Vectors::DataFrame(
#'     feature_id = c("a", "b"), gene_symbol = c("A", "B"), biotype = "lncRNA", chromosome = "1"
#'   ),
#'   colData = S4Vectors::DataFrame(sample_id = c("s1", "s2"), row.names = c("s1", "s2"))
#' )
#' res <- run_qc(list(se = se), list(top_n_list = 1), list(cfg = list(), gene_lists = list(), organism = "human"))
#' res$tables$qc_summary
#' @export
run_qc <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "qc"))
  cfg <- context$cfg
  se <- inputs$se
  has_counts <- !is.null(se) && "counts" %in% SummarizedExperiment::assayNames(se)
  tables <- list()
  meta <- list()
  qc <- inputs$qc
  if (has_counts) {
    eq <- expression_qc_from_params(se, params, context)
    tables$expression_qc <- eq
    tables$top_genes <- top_genes(se, params$top_n_list)
    qc <- rbind(qc, eq)
    meta$expression_qc <- list(
      description = "Expression-level QC per sample, on all genes of the count matrix (no filter).",
      units = list(
        library_size = "reads assigned to genes", genes_count_ge1 = "genes", pct_counts_mito = "percent of library",
        pct_tpm_mito = "percent of TPM"
      )
    )
    meta$top_genes <- list(
      description = "Most expressed genes per sample.", units = list(pct_of_library = "percent of assigned reads")
    )
  }
  if (is.null(qc) || nrow(qc) == 0) abort("qc module: no QC metrics available (no counts and no imported QC)")
  samples <- if (has_counts) colnames(se) else unique(qc$sample_id)
  wide <- qc_wide(qc, union(samples, unique(qc$sample_id)))
  tables$qc_summary <- wide
  meta$qc_summary <- list(
    description = "QC metrics, one row per sample: imported QC (STAR, MultiQC, extra tables) and expression QC.",
    sources = unique(qc$source)
  )
  figures <- list()
  if (has_counts) {
    f <- sample_factors(se, cfg)
    fill <- if (length(f)) SummarizedExperiment::colData(se)[[f[[1]]]]
    w <- wide[match(colnames(se), wide$sample_id), , drop = FALSE]
    figures$qc_metrics <- plot_qc_metrics(w, unlist(params$figure_metrics),
      labels = sample_labels(se, cfg), fill = fill
    )
  }
  list(tables = tables, figures = figures, meta = meta, figure_sizes = list(qc_metrics = c(10, 8)))
}

# Expression QC with the parameters of modules.qc (shared by the qc and pca modules).
expression_qc_from_params <- function(se, params, context) {
  params <- with_defaults(params, c("modules", "qc"))
  expression_qc(se,
    ribosomal_protein_regex = params$ribosomal_protein_regex,
    hemoglobin_genes = resolve_gene_list(params$hemoglobin_genes, context$gene_lists),
    rrna_biotypes = unlist(params$rrna_biotypes), top_n_share = params$top_n_share,
    organism = context$organism %||% "human"
  )
}
