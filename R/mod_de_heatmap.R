#' DE heatmap module
#'
#' One heatmap per contrast with a significance metric: the `n_top` most significant genes among those classified
#' UP or DOWN by the project thresholds (lowest padj, or highest posterior probability, first; ties by larger
#' |log2 fold change|, then feature_id; with `per_direction`, up to half from each direction), shown as per-gene
#' z-scores of log2(normalised count + 1) across all samples with counts. Normalised counts are the `normalized`
#' input when present (`values: auto` or `normalized`), else counts divided by DESeq2 median-of-ratios size factors
#' computed on the genes of the contrast table (`median_of_ratios`). Descriptive: the genes were selected on the
#' same samples.
#'
#' @param inputs `list` with `contrast_stats` and `se` (counts; optional `normalized` assay).
#' @param params `list`. `modules.de_heatmap` from the config.
#' @param context `list` with `cfg` (thresholds, sample factors and labels).
#' @return `list(tables, figures, meta, figure_sizes)`; tables `de_heatmap_genes` (selected genes, rank and row order)
#'   and `de_heatmap_values` (log2 value and z-score per gene and sample); one figure `de_heatmap_<contrast_id>` per
#'   contrast.
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- matrix(rnbinom(600, mu = 100, size = 10), 100, dimnames = list(paste0("G", 1:100), paste0("s", 1:6)))
#' m[1:10, 4:6] <- m[1:10, 4:6] * 8L
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   colData = S4Vectors::DataFrame(sample_id = colnames(m), row.names = colnames(m))
#' )
#' ct <- new_contrast_table("b_vs_a", rownames(m), c(rep(3, 10), rep(0, 90)),
#'   padj = c(rep(0.001, 10), rep(0.9, 90)), significance_metric = "padj", method = "m", source = "computed"
#' )
#' res <- run_de_heatmap(list(se = se, contrast_stats = ct), list(), list(cfg = list()))
#' head(res$tables$de_heatmap_genes)
#' }
#' @export
run_de_heatmap <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "de_heatmap"))
  cfg <- context$cfg %||% list()
  se <- inputs$se
  if (is.null(se) || !"counts" %in% SummarizedExperiment::assayNames(se)) abort("de_heatmap: needs counts")
  ct <- inputs$contrast_stats
  if (is.null(ct)) abort("de_heatmap: no contrast tables")
  ct <- ct[ct$significance_metric != "none", , drop = FALSE]
  if (nrow(ct) == 0) abort("de_heatmap: no contrast with a significance metric")
  wanted <- unlist(params$contrasts)
  if (length(wanted) > 0) {
    missing <- setdiff(wanted, ct$contrast_id)
    if (length(missing) > 0) abort("de_heatmap: contrast(s) without a significance metric: ", toString(missing))
    ct <- ct[ct$contrast_id %in% wanted, , drop = FALSE]
  }
  thresholds <- with_defaults(cfg$thresholds, "thresholds")
  keys <- unique(ct[, c("contrast_id", "method")])
  if (length(wanted) > 0) keys <- keys[order(match(keys$contrast_id, wanted)), , drop = FALSE]
  ambiguous <- keys$contrast_id[duplicated(keys$contrast_id)]
  samples <- de_heatmap_sample_order(se, cfg)
  labels <- sample_labels(se, cfg)[match(samples, colnames(se))]
  f <- sample_factors(se, cfg)
  group <- if (length(f)) as.character(SummarizedExperiment::colData(se)[samples, f[[1]]])
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  genes_out <- list()
  values_out <- list()
  figures <- list()
  sizes <- list()
  sources <- character()
  for (i in seq_len(nrow(keys))) {
    k <- keys$contrast_id[[i]]
    x <- ct[ct$contrast_id == k & ct$method == keys$method[[i]], , drop = FALSE]
    norm <- de_heatmap_matrix(se, x$feature_id, params$values)
    sources[[paste(k, keys$method[[i]])]] <- attr(norm, "source")
    sel <- de_heatmap_select(x, thresholds, params$n_top, rownames(norm), isTRUE(params$per_direction))
    if (nrow(sel) == 0) next
    lv <- log2(norm[sel$feature_id, samples, drop = FALSE] + 1)
    z <- t(scale(t(lv)))
    z[is.na(z)] <- 0
    ord <- if (isTRUE(params$cluster_rows) && nrow(z) > 1) {
      stats::hclust(stats::dist(z), method = "complete")$order
    } else {
      seq_len(nrow(z))
    }
    symbol <- rd$gene_symbol_unique[match(sel$feature_id, rd$feature_id)] %||% sel$gene_symbol
    symbol[is.na(symbol)] <- sel$feature_id[is.na(symbol)]
    sel$row_order <- match(seq_len(nrow(sel)), ord)
    genes_out[[length(genes_out) + 1]] <- sel
    values_out[[length(values_out) + 1]] <- data.frame(
      contrast_id = k, method = keys$method[[i]], feature_id = rep(sel$feature_id, times = length(samples)),
      gene_symbol = rep(symbol, times = length(samples)), sample_id = rep(samples, each = nrow(sel)),
      log2_value = as.vector(lv), z = as.vector(z), stringsAsFactors = FALSE
    )
    name <- paste0("de_heatmap_", k, if (k %in% ambiguous) paste0("_", keys$method[[i]]))
    figures[[name]] <- plot_de_heatmap(z, symbol, ord, labels, group, params$z_limit,
      title = paste0("Top ", nrow(sel), " DE genes: ", k), subtitle = paste0(
        "z-score of log2(", attr(norm, "source"), " + 1) per gene; genes ranked by ", x$significance_metric[[1]],
        " among UP/DOWN", if (isTRUE(params$per_direction)) " (half per direction)"
      )
    )
    sizes[[name]] <- c(3.5 + 0.5 * length(samples), max(4, 1.8 + 0.14 * nrow(sel)))
  }
  genes <- do.call(rbind, c(list(de_heatmap_empty_genes()), genes_out))
  values <- do.call(rbind, c(list(de_heatmap_empty_values()), values_out))
  rownames(genes) <- NULL
  rownames(values) <- NULL
  list(
    tables = list(de_heatmap_genes = genes, de_heatmap_values = values),
    figures = figures, figure_sizes = sizes,
    meta = list(
      de_heatmap_genes = list(
        description = paste0(
          "Per contrast, the ", params$n_top, " most significant genes classified UP or DOWN (thresholds below); ",
          "rank 1 = most significant; row_order = position in the heatmap (row clustering: ",
          if (isTRUE(params$cluster_rows)) "complete linkage on Euclidean distances of the z-scores" else "none", ")."
        ),
        thresholds = thresholds, values = as.list(sources),
        units = list(lfc = "log2_ratio", significance = "as named in significance_metric")
      ),
      de_heatmap_values = list(
        description = paste0(
          "log2(normalised count + 1) and its z-score across the samples shown, per selected gene (normalised ",
          "counts: ", paste(unique(sources), collapse = "; "), ")."
        ),
        units = list(log2_value = "log2(normalised count + 1)", z = "standard deviations across samples"),
        caveats = "Genes were selected on the same samples: the separation of the groups in the heatmap is expected."
      )
    )
  )
}

de_heatmap_empty_genes <- function() {
  data.frame(
    contrast_id = character(), method = character(), rank = integer(), feature_id = character(),
    gene_symbol = character(), lfc = numeric(), significance_metric = character(), significance = numeric(),
    class = character(), row_order = integer(), stringsAsFactors = FALSE
  )
}

de_heatmap_empty_values <- function() {
  data.frame(
    contrast_id = character(), method = character(), feature_id = character(), gene_symbol = character(),
    sample_id = character(), log2_value = numeric(), z = numeric(), stringsAsFactors = FALSE
  )
}

# Samples ordered by the levels of the first experimental factor, then sample sheet order.
de_heatmap_sample_order <- function(se, cfg) {
  f <- sample_factors(se, cfg)
  s <- colnames(se)
  if (length(f) == 0) {
    return(s)
  }
  v <- as.character(SummarizedExperiment::colData(se)[[f[[1]]]])
  lv <- unlist(cfg$samples$factors[[f[[1]]]])
  s[order(match(v, lv), seq_along(s))]
}

# Normalised counts (genes x samples) with attribute "source".
de_heatmap_matrix <- function(se, tested, values) {
  has_norm <- "normalized" %in% SummarizedExperiment::assayNames(se)
  if (values == "normalized" && !has_norm) abort("de_heatmap: values 'normalized' needs the normalized input")
  if (has_norm && values %in% c("auto", "normalized")) {
    m <- SummarizedExperiment::assay(se, "normalized")
    m <- m[stats::complete.cases(m), , drop = FALSE]
    attr(m, "source") <- "normalized input"
    return(m)
  }
  counts <- SummarizedExperiment::assay(se, "counts")
  ids <- intersect(tested, rownames(counts))
  if (length(ids) == 0) abort("de_heatmap: no gene of the contrast table is in the count matrix")
  sf <- DESeq2::estimateSizeFactorsForMatrix(counts[ids, , drop = FALSE])
  m <- sweep(counts, 2, sf, "/")
  attr(m, "source") <- "counts / median-of-ratios size factors"
  m
}

# Most significant UP/DOWN genes; with per_direction, up to ceiling(n_top / 2) of each class (UP first), kept in
# significance order.
de_heatmap_select <- function(x, thresholds, n_top, available, per_direction = FALSE) {
  sig <- significance_values(x)
  cls <- classify_de(x, thresholds)
  keep <- cls != "NS" & x$feature_id %in% available
  s <- if (sig$metric == "padj") sig$values else -sig$values
  o <- order(s[keep], -abs(x$lfc[keep]), x$feature_id[keep], method = "radix")
  ranked <- which(keep)[o]
  idx <- if (per_direction) {
    half <- ceiling(n_top / 2)
    up <- ranked[cls[ranked] == "UP"]
    down <- ranked[cls[ranked] == "DOWN"]
    sort_by_rank <- function(v) v[order(match(v, ranked))]
    sort_by_rank(c(up[seq_len(min(half, length(up)))], down[seq_len(min(half, length(down)))]))
  } else {
    ranked[seq_len(min(n_top, length(ranked)))]
  }
  data.frame(
    contrast_id = x$contrast_id[idx], method = x$method[idx], rank = seq_along(idx), feature_id = x$feature_id[idx],
    gene_symbol = x$gene_symbol[idx], lfc = x$lfc[idx], significance_metric = rep(sig$metric, length(idx)),
    significance = sig$values[idx], class = cls[idx], stringsAsFactors = FALSE
  )
}

#' Heatmap of the top differentially expressed genes
#'
#' @param z Numeric matrix, genes x samples (z-scores), rows in selection order.
#' @param symbols `character`. Row labels, in the order of `z`.
#' @param order `integer`. Display order of the rows (e.g. from clustering), top to bottom.
#' @param labels `character`. Column labels, in the order of `z`'s columns.
#' @param group `character` or `NULL`. Sample grouping (columns are split by it).
#' @param z_limit `numeric(1)`. Colour scale limit; larger |z| are drawn at the limit.
#' @param title,subtitle `character(1)`. Plot titles.
#' @return A `ggplot` object.
#' @examples
#' z <- matrix(c(-1, 1, 1, -1), 2, dimnames = list(c("E1", "E2"), c("s1", "s2")))
#' plot_de_heatmap(z, c("A", "B"), 1:2, c("S1", "S2"), c("a", "b"))
#' @export
plot_de_heatmap <- function(z, symbols, order, labels, group = NULL, z_limit = 2.5, title = NULL, subtitle = NULL) {
  d <- data.frame(
    row = factor(rep(symbols, times = ncol(z)), levels = rev(unique(symbols[order]))),
    label = factor(rep(labels, each = nrow(z)), levels = unique(labels)),
    group = if (is.null(group)) "all" else factor(rep(group, each = nrow(z)), levels = unique(group)),
    z = pmax(pmin(as.vector(z), z_limit), -z_limit), stringsAsFactors = FALSE
  )
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$label, y = .data$row, fill = .data$z)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low = "#2a78d6", mid = "#f7f7f5", high = "#eb6834", midpoint = 0, limits = c(-z_limit, z_limit),
      name = "z-score"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_bulkrnaseq(base_size = 9) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      axis.text.y = ggplot2::element_text(size = if (nrow(z) > 60) 5 else 7)
    )
  if (!is.null(group)) p <- p + ggplot2::facet_grid(. ~ group, scales = "free_x", space = "free_x")
  p
}
