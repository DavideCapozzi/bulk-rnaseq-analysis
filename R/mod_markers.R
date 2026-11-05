#' Markers module
#'
#' Per-sample expression of chosen marker genes, grouped in panels, with panel scores and log2 ratios between
#' chosen genes and samples. All values are descriptive.
#' * `marker_expression`: one row per panel, gene and sample: counts, CPM on the library size (all genes of the
#'   count matrix, no normalisation factors), and, when the inputs provide them, TPM (`abundance`) and the
#'   provider's normalised counts (`normalized`). A panel gene matches every feature with that symbol whose
#'   biotype is in `biotypes` (`null`: any); genes without a match are listed in the table's metadata.
#' * `marker_panel_scores`: per panel and sample, the sum of the CPM of the panel's genes.
#' * `marker_ratios`: per configured ratio, `log2((CPM numerator + pseudocount) / (CPM denominator + pseudocount))`
#'   of a gene in a sample over a gene in a sample (e.g. a sort marker in the positive over the negative fraction).
#'
#' @param inputs `list` with `se` (counts; optional `tpm` and `normalized` assays).
#' @param params `list`. `modules.markers` from the config.
#' @param context `list` with `cfg` and `gene_lists`.
#' @return `list(tables, figures, meta, figure_sizes)`; figures `markers_heatmap` (log2(CPM + 1) per gene and
#'   sample, one block per panel) and `marker_panels` (panel scores per sample).
#' @examples
#' m <- matrix(c(100L, 5L, 900L, 10L, 50L, 940L), 3, dimnames = list(c("E1", "E2", "E3"), c("s1", "s2")))
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   rowData = S4Vectors::DataFrame(
#'     feature_id = rownames(m), gene_symbol = c("FOXP3", "IL7R", "ACTB"), biotype = "protein_coding"
#'   ),
#'   colData = S4Vectors::DataFrame(sample_id = colnames(m), row.names = colnames(m))
#' )
#' res <- run_markers(list(se = se), list(panels = list(list(id = "treg", genes = c("FOXP3", "IL7R")))), list())
#' res$tables$marker_panel_scores
#' @export
run_markers <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "markers"))
  se <- inputs$se
  if (is.null(se) || !"counts" %in% SummarizedExperiment::assayNames(se)) abort("markers: needs counts")
  if (length(params$panels) == 0) abort("markers: no panels configured (modules.markers.panels)")
  ids <- vapply(params$panels, `[[`, "", "id")
  if (anyDuplicated(ids)) abort("markers: duplicated panel id(s): ", toString(unique(ids[duplicated(ids)])))
  cfg <- context$cfg %||% list()
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  counts <- SummarizedExperiment::assay(se, "counts")
  cpm <- edgeR::cpm(counts)
  biotypes <- unlist(params$biotypes)
  if (!is.null(biotypes) && is.null(rd$biotype)) abort("markers: biotypes needs gene biotypes (annotation.gtf)")
  eligible <- if (is.null(biotypes)) rep(TRUE, nrow(rd)) else rd$biotype %in% biotypes
  panel_genes <- markers_panel_genes(params$panels, rd, eligible, context$gene_lists %||% list())
  found <- panel_genes[!is.na(panel_genes$row), , drop = FALSE]
  if (nrow(found) == 0) abort("markers: no panel gene matches a feature of the count matrix")
  samples <- colnames(se)
  rows <- rep(found$row, each = length(samples))
  cols <- rep(seq_along(samples), times = nrow(found))
  expr <- data.frame(
    panel_id = rep(found$panel_id, each = length(samples)),
    gene_symbol = rep(found$gene_symbol, each = length(samples)),
    feature_id = rep(rd$feature_id[found$row], each = length(samples)),
    sample_id = rep(samples, times = nrow(found)),
    count = counts[cbind(rows, cols)], cpm = cpm[cbind(rows, cols)], stringsAsFactors = FALSE
  )
  assays <- SummarizedExperiment::assayNames(se)
  if ("tpm" %in% assays) expr$tpm <- SummarizedExperiment::assay(se, "tpm")[cbind(rows, cols)]
  if ("normalized" %in% assays) expr$normalized <- SummarizedExperiment::assay(se, "normalized")[cbind(rows, cols)]
  pk <- paste(expr$panel_id, expr$sample_id, sep = "\r")
  key <- factor(pk, levels = unique(pk))
  first <- !duplicated(pk)
  scores <- data.frame(
    panel_id = expr$panel_id[first], sample_id = expr$sample_id[first],
    n_genes = as.integer(tapply(expr$cpm, key, length)), cpm_sum = as.numeric(tapply(expr$cpm, key, sum)),
    stringsAsFactors = FALSE
  )
  ratios <- markers_ratios(params$ratios, params$pseudocount, rd, eligible, cpm)
  missing <- panel_genes[is.na(panel_genes$row), c("panel_id", "gene_symbol"), drop = FALSE]
  labels <- sample_labels(se, cfg)
  panel_labels <- stats::setNames(vapply(params$panels, function(p) p$label %||% p$id, ""), ids)
  figures <- list(
    markers_heatmap = plot_markers_heatmap(expr, panel_labels, samples, labels),
    marker_panels = plot_marker_panels(scores, panel_labels, samples, labels, markers_fill(se, cfg))
  )
  n_rows <- length(unique(paste(found$panel_id, found$row)))
  list(
    tables = list(marker_expression = expr, marker_panel_scores = scores, marker_ratios = ratios),
    figures = figures,
    figure_sizes = list(
      markers_heatmap = c(4.5 + 0.6 * length(samples), max(4, 1.5 + 0.17 * n_rows + 0.15 * length(ids))),
      marker_panels = c(11, 2 + 2.2 * ceiling(length(ids) / 4))
    ),
    meta = list(
      marker_expression = list(
        description = paste0(
          "Expression of the panel genes per sample. cpm: counts per million of the library size (all genes of ",
          "the count matrix, no normalisation factors). Genes are matched by symbol",
          if (is.null(biotypes)) "" else paste0(" among features of biotype ", toString(biotypes)), "."
        ),
        panels = lapply(params$panels, function(p) {
          list(id = p$id, label = p$label %||% p$id, genes = resolve_gene_list(p$genes, context$gene_lists %||% list()))
        }),
        biotypes = biotypes %||% NA_character_,
        genes_not_found = if (nrow(missing)) paste(missing$panel_id, missing$gene_symbol, sep = ":") else character(),
        units = list(count = "reads", cpm = "counts per million", tpm = "TPM", normalized = "provider units"),
        caveats = "Panels are starting markers, not validated signatures; CPM are library-size only (no TMM)."
      ),
      marker_panel_scores = list(
        description = "Sum of the CPM of the panel genes found, per sample (descriptive).",
        units = list(cpm_sum = "counts per million")
      ),
      marker_ratios = list(
        description = paste0(
          "log2((CPM numerator + pseudocount) / (CPM denominator + pseudocount)) of a gene in a sample over a gene ",
          "in a sample; pseudocount ", params$pseudocount, "."
        ),
        pseudocount = params$pseudocount, units = list(log2_ratio = "log2 ratio")
      )
    )
  )
}

# One row per (panel, gene, matching feature), in panel and gene order; row = NA when the gene has no feature.
markers_panel_genes <- function(panels, rd, eligible, gene_lists) {
  out <- lapply(panels, function(p) {
    genes <- unique(resolve_gene_list(p$genes, gene_lists))
    if (length(genes) == 0) abort("markers: panel '", p$id, "' has no genes")
    do.call(rbind, lapply(genes, function(g) {
      hit <- which(rd$gene_symbol == g & eligible)
      data.frame(
        panel_id = p$id, gene_symbol = g, row = if (length(hit)) hit else NA_integer_, stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, out)
}

markers_ratios <- function(ratios, pseudocount, rd, eligible, cpm) {
  empty <- data.frame(
    ratio_id = character(), label = character(), numerator_gene = character(), numerator_sample = character(),
    denominator_gene = character(), denominator_sample = character(), numerator_cpm = numeric(),
    denominator_cpm = numeric(), log2_ratio = numeric(), stringsAsFactors = FALSE
  )
  ids <- vapply(ratios, `[[`, "", "id")
  if (anyDuplicated(ids)) abort("markers: duplicated ratio id(s): ", toString(unique(ids[duplicated(ids)])))
  value <- function(r, side) {
    g <- r[[side]]$gene
    s <- r[[side]]$sample
    hit <- which(rd$gene_symbol == g & eligible)
    if (length(hit) != 1) {
      abort(
        "markers: ratio '", r$id, "': gene '", g, "' matches ", length(hit), " features (needs exactly one; ",
        "restrict biotypes)"
      )
    }
    if (!s %in% colnames(cpm)) abort("markers: ratio '", r$id, "': unknown sample_id '", s, "'")
    cpm[hit, s]
  }
  out <- lapply(ratios, function(r) {
    num <- value(r, "numerator")
    den <- value(r, "denominator")
    data.frame(
      ratio_id = r$id, label = r$label %||% r$id, numerator_gene = r$numerator$gene,
      numerator_sample = r$numerator$sample, denominator_gene = r$denominator$gene,
      denominator_sample = r$denominator$sample, numerator_cpm = num, denominator_cpm = den,
      log2_ratio = log2((num + pseudocount) / (den + pseudocount)), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(list(empty), out))
}

# Bar colour of the panel figure: the first experimental factor of the sample sheet, if any.
markers_fill <- function(se, cfg) {
  f <- sample_factors(se, cfg)
  if (length(f)) as.character(SummarizedExperiment::colData(se)[[f[[1]]]])
}

#' Heatmap of marker expression
#'
#' @param expr `marker_expression` table of [run_markers()].
#' @param panels `character`. Panel labels in display order, named by panel id.
#' @param samples,labels `character`. Sample ids in display order and their labels.
#' @return A `ggplot` object: log2(CPM + 1) per gene and sample, one block per panel, values printed in the cells.
#' @examples
#' e <- data.frame(
#'   panel_id = "p", gene_symbol = rep(c("A", "B"), each = 2), feature_id = rep(c("E1", "E2"), each = 2),
#'   sample_id = c("s1", "s2"), cpm = c(10, 100, 1, 0)
#' )
#' plot_markers_heatmap(e, c(p = "Panel"), c("s1", "s2"), c("S1", "S2"))
#' @export
plot_markers_heatmap <- function(expr, panels, samples, labels) {
  d <- expr
  d$value <- log2(d$cpm + 1)
  d$label <- factor(labels[match(d$sample_id, samples)], levels = unique(labels))
  d$panel <- factor(panels[d$panel_id], levels = unique(panels))
  # One row per panel and feature; a symbol matching several features keeps its feature ID in the label.
  pg <- paste(d$panel_id, d$gene_symbol, sep = "\r")
  feats <- unique(data.frame(pg = pg, feature_id = d$feature_id))
  multi <- pg %in% feats$pg[duplicated(feats$pg)]
  d$row <- ifelse(multi, paste0(d$gene_symbol, " (", d$feature_id, ")"), d$gene_symbol)
  d$row <- factor(d$row, levels = rev(unique(d$row)))
  d$text_colour <- ifelse(d$value > 7, "white", "#1b1b1b")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$label, y = .data$row, fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", .data$value), colour = .data$text_colour), size = 2.3) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradientn(
      colours = c("#f7f7f5", "#9ec5f4", "#2a78d6", "#0d366b"), values = c(0, 1 / 3, 2 / 3, 1),
      limits = c(0, max(12, max(d$value, na.rm = TRUE))), name = "log2(CPM+1)"
    ) +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::labs(
      title = "Marker genes per sample", subtitle = "log2(CPM + 1), CPM of the library size", x = NULL, y = NULL
    ) +
    theme_bulkrnaseq(base_size = 9) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(), strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
}

#' Bar chart of marker panel scores
#'
#' @param scores `marker_panel_scores` table of [run_markers()].
#' @param panels `character`. Panel labels in display order, named by panel id.
#' @param samples,labels `character`. Sample ids in display order and their labels.
#' @param fill `character` or `NULL`. Sample grouping used as bar colour, in the order of `samples`.
#' @return A `ggplot` object, one facet per panel with its own y scale.
#' @examples
#' s <- data.frame(panel_id = "p", sample_id = c("s1", "s2"), n_genes = 2L, cpm_sum = c(10, 30))
#' plot_marker_panels(s, c(p = "Panel"), c("s1", "s2"), c("S1", "S2"))
#' @export
plot_marker_panels <- function(scores, panels, samples, labels, fill = NULL) {
  d <- scores
  i <- match(d$sample_id, samples)
  d$label <- factor(labels[i], levels = unique(labels))
  d$group <- if (is.null(fill)) "all" else fill[i]
  d$panel <- factor(panels[d$panel_id], levels = unique(panels))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$label, y = .data$cpm_sum, fill = .data$group)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 4) +
    ggplot2::scale_fill_manual(values = palette_bulkrnaseq(unique(d$group)), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = "Sum of the CPM per marker panel", x = NULL, y = "CPM (sum)") +
    theme_bulkrnaseq(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = if (is.null(fill)) "none" else "right"
    )
}
