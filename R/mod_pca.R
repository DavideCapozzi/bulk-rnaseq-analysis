#' PCA on the most variable genes
#'
#' Genes are ranked by variance across samples and the top `ntop` are kept (all genes when `ntop` is `Inf`;
#' they are reordered by variance in both cases), then `prcomp()` with centring and no scaling.
#'
#' @param x Numeric matrix, genes x samples (e.g. VST or log CPM).
#' @param ntop `numeric(1)`. Number of genes, `Inf` for all.
#' @return `list(x = scores, ve = percent variance per PC, rotation = loadings, n_genes)`.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(40), 10, dimnames = list(paste0("g", 1:10), paste0("s", 1:4)))
#' pca_top_variance(x, 5)$ve
#' @export
pca_top_variance <- function(x, ntop = Inf) {
  v <- apply(x, 1, stats::var)
  sel <- order(v, decreasing = TRUE)[seq_len(min(ntop, nrow(x)))]
  p <- stats::prcomp(t(x[sel, , drop = FALSE]), center = TRUE, scale. = FALSE)
  list(x = p$x, ve = 100 * p$sdev^2 / sum(p$sdev^2), rotation = p$rotation, n_genes = length(sel))
}

#' Which factor separates samples along a principal component
#'
#' For each factor, the range of the per-level mean scores (the absolute mean difference for two levels);
#' with exactly two two-level factors, also their interaction (samples where the two factors "agree" versus
#' the others). The label is the factor with the largest value, ties going to the first in `factors` order.
#' A descriptive label, not a test.
#'
#' @param scores `numeric`. PC scores, one per sample.
#' @param samples `data.frame` of sample metadata in the same order.
#' @param factors `character`. Sample sheet columns, in priority order.
#' @return `character(1)`, a factor name, `"interaction"`, or `NA` when `factors` is empty.
#' @examples
#' ss <- data.frame(f = c("a", "a", "b", "b"), r = c("x", "y", "x", "y"))
#' pc_separation(c(-2, -1.8, 2, 1.9), ss, c("f", "r"))
#' @export
pc_separation <- function(scores, samples, factors) {
  if (length(factors) == 0) {
    return(NA_character_)
  }
  sep <- vapply(factors, function(f) {
    lv <- as.character(samples[[f]])
    means <- vapply(unique(lv), function(l) mean(scores[lv == l]), numeric(1))
    max(means) - min(means)
  }, numeric(1))
  codes <- lapply(factors, function(f) {
    lv <- factor(as.character(samples[[f]]), levels = levels_of(samples[[f]]))
    if (nlevels(droplevels(lv)) == 2) ifelse(lv == levels(droplevels(lv))[[2]], 1, -1)
  })
  if (length(factors) == 2 && !any(vapply(codes, is.null, logical(1)))) {
    prod_code <- codes[[1]] * codes[[2]]
    sep <- c(sep, interaction = abs(mean(scores[prod_code < 0]) - mean(scores[prod_code > 0])))
  }
  names(sep)[which.max(sep)]
}

levels_of <- function(x) if (is.factor(x)) levels(x) else unique(as.character(x))

#' Variance-stabilising transformation of counts
#'
#' DESeq2 `varianceStabilizingTransformation(blind = TRUE)` with an intercept-only design, so the transform
#' does not see the experimental groups.
#'
#' @param counts Integer matrix, genes x samples.
#' @return Numeric matrix on the VST scale (approximately log2).
#' @examples
#' set.seed(1)
#' m <- matrix(rnbinom(4000, mu = 100, size = 5), 1000, dimnames = list(NULL, paste0("s", 1:4)))
#' dim(vst_matrix(m))
#' @export
vst_matrix <- function(counts) {
  ids <- colnames(counts) %||% paste0("s", seq_len(ncol(counts)))
  colnames(counts) <- ids
  # Thinned counts are whole doubles; DESeq2 would convert them itself, with a message.
  if (!is.integer(counts) && is_whole_number(counts)) storage.mode(counts) <- "integer"
  cd <- data.frame(sample = ids, row.names = ids)
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, cd, design = ~1)
  SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
}

# Sign that orients PC1 so that the samples at `orient$level` of `orient$factor` have the higher mean score.
orientation_sign <- function(scores, samples, orient) {
  if (is.null(orient) || is.null(orient$factor)) {
    return(1)
  }
  if (!orient$factor %in% names(samples)) abort("modules.pca.orient_pc1: no sample sheet column '", orient$factor, "'")
  at <- as.character(samples[[orient$factor]]) == orient$level
  if (!any(at) || all(at)) abort("modules.pca.orient_pc1: level '", orient$level, "' must select some samples, not all")
  if (mean(scores[at]) < mean(scores[!at])) -1 else 1
}

# Direction labels of positive and negative PC1 loadings after orientation.
loading_directions <- function(samples, orient) {
  if (is.null(orient) || is.null(orient$factor)) {
    return(c("positive", "negative"))
  }
  other <- setdiff(levels_of(samples[[orient$factor]]), orient$level)
  negative <- if (length(other) == 1) paste0("towards_", other) else paste0("away_from_", orient$level)
  c(paste0("towards_", orient$level), negative)
}

# All ordered pairs of a sample x sample matrix as a long table.
pair_table <- function(mat, transform, method) {
  data.frame(
    transform = transform, method = method, sample_a = rep(rownames(mat), times = ncol(mat)),
    sample_b = rep(colnames(mat), each = nrow(mat)), value = as.vector(mat), stringsAsFactors = FALSE
  )
}

#' PC1 drivers: technical metrics and gene features
#'
#' PCA on all filtered genes (VST), PC1 oriented by `orient`; correlates PC1 scores with per-sample QC
#' metrics and PC1 loadings with gene span, lncRNA status and mean expression. With few samples the
#' correlations are descriptive.
#'
#' @param x VST matrix of the filtered genes.
#' @param features `data.frame` (rowData) of the same genes with `feature_id`, `gene_symbol`, `biotype`,
#'   `start`, `end`.
#' @param samples `data.frame` of sample metadata in column order.
#' @param qc `data.frame` from [qc_wide()] or `NULL`.
#' @param metrics `character`. QC metrics to correlate (absent ones are skipped).
#' @param orient `list(factor, level)` or `NULL`.
#' @param top_n `integer(1)`. Genes listed per direction in `pc1_top_loadings`.
#' @return `list(pc1_vs_qc, pc1_loading_vs_gene_features, pc1_top_loadings, genes)`; `genes` holds the
#'   per-gene loadings for the figure.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(400, 8), 100, dimnames = list(paste0("g", 1:100), paste0("s", 1:4)))
#' f <- data.frame(
#'   feature_id = rownames(x), gene_symbol = rownames(x), biotype = "protein_coding",
#'   start = 1, end = seq(1000, 1e5, length.out = 100)
#' )
#' d <- pc_drivers(x, f, data.frame(g = c("a", "a", "b", "b")), NULL, character(), top_n = 3)
#' d$pc1_top_loadings
#' @export
pc_drivers <- function(x, features, samples, qc = NULL, metrics = character(), orient = NULL, top_n = 40) {
  p <- stats::prcomp(t(x), center = TRUE)
  sc <- p$x[, seq_len(min(3, ncol(p$x))), drop = FALSE]
  s <- orientation_sign(sc[, 1], samples, orient)
  sc[, 1] <- s * sc[, 1]
  loading <- s * p$rotation[, 1]
  pcqc <- NULL
  for (v in intersect(metrics, names(qc))) {
    q <- qc[[v]]
    pc_cor <- function(k, method = "pearson") {
      if (k <= ncol(sc)) suppressWarnings(stats::cor(sc[, k], q, method = method)) else NA
    }
    pcqc <- rbind(pcqc, data.frame(
      metric = v, values = paste(signif(q, 4), collapse = ";"), spearman_pc1 = pc_cor(1, "spearman"),
      pearson_pc1 = pc_cor(1), pearson_pc2 = pc_cor(2), pearson_pc3 = pc_cor(3), stringsAsFactors = FALSE
    ))
  }
  if (!is.null(pcqc)) pcqc$pc1_scores <- paste(round(sc[, 1], 1), collapse = ";")
  span_kb <- (features$end - features$start + 1) / 1e3
  is_lnc <- features$biotype %in% "lncRNA"
  mean_vst <- rowMeans(x)
  ok <- !is.na(span_kb)
  q <- stats::quantile(loading, c(0.05, 0.95))
  top_pos <- loading >= q[[2]]
  top_neg <- loading <= q[[1]]
  dirs <- loading_directions(samples, orient)
  feat <- data.frame(
    feature = c(
      "spearman_loading_vs_log10_span", "spearman_loading_vs_lncrna", "spearman_loading_vs_mean_vst",
      paste0("median_span_kb_top5pct_", dirs), "median_span_kb_all", paste0("pct_lncrna_top5pct_", dirs),
      "pct_lncrna_all"
    ),
    value = c(
      stats::cor(loading[ok], log10(span_kb[ok]), method = "spearman"),
      stats::cor(loading, as.numeric(is_lnc), method = "spearman"),
      stats::cor(loading, mean_vst, method = "spearman"),
      stats::median(span_kb[top_pos & ok]), stats::median(span_kb[top_neg & ok]), stats::median(span_kb[ok]),
      100 * mean(is_lnc[top_pos]), 100 * mean(is_lnc[top_neg]), 100 * mean(is_lnc)
    ),
    stringsAsFactors = FALSE
  )
  side <- function(o, direction) {
    o <- o[seq_len(min(top_n, length(o)))]
    data.frame(
      direction = direction, rank = seq_along(o), feature_id = features$feature_id[o],
      gene_symbol = features$gene_symbol[o], biotype = features$biotype[o], span_kb = span_kb[o],
      loading_pc1 = loading[o], stringsAsFactors = FALSE
    )
  }
  top <- rbind(side(order(-loading), dirs[[1]]), side(order(loading), dirs[[2]]))
  list(
    pc1_vs_qc = pcqc, pc1_loading_vs_gene_features = feat, pc1_top_loadings = top,
    genes = data.frame(span_kb = span_kb, loading_pc1 = loading)
  )
}

#' PCA module
#'
#' Runs a grid of PCAs (sex-gene variants x symbol exclusions x original/thinned counts x gene subsets x
#' transforms x numbers of top genes) and summarises the variance of PC1-3 and the factor each one separates;
#' writes scores for one gene number, PCA on log2(TPM + 1) when abundances are available, sample correlations
#' and distances on the main gene filter, and the PC1 drivers.
#'
#' @param inputs `list` with `se` (counts; optionally tpm), optional `qc` and `gene_flags`.
#' @param params `list`. `modules.pca` from the config.
#' @param context `list` with `cfg`, `gene_lists`, `organism`, `seed`.
#' @return `list(tables, figures, meta)`.
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- matrix(rnbinom(8000, mu = 200, size = 5), 2000, dimnames = list(paste0("g", 1:2000), paste0("s", 1:4)))
#' rd <- S4Vectors::DataFrame(
#'   feature_id = rownames(m), gene_symbol = rownames(m), biotype = "protein_coding",
#'   chromosome = "1", start = 1L, end = 5000L
#' )
#' cd <- S4Vectors::DataFrame(sample_id = colnames(m), g = c("a", "a", "b", "b"), row.names = colnames(m))
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m), rowData = rd, colData = cd)
#' cfg <- list(samples = list(factors = list(g = c("a", "b"))), gene_filter = list(sex_genes = "keep"))
#' res <- run_pca(
#'   list(se = se), list(ntop = list(500), sex_gene_variants = "keep", separation_factors = "g"),
#'   list(cfg = cfg, gene_lists = list(), organism = "human", seed = 1)
#' )
#' res$tables$pca_summary
#' }
#' @export
run_pca <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "pca"))
  cfg <- context$cfg
  se <- inputs$se
  counts <- SummarizedExperiment::assay(se, "counts")
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  samples <- as.data.frame(SummarizedExperiment::colData(se))
  labels <- sample_labels(se, cfg)
  # An empty list (the schema default) means every factor of samples.factors.
  factors <- unlist(params$separation_factors)
  if (length(factors) == 0) factors <- sample_factors(se, cfg)
  main_sex <- with_defaults(cfg$gene_filter, "gene_filter")$sex_genes
  keep_for <- function(sex) config_filter_genes(se, cfg, context, sex)
  thinned <- if (isTRUE(params$thinning)) {
    withr::with_seed(context$seed %||% 1, edgeR::thinCounts(counts, target.size = min(colSums(counts))))
  }
  grid <- pca_grid(counts, thinned, rd, samples, labels, factors, params, keep_for, main_sex,
    gene_flags = inputs$gene_flags, organism = context$organism %||% "human"
  )
  main <- grid$main
  if (is.null(main)) {
    # The main gene filter is not among the grid variants: correlations and drivers still use it.
    keep <- keep_for(main_sex)
    xk <- counts[keep, , drop = FALSE]
    main <- list(vst = vst_matrix(xk), logcpm_tmm = log_cpm(xk), keep = keep)
  }
  tables <- list(pca_summary = grid$summary, pca_scores = grid$scores)
  figures <- list()
  if (!is.null(tables$pca_scores)) figures$pca <- plot_pca_scores(tables$pca_scores, factors)
  if ("tpm" %in% SummarizedExperiment::assayNames(se)) {
    tpm <- SummarizedExperiment::assay(se, "tpm")[keep_for(params$abundance_sex_genes), , drop = FALSE]
    p <- pca_top_variance(log2(tpm + 1), Inf)
    tables$pca_scores_log2tpm <- pca_scores_table(p, samples, labels, factors, NULL, "log2_tpm_plus1", n_pc = 2)
  }
  extra <- pca_main_outputs(main, se, rd, samples, labels, inputs$qc, params, context)
  list(
    tables = c(tables, extra$tables), figures = c(figures, extra$figures),
    meta = list(
      pca_summary = list(
        description = paste(
          "Percent variance of PC1-3 and the factor with the largest mean difference on each (descriptive label)",
          "per sex-gene variant, symbol exclusion, data (original or thinned to the smallest library), gene subset,",
          "transform and number of most variable genes."
        ),
        units = list(pc1_var = "percent of total variance"), seed = if (isTRUE(params$thinning)) context$seed %||% 1,
        separation_factors = factors
      ),
      pca_scores = list(description = paste0("PC scores, ", params$scores_ntop, " most variable genes.")),
      sample_correlations = list(description = paste0("Sample correlations, main gene filter (", main_sex, ").")),
      sample_distances = list(description = paste0("Euclidean distances on VST, main gene filter (", main_sex, ").")),
      pc1_vs_qc = list(
        description = "PC1-3 scores (VST, all filtered genes) vs QC metrics; descriptive.", n_samples = ncol(se),
        min_two_sided_p_spearman = 2 / factorial(ncol(se))
      ),
      pc1_top_loadings = list(units = list(span_kb = "kb (gene start to end)"))
    ),
    figure_sizes = list(pca = c(11, 9), sample_correlation = c(6, 5), pc1_loading_vs_gene_span = c(7, 4.5))
  )
}

# The PCA grid: one pca_summary row per sex-gene variant x symbol exclusion x data (original, thinned) x gene
# subset x transform x number of genes; PC scores for scores_ntop on the unexcluded original data and first
# subset; and the matrices of the main gene filter (NULL when its sex-gene variant is not in the grid).
pca_grid <- function(counts, thinned, rd, samples, labels, factors, params, keep_for, main_sex, gene_flags,
                     organism) {
  ntops <- vapply(unlist(params$ntop), as.character, "")
  scores_ntop <- as.character(params$scores_ntop)
  transforms <- unlist(params$transforms)
  subsets <- unlist(params$gene_subsets)
  data_variants <- c("original", if (!is.null(thinned)) "thinned")
  summary_rows <- list()
  scores_rows <- list()
  main <- NULL
  for (sex in unlist(params$sex_gene_variants)) {
    keep_sex <- keep_for(sex)
    for (excl in c("none", unlist(params$symbol_exclusion_variants))) {
      keep <- keep_sex
      if (excl != "none") keep <- keep & !grepl(excl, rd$gene_symbol)
      genes <- merge_gene_flags(
        data.frame(
          feature_id = rd$feature_id[keep], chromosome = rd$chromosome[keep],
          expression = rowMeans(log_cpm(counts[keep, , drop = FALSE]))
        ),
        gene_flags
      )
      for (dat in data_variants) {
        x0 <- if (dat == "original") counts[keep, , drop = FALSE] else thinned[keep, , drop = FALSE]
        for (sub in subsets) {
          xs <- x0[gene_subset(genes, sub, organism), , drop = FALSE]
          mats <- list()
          if ("vst" %in% transforms) mats$vst <- vst_matrix(xs)
          if ("logcpm_tmm" %in% transforms) mats$logcpm_tmm <- log_cpm(xs)
          first <- excl == "none" && dat == "original" && sub == subsets[[1]]
          if (first && sex == main_sex) main <- c(mats, list(keep = keep))
          for (tr in names(mats)) {
            for (nt in ntops) {
              p <- pca_top_variance(mats[[tr]], if (identical(nt, "all")) Inf else as.numeric(nt))
              variant <- data.frame(
                sex_genes = sex, symbol_exclusion = excl, data = dat, gene_subset = sub, transform = tr,
                ntop = nt, n_genes = nrow(xs), stringsAsFactors = FALSE
              )
              summary_rows[[length(summary_rows) + 1]] <- cbind(
                variant, pca_summary_stats(p, samples, labels, factors, params$orient_pc1)
              )
              if (first && nt == scores_ntop) {
                scores_rows[[length(scores_rows) + 1]] <- pca_scores_table(p, samples, labels, factors, sex, tr)
              }
            }
          }
        }
      }
    }
  }
  list(summary = do.call(rbind, summary_rows), scores = do.call(rbind, scores_rows), main = main)
}

# Variance of PC1-3, the factor each one separates and the sample order along the oriented PC1.
pca_summary_stats <- function(p, samples, labels, factors, orient) {
  ve <- c(p$ve, NA, NA)[1:3]
  sep <- vapply(1:3, function(k) {
    if (k <= ncol(p$x)) pc_separation(p$x[, k], samples, factors) else NA_character_
  }, "")
  o <- order(p$x[, 1] * orientation_sign(p$x[, 1], samples, orient))
  data.frame(
    n_genes_used = p$n_genes, pc1_var = ve[[1]], pc2_var = ve[[2]], pc3_var = ve[[3]],
    pc1_separates = sep[[1]], pc2_separates = sep[[2]], pc3_separates = sep[[3]],
    pc1_order = paste(labels[o], collapse = " < "), stringsAsFactors = FALSE
  )
}

# Sample correlations and distances, and the PC1 drivers, on the main gene filter.
pca_main_outputs <- function(main, se, rd, samples, labels, qc, params, context) {
  tables <- list()
  figures <- list()
  pairs <- list()
  if (!is.null(main$logcpm_tmm)) {
    cs <- stats::cor(main$logcpm_tmm, method = "spearman")
    pairs$spearman <- pair_table(cs, "logcpm_tmm", "spearman")
    figures$sample_correlation <- plot_sample_correlation(cs, labels)
  }
  if (!is.null(main$vst)) {
    pairs$pearson <- pair_table(stats::cor(main$vst, method = "pearson"), "vst", "pearson")
    tables$sample_distances <- pair_table(as.matrix(stats::dist(t(main$vst))), "vst", "euclidean")
  }
  if (length(pairs)) tables$sample_correlations <- do.call(rbind, unname(pairs))
  if (!isFALSE(params$pc_drivers) && !is.null(main$vst)) {
    qc_all <- rbind(qc, expression_qc_from_params(se, context$cfg$modules$qc, context))
    d <- pc_drivers(main$vst, rd[main$keep, , drop = FALSE], samples, qc_wide(qc_all, colnames(se)),
      metrics = unlist(params$qc_metrics_for_drivers), orient = params$orient_pc1, top_n = params$top_loadings
    )
    if (!is.null(d$pc1_vs_qc)) tables$pc1_vs_qc <- d$pc1_vs_qc
    tables$pc1_loading_vs_gene_features <- d$pc1_loading_vs_gene_features
    tables$pc1_top_loadings <- d$pc1_top_loadings
    figures$pc1_loading_vs_gene_span <- plot_loading_vs_span(d$genes, loading_directions(samples, params$orient_pc1))
  }
  list(tables = tables, figures = figures)
}

pca_scores_table <- function(p, samples, labels, factors, sex, transform, n_pc = 3) {
  out <- data.frame(sample_id = rownames(samples), label = labels, stringsAsFactors = FALSE)
  for (f in factors) out[[f]] <- as.character(samples[[f]])
  if (!is.null(sex)) out$sex_genes <- sex
  out$transform <- transform
  for (k in seq_len(n_pc)) out[[paste0("pc", k)]] <- if (k <= ncol(p$x)) unname(p$x[, k]) else NA_real_
  for (k in seq_len(n_pc)) out[[paste0("ve", k)]] <- c(p$ve, rep(NA, n_pc))[[k]]
  out
}

#' Plot PC1 against PC2 for each variant
#'
#' @param scores `data.frame` from the `pca_scores` table.
#' @param factors `character`. First factor is the colour, second the shape.
#' @return A `ggplot` object, one panel per sex-gene variant and transform.
#' @examples
#' s <- data.frame(
#'   label = c("A", "B", "C"), g = c("x", "y", "y"), sex_genes = "keep", transform = "vst",
#'   pc1 = c(-1, 0.5, 0.5), pc2 = c(0, 1, -1), ve1 = 60, ve2 = 40
#' )
#' plot_pca_scores(s, "g")
#' @export
plot_pca_scores <- function(scores, factors = character()) {
  d <- scores
  d$panel <- sprintf("%s, %s\nPC1 %.1f%%, PC2 %.1f%%", d$transform, d$sex_genes, d$ve1, d$ve2)
  colour <- if (length(factors) >= 1) factors[[1]]
  shape <- if (length(factors) >= 2) factors[[2]]
  aes_pt <- ggplot2::aes(
    colour = if (!is.null(colour)) .data[[colour]], shape = if (!is.null(shape)) .data[[shape]]
  )
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$pc1, y = .data$pc2)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#d6d5d0", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#d6d5d0", linewidth = 0.3) +
    ggplot2::geom_point(aes_pt, size = 3.5) +
    ggrepel::geom_text_repel(ggplot2::aes(label = .data$label),
      colour = "#2b2b2b", size = 3, seed = 1,
      box.padding = 0.5, min.segment.length = 0
    ) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    ggplot2::labs(x = "PC1", y = "PC2", colour = colour, shape = shape) +
    theme_bulkrnaseq()
  if (!is.null(colour)) p <- p + ggplot2::scale_colour_manual(values = palette_bulkrnaseq(unique(d[[colour]])))
  p
}

#' Heatmap of sample correlations
#'
#' @param cmat Correlation matrix (samples x samples).
#' @param labels `character`. Sample labels in matrix order.
#' @return A `ggplot` object.
#' @examples
#' plot_sample_correlation(matrix(c(1, 0.9, 0.9, 1), 2, dimnames = list(c("a", "b"), c("a", "b"))), c("A", "B"))
#' @export
plot_sample_correlation <- function(cmat, labels) {
  d <- data.frame(
    a = factor(rep(labels, times = ncol(cmat)), levels = labels),
    b = factor(rep(labels, each = nrow(cmat)), levels = labels), rho = as.vector(cmat)
  )
  ggplot2::ggplot(d, ggplot2::aes(x = .data$a, y = .data$b, fill = .data$rho)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$rho)), size = 3.4, colour = "#1b1b1b") +
    ggplot2::scale_fill_gradient(low = "#cde2fb", high = "#256abf", name = "Spearman") +
    ggplot2::labs(title = "Sample correlation (log CPM, main gene filter)", x = NULL, y = NULL) +
    ggplot2::coord_fixed() +
    theme_bulkrnaseq()
}

#' PC1 loadings by decile of gene span
#'
#' @param genes `data.frame` with `span_kb` and `loading_pc1`.
#' @param directions `character(2)`. Labels of positive and negative loadings.
#' @return A `ggplot` object.
#' @examples
#' set.seed(1)
#' plot_loading_vs_span(data.frame(span_kb = rexp(200), loading_pc1 = rnorm(200)), c("positive", "negative"))
#' @export
plot_loading_vs_span <- function(genes, directions) {
  g <- genes[!is.na(genes$span_kb), , drop = FALSE]
  br <- unique(stats::quantile(log10(g$span_kb), 0:10 / 10))
  g$span_bin <- cut(log10(g$span_kb), breaks = br, include.lowest = TRUE, labels = paste0("D", seq_len(length(br) - 1)))
  ggplot2::ggplot(g, ggplot2::aes(x = .data$span_bin, y = .data$loading_pc1)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#b5b4af", linewidth = 0.4) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "#cde2fb", colour = "#256abf", width = 0.6, linewidth = 0.4) +
    ggplot2::coord_cartesian(ylim = stats::quantile(g$loading_pc1, c(0.01, 0.99))) +
    ggplot2::labs(
      title = "PC1 loading by decile of gene span",
      subtitle = paste0("Positive PC1 = ", directions[[1]]),
      x = "Decile of genomic span (D1 = shortest genes)", y = "PC1 loading (VST, all filtered genes)"
    ) +
    theme_bulkrnaseq()
}
