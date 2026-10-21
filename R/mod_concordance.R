#' Concordance between two vectors of fold changes
#'
#' Counts genes whose fold change exceeds a threshold in both contrasts and how many of them agree in sign.
#' Among genes above the threshold in both, 50% concordant signs is what independent random signs would give.
#'
#' @param x,y `numeric`. Fold changes (log2) of the same genes in two contrasts.
#' @param threshold `numeric(1)`. |log2 fold change| threshold (strict inequality).
#' @return Named `numeric`: `pearson`, `spearman`, `n_both`, `n_concordant`, `n_discordant`, `n_conc_up`,
#'   `n_conc_down`, `n_x_gt`, `n_y_gt`, `pct_concordant` (percentage; denominator at least 1).
#' @examples
#' concordance_metrics(c(2, -2, 0.1, 3), c(1.5, 2, 0, 2), threshold = 1)
#' @export
concordance_metrics <- function(x, y, threshold = 1) {
  both <- abs(x) > threshold & abs(y) > threshold
  same <- both & sign(x) == sign(y)
  c(
    pearson = stats::cor(x, y), spearman = stats::cor(x, y, method = "spearman"),
    n_both = sum(both), n_concordant = sum(same), n_discordant = sum(both & sign(x) != sign(y)),
    n_conc_up = sum(x > threshold & y > threshold), n_conc_down = sum(x < -threshold & y < -threshold),
    n_x_gt = sum(abs(x) > threshold), n_y_gt = sum(abs(y) > threshold),
    pct_concordant = 100 * sum(same) / max(1, sum(both))
  )
}

#' Gene table used for subsets and strata
#'
#' Builds one row per feature of a contrast table with an `expression` value (`ave_log2cpm` when available,
#' otherwise `log2(base_mean + 1)`), `chromosome` when known, and any extra logical flags.
#'
#' @param ct Canonical contrast table.
#' @param features `data.frame` with `feature_id` and optionally `chromosome` (e.g. the ingest rowData).
#' @param flags `data.frame` with `feature_id` and logical columns, or `NULL`.
#' @return `data.frame` with `feature_id`, `expression`, `chromosome` and flag columns.
#' @examples
#' ct <- new_contrast_table("c", c("a", "b"), c(1, 2),
#'   base_mean = c(10, 100), method = "m", source = "computed"
#' )
#' contrast_gene_table(ct)
#' @export
contrast_gene_table <- function(ct, features = NULL, flags = NULL) {
  ids <- unique(ct$feature_id)
  expr <- rep(NA_real_, length(ids))
  if ("ave_log2cpm" %in% names(ct)) {
    v <- ct[!is.na(ct$ave_log2cpm), c("feature_id", "ave_log2cpm")]
    expr <- v$ave_log2cpm[match(ids, v$feature_id)]
  }
  if (all(is.na(expr))) {
    v <- ct[!is.na(ct$base_mean), c("feature_id", "base_mean")]
    expr <- log2(v$base_mean[match(ids, v$feature_id)] + 1)
  }
  genes <- data.frame(feature_id = ids, expression = expr, stringsAsFactors = FALSE)
  genes$chromosome <- if (!is.null(features) && "chromosome" %in% names(features)) {
    features$chromosome[match(ids, features$feature_id)]
  } else {
    NA_character_
  }
  merge_gene_flags(genes, flags)
}

#' Concordance module
#'
#' For each pair of contrasts expected to share an effect, each estimator, gene subset and threshold, reports
#' correlation and sign concordance of the fold changes. This is a descriptive comparison: genes are not
#' independent and no p-value is attached.
#'
#' @param inputs `list` with `contrast_stats` (canonical table, all available contrasts), optional `features`
#'   and `gene_flags`.
#' @param params `list`. `modules.concordance` (`thresholds`, `subsets`).
#' @param context `list` with `cfg` (needs `contrast_pairs`) and `organism`.
#' @return `list(tables, figures, meta)`; table `concordance_summary`, figure `concordance_scatter`.
#' @examples
#' set.seed(1)
#' ids <- paste0("g", 1:50)
#' ct <- rbind(
#'   new_contrast_table("a", ids, rnorm(50), method = "m", source = "computed"),
#'   new_contrast_table("b", ids, rnorm(50), method = "m", source = "computed")
#' )
#' cfg <- list(contrast_pairs = list(list(id = "ab", first = "a", second = "b")))
#' run_concordance(
#'   list(contrast_stats = ct), list(thresholds = 1, subsets = "all"),
#'   list(cfg = cfg, organism = "human")
#' )$tables$concordance_summary
#' @export
run_concordance <- function(inputs, params, context) {
  ct <- inputs$contrast_stats
  pairs <- context$cfg$contrast_pairs %||% list()
  if (length(pairs) == 0) abort("concordance module: no contrast_pairs in the config")
  genes <- contrast_gene_table(ct, inputs$features, inputs$gene_flags)
  thresholds <- unlist(params$thresholds %||% c(0.58, 1))
  subsets <- unlist(params$subsets %||% "all")
  rows <- list()
  for (m in unique(ct$method)) {
    for (sub in subsets) {
      keep_ids <- genes$feature_id[gene_subset(genes, sub, context$organism)]
      for (thr in thresholds) {
        for (p in pairs) {
          x <- ct[ct$method == m & ct$contrast_id == p$first, c("feature_id", "lfc")]
          y <- ct[ct$method == m & ct$contrast_id == p$second, c("feature_id", "lfc")]
          if (nrow(x) == 0 || nrow(y) == 0) {
            abort(
              "concordance module: pair '", p$id, "' refers to contrast(s) without estimates for method '", m,
              "': ", p$first, ", ", p$second
            )
          }
          ids <- intersect(intersect(x$feature_id, y$feature_id), keep_ids)
          ids <- x$feature_id[x$feature_id %in% ids]
          met <- concordance_metrics(x$lfc[match(ids, x$feature_id)], y$lfc[match(ids, y$feature_id)], thr)
          rows[[length(rows) + 1]] <- data.frame(
            method = m, subset = sub, n_genes = length(ids), threshold = thr,
            pair_id = p$id, first = p$first, second = p$second, as.list(met)
          )
        }
      }
    }
  }
  summary <- do.call(rbind, rows)
  figs <- list()
  first_method <- unique(ct$method)[[1]]
  figs$concordance_scatter <- plot_concordance_scatter(
    ct[ct$method == first_method, ], pairs, thresholds[[length(thresholds)]]
  )
  list(
    tables = list(concordance_summary = summary),
    figures = figs,
    meta = list(concordance_summary = list(
      description = "Correlation and sign concordance of fold changes between contrast pairs (descriptive).",
      units = list(
        pearson = "correlation", spearman = "correlation", pct_concordant = "percent", threshold = "log2_ratio"
      )
    ))
  )
}

#' Scatter plot of two contrasts' fold changes
#'
#' @param ct Canonical contrast table of one method.
#' @param pairs `list` of contrast pairs (`id`, `first`, `second`).
#' @param threshold `numeric(1)`. Threshold drawn and used to colour concordant genes.
#' @return A `ggplot` object faceted by pair.
#' @examples
#' ct <- rbind(
#'   new_contrast_table("a", c("g1", "g2"), c(1, -2), method = "m", source = "computed"),
#'   new_contrast_table("b", c("g1", "g2"), c(2, 1), method = "m", source = "computed")
#' )
#' plot_concordance_scatter(ct, list(list(id = "ab", first = "a", second = "b")), 1)
#' @export
plot_concordance_scatter <- function(ct, pairs, threshold) {
  d <- do.call(rbind, lapply(pairs, function(p) {
    x <- ct[ct$contrast_id == p$first, c("feature_id", "lfc")]
    y <- ct[ct$contrast_id == p$second, c("feature_id", "lfc")]
    ids <- intersect(x$feature_id, y$feature_id)
    data.frame(pair = p$id, x = x$lfc[match(ids, x$feature_id)], y = y$lfc[match(ids, y$feature_id)])
  }))
  concordant <- abs(d$x) > threshold & abs(d$y) > threshold & sign(d$x) == sign(d$y)
  d$class <- factor(ifelse(concordant, "concordant", "other"), levels = c("other", "concordant"))
  d <- d[order(d$class), ]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$class)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#c9c8c3", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#c9c8c3", linewidth = 0.3) +
    ggplot2::geom_point(size = 0.6, alpha = 0.7) +
    ggplot2::scale_colour_manual(
      values = c(other = "#b5b4af", concordant = "#2a78d6"), drop = FALSE,
      name = paste0("|log2FC| > ", threshold, " in both, same sign")
    ) +
    ggplot2::facet_wrap(~pair) +
    ggplot2::labs(x = "log2 fold change (first contrast)", y = "log2 fold change (second contrast)") +
    theme_bulkrnaseq() +
    ggplot2::theme(legend.position = "bottom")
}
