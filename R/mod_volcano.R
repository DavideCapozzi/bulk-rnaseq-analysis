#' Classify genes as up, down or not significant
#'
#' Uses the contrast's own significance metric: adjusted p-value at most `padj` or posterior probability at
#' least `posterior_prob`, combined with a linear fold-change threshold.
#'
#' @param ct Canonical contrast table of one contrast.
#' @param thresholds `list` with `padj`, `posterior_prob`, `abs_fold_change` (linear, e.g. 1.5).
#' @return `character` vector: `"UP"`, `"DOWN"` or `"NS"`.
#' @examples
#' ct <- new_contrast_table("c", c("a", "b", "c"), c(1, -1, 2),
#'   padj = c(0.01, 0.01, 0.5),
#'   significance_metric = "padj", method = "m", source = "imported:x",
#'   inference_scope = "sample"
#' )
#' classify_de(ct, list(padj = 0.05, abs_fold_change = 1.5))
#' @export
classify_de <- function(ct, thresholds) {
  thresholds <- with_defaults(thresholds, "thresholds")
  sig <- significance_values(ct)
  if (sig$metric == "none") abort("classify_de(): contrast '", ct$contrast_id[[1]], "' has no significance metric")
  passes <- if (sig$metric == "padj") {
    sig$values <= thresholds$padj
  } else {
    sig$values >= thresholds$posterior_prob
  }
  lfc_thr <- log2(thresholds$abs_fold_change)
  out <- ifelse(passes & ct$lfc >= lfc_thr, "UP", ifelse(passes & ct$lfc <= -lfc_thr, "DOWN", "NS"))
  out[is.na(out)] <- "NS"
  out
}

#' Volcano module
#'
#' One volcano plot per contrast with a significance metric. The y axis is labelled with the metric's own name
#' (`-log10(padj)` or posterior probability), never with a relabelled quantity.
#'
#' @param inputs `list` with `contrast_stats`.
#' @param params `list`. `modules.volcano` (`label_top`).
#' @param context `list` with `cfg` (for `thresholds`).
#' @return `list(tables, figures, meta)`; table `volcano_points`, one figure per contrast.
#' @examples
#' ct <- new_contrast_table("c", paste0("g", 1:3), c(1, -1, 2),
#'   padj = c(0.01, 0.01, 0.5),
#'   significance_metric = "padj", method = "m", source = "imported:x",
#'   inference_scope = "sample"
#' )
#' res <- run_volcano(list(contrast_stats = ct), list(label_top = 2), list(cfg = list()))
#' res$tables$volcano_points
#' @export
run_volcano <- function(inputs, params, context) {
  ct <- inputs$contrast_stats
  ct <- ct[ct$significance_metric != "none", , drop = FALSE]
  if (nrow(ct) == 0) abort("volcano module: no contrast with a significance metric")
  params <- with_defaults(params, c("modules", "volcano"))
  thresholds <- with_defaults(context$cfg$thresholds, "thresholds")
  pts <- list()
  figs <- list()
  for (key in unique(paste(ct$contrast_id, ct$method, sep = "\r"))) {
    x <- ct[paste(ct$contrast_id, ct$method, sep = "\r") == key, , drop = FALSE]
    sig <- significance_values(x)
    d <- data.frame(
      contrast_id = x$contrast_id, method = x$method, feature_id = x$feature_id,
      gene_symbol = x$gene_symbol, lfc = x$lfc, significance_metric = sig$metric,
      significance = sig$values, class = classify_de(x, thresholds),
      inference_scope = x$inference_scope, exploratory = contrast_exploratory(x), stringsAsFactors = FALSE
    )
    pts[[key]] <- d
    figs[[paste0("volcano_", x$contrast_id[[1]])]] <- plot_volcano(d, thresholds, params$label_top)
  }
  tab <- do.call(rbind, unname(pts))
  rownames(tab) <- NULL
  list(
    tables = list(volcano_points = tab), figures = figs,
    meta = list(volcano_points = list(
      description = "Genes classified by significance metric and linear fold-change threshold.",
      units = list(lfc = "log2_ratio", significance = "as named in significance_metric"),
      thresholds = thresholds
    ))
  )
}

#' Draw a volcano plot
#'
#' @param d `data.frame` from [run_volcano()] for one contrast.
#' @param thresholds `list` of thresholds.
#' @param label_top `integer(1)`. Number of most significant UP/DOWN genes labelled.
#' @return A `ggplot` object.
#' @examples
#' d <- data.frame(
#'   gene_symbol = c("A", "B"), lfc = c(2, -2), significance = c(0.001, 0.01),
#'   significance_metric = "padj", class = c("UP", "DOWN"), contrast_id = "c"
#' )
#' plot_volcano(d, list(padj = 0.05, abs_fold_change = 1.5), 2)
#' @export
plot_volcano <- function(d, thresholds, label_top = 15) {
  thresholds <- with_defaults(thresholds, "thresholds")
  metric <- d$significance_metric[[1]]
  if (metric == "padj") {
    d$y <- -log10(d$significance)
    ylab <- "-log10(adjusted p-value)"
    yline <- -log10(thresholds$padj)
  } else {
    d$y <- d$significance
    ylab <- "posterior probability"
    yline <- thresholds$posterior_prob
  }
  d$class <- factor(d$class, levels = c("NS", "DOWN", "UP"))
  d <- d[order(d$class), ]
  sig <- d[d$class != "NS", , drop = FALSE]
  lab <- utils::head(sig[order(-sig$y, -abs(sig$lfc)), , drop = FALSE], label_top)
  lfc_thr <- log2(thresholds$abs_fold_change)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$lfc, y = .data$y, colour = .data$class)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.7) +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = yline, linetype = "dashed", linewidth = 0.3) +
    ggrepel::geom_text_repel(
      data = lab, ggplot2::aes(label = .data$gene_symbol), size = 2.8, colour = "#1b1b1b",
      seed = 1, max.overlaps = 30, min.segment.length = 0
    ) +
    ggplot2::scale_colour_manual(
      values = c(NS = "#c9c8c3", DOWN = "#2a78d6", UP = "#e34948"), drop = FALSE, name = NULL
    ) +
    ggplot2::labs(title = d$contrast_id[[1]], x = "log2 fold change", y = ylab) +
    theme_bulkrnaseq()
}
