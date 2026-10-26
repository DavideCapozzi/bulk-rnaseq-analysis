#' DESeq2 Wald test for one design
#'
#' Keeps genes with at least `min_count` reads in at least `min_samples` samples, fits the design with
#' `DESeq2::DESeq()` and extracts the contrast with `DESeq2::results()` (maximum-likelihood fold changes, not
#' shrunken). The design must have residual degrees of freedom: [assert_df()] stops at 0 and flags 1 as
#' exploratory.
#'
#' @param se `SummarizedExperiment` with a `counts` assay and the design variables in colData.
#' @param design `list` with `id`, `formula` and `contrast` (`c(factor, numerator, denominator)`), as in the
#'   config's `designs`.
#' @param min_count,min_samples Pre-filter: genes with at least `min_count` reads in at least `min_samples`
#'   samples are tested.
#' @param alpha `numeric(1)`. `alpha` of `DESeq2::results()`, the FDR target of independent filtering.
#' @return Canonical contrast table (method `"deseq2_wald"`, `inference_scope = "sample"`) with attributes
#'   `df_residual`, `exploratory`, `design_rank`, `n_samples`.
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- matrix(rnbinom(6000, mu = 100, size = 10), 1000,
#'   dimnames = list(paste0("g", 1:1000), paste0("s", 1:6))
#' )
#' cd <- S4Vectors::DataFrame(sample_id = colnames(m), g = rep(c("a", "b"), each = 3), row.names = colnames(m))
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   rowData = S4Vectors::DataFrame(feature_id = rownames(m), gene_symbol = rownames(m)), colData = cd
#' )
#' ct <- deseq2_wald(se, list(id = "b_vs_a", formula = "~ g", contrast = c("g", "b", "a")))
#' head(ct)
#' }
#' @export
deseq2_wald <- function(se, design, min_count = 10, min_samples = 3, alpha = 0.1) {
  counts <- SummarizedExperiment::assay(se, "counts")
  f <- stats::as.formula(design$formula)
  vars <- all.vars(f)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  assert_columns(cd, vars, paste0("de_test design '", design$id, "' sample metadata"))
  cd <- cd[, vars, drop = FALSE]
  for (v in vars) if (!is.factor(cd[[v]])) cd[[v]] <- factor(cd[[v]])
  contrast <- as.character(unlist(design$contrast))
  if (length(contrast) != 3 || !contrast[[1]] %in% vars) {
    abort("de_test design '", design$id, "': contrast must be c(<factor in the formula>, <numerator>, <denominator>)")
  }
  bad <- setdiff(contrast[2:3], levels(cd[[contrast[[1]]]]))
  if (length(bad) > 0) abort("de_test design '", design$id, "': unknown level(s) ", paste(bad, collapse = ", "))
  guard <- assert_df(f, cd, paste0("de_test design '", design$id, "'"))
  keep <- rowSums(counts >= min_count) >= min_samples
  if (!any(keep)) abort("de_test design '", design$id, "': no gene passes the pre-filter")
  dds <- DESeq2::DESeqDataSetFromMatrix(counts[keep, , drop = FALSE], cd, design = f)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  ct <- new_contrast_table(design$id, rownames(res), res$log2FoldChange,
    gene_symbol = SummarizedExperiment::rowData(se)$gene_symbol[keep], lfc_se = res$lfcSE, stat = res$stat,
    pvalue = res$pvalue, padj = res$padj, base_mean = res$baseMean, significance_metric = "padj",
    method = "deseq2_wald", df_residual = guard$df, source = "computed", inference_scope = "sample"
  )
  attr(ct, "exploratory") <- guard$exploratory
  attr(ct, "design_rank") <- as.integer(qr(stats::model.matrix(f, cd))$rank)
  attr(ct, "n_samples") <- ncol(counts)
  ct
}

#' Compare a computed contrast with an imported one
#'
#' Matches genes by `feature_id` and reports agreement of fold changes, p-values and up/down calls, e.g. to
#' check that a toolkit fit reproduces a provider's table before using it.
#'
#' @param computed,imported Canonical contrast tables of one contrast each, with significance metric `padj`.
#' @param thresholds `list` of thresholds for [classify_de()].
#' @return One-row `data.frame`.
#' @examples
#' a <- new_contrast_table("a", c("g1", "g2"), c(1, -1),
#'   pvalue = c(0.01, 0.2), padj = c(0.02, 0.3),
#'   significance_metric = "padj", method = "m", source = "computed", inference_scope = "sample"
#' )
#' b <- transform(a, contrast_id = "b", source = "imported:x")
#' compare_contrasts(a, b, list(padj = 0.05, abs_fold_change = 1.5))
#' @export
compare_contrasts <- function(computed, imported, thresholds) {
  i <- match(imported$feature_id, computed$feature_id)
  shared <- !is.na(i)
  a <- computed[i[shared], , drop = FALSE]
  b <- imported[shared, , drop = FALSE]
  ca <- classify_de(a, thresholds)
  cb <- classify_de(b, thresholds)
  rel <- function(x, y) {
    ok <- !is.na(x) & !is.na(y) & y != 0
    if (any(ok)) max(abs(x[ok] - y[ok]) / abs(y[ok])) else NA_real_
  }
  absd <- function(x, y) if (any(!is.na(x - y))) max(abs(x - y), na.rm = TRUE) else NA_real_
  data.frame(
    computed = computed$contrast_id[[1]], imported = imported$contrast_id[[1]],
    n_computed = nrow(computed), n_imported = nrow(imported), n_shared = sum(shared),
    lfc_pearson = suppressWarnings(stats::cor(a$lfc, b$lfc, use = "complete.obs")),
    lfc_max_abs_diff = absd(a$lfc, b$lfc), pvalue_max_rel_diff = rel(a$pvalue, b$pvalue),
    padj_max_abs_diff = absd(a$padj, b$padj), padj_na_computed = sum(is.na(a$padj)),
    n_same_class = sum(ca == cb), n_up_both = sum(ca == "UP" & cb == "UP"),
    n_down_both = sum(ca == "DOWN" & cb == "DOWN"), n_sig_only_computed = sum(ca != "NS" & cb == "NS"),
    n_sig_only_imported = sum(ca == "NS" & cb != "NS"), n_opposite = sum(ca != "NS" & cb != "NS" & ca != cb),
    stringsAsFactors = FALSE
  )
}

#' Inferential differential-expression module
#'
#' Runs a sample-level test for every design of the config that declares a `contrast`. Only engine
#' `deseq2` is implemented. Results feed volcano and ORA as any other contrast, with
#' `inference_scope = "sample"`.
#'
#' @param inputs `list` with `se` (counts) and optional `contrast_stats` (imported tables, for
#'   `compare_with`).
#' @param params `list`. `modules.de_test` from the config.
#' @param context `list` with `cfg`.
#' @return `list(tables, figures, meta)`; tables `contrast_stats`, `design_summary` and, when
#'   `compare_with` is set, `provider_comparison`.
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- matrix(rnbinom(6000, mu = 100, size = 10), 1000,
#'   dimnames = list(paste0("g", 1:1000), paste0("s", 1:6))
#' )
#' cd <- S4Vectors::DataFrame(sample_id = colnames(m), g = rep(c("a", "b"), each = 3), row.names = colnames(m))
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   rowData = S4Vectors::DataFrame(feature_id = rownames(m), gene_symbol = rownames(m)), colData = cd
#' )
#' cfg <- list(designs = list(list(id = "b_vs_a", formula = "~ g", contrast = c("g", "b", "a"))))
#' run_de_test(list(se = se), list(), list(cfg = cfg))$tables$design_summary
#' }
#' @export
run_de_test <- function(inputs, params, context) {
  cfg <- context$cfg
  designs <- Filter(function(d) !is.null(d$contrast), cfg$designs %||% list())
  if (length(designs) == 0) abort("de_test module: no design with a contrast in the config (designs[].contrast)")
  thresholds <- cfg$thresholds %||% list(padj = 0.05, abs_fold_change = 1.5)
  pf <- params$prefilter %||% list()
  alpha <- params$results_alpha %||% 0.1
  pieces <- list()
  summary_rows <- list()
  for (d in designs) {
    engine <- d$engine %||% "deseq2"
    if (engine != "deseq2") abort("de_test design '", d$id, "': engine '", engine, "' is not implemented (use deseq2)")
    ct <- deseq2_wald(inputs$se, d, pf$min_count %||% 10, pf$min_samples %||% 3, alpha)
    cls <- classify_de(ct, thresholds)
    summary_rows[[d$id]] <- data.frame(
      design_id = d$id, formula = d$formula, engine = engine, contrast = paste(unlist(d$contrast), collapse = ","),
      n_samples = attr(ct, "n_samples"), design_rank = attr(ct, "design_rank"),
      df_residual = unique(ct$df_residual), exploratory = attr(ct, "exploratory"), n_genes_tested = nrow(ct),
      n_padj_na = sum(is.na(ct$padj)), n_up = sum(cls == "UP"), n_down = sum(cls == "DOWN"),
      stringsAsFactors = FALSE
    )
    pieces[[d$id]] <- ct
  }
  ct <- bind_contrast_tables(pieces)
  tables <- list(contrast_stats = ct, design_summary = do.call(rbind, unname(summary_rows)))
  cmp <- list()
  for (w in params$compare_with %||% list()) {
    imp <- inputs$contrast_stats
    imp <- if (!is.null(imp)) imp[imp$contrast_id == w$contrast_table, , drop = FALSE]
    if (is.null(imp) || nrow(imp) == 0) {
      abort("de_test compare_with: imported contrast table '", w$contrast_table, "' not found")
    }
    if (!w$design %in% names(pieces)) abort("de_test compare_with: design '", w$design, "' has no contrast")
    cmp[[length(cmp) + 1]] <- compare_contrasts(pieces[[w$design]], imp, thresholds)
  }
  if (length(cmp)) tables$provider_comparison <- do.call(rbind, cmp)
  list(
    tables = tables, figures = list(),
    meta = list(
      contrast_stats = list(
        description = "DESeq2 Wald tests (maximum-likelihood log2 fold changes, not shrunken).",
        units = list(lfc = "log2_ratio", base_mean = "normalized_counts"),
        filters = list(prefilter = pf, results_alpha = alpha), inference_scope = "sample"
      ),
      design_summary = list(
        description = "Designs, residual degrees of freedom (n samples - rank) and up/down counts.",
        thresholds = thresholds
      ),
      provider_comparison = list(
        description = "Agreement between computed and imported contrasts on shared genes.",
        thresholds = thresholds
      )
    )
  )
}
