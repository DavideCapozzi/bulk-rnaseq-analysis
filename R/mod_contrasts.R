#' Resolve the samples of a contrast side
#'
#' @param spec `character`. Sample IDs, or `"factor=level"` terms selecting all samples with that level.
#' @param samples `data.frame` of sample metadata with `sample_id`.
#' @return `character` sample IDs.
#' @examples
#' ss <- data.frame(sample_id = c("a", "b", "c"), g = c("x", "y", "x"))
#' contrast_samples("g=x", ss)
#' @export
contrast_samples <- function(spec, samples) {
  out <- unlist(lapply(unlist(spec), function(s) {
    if (grepl("=", s, fixed = TRUE)) {
      kv <- strsplit(s, "=", fixed = TRUE)[[1]]
      if (!kv[[1]] %in% names(samples)) abort("contrast term '", s, "': no sample sheet column '", kv[[1]], "'")
      as.character(samples$sample_id[as.character(samples[[kv[[1]]]]) == kv[[2]]])
    } else {
      if (!s %in% samples$sample_id) abort("contrast term '", s, "': unknown sample_id")
      s
    }
  }))
  if (length(out) == 0) abort("contrast side '", paste(spec, collapse = ","), "' selects no sample")
  unique(out)
}

#' Log-CPM difference estimator of fold changes
#'
#' A descriptive estimator: the difference between mean log2 CPM of the numerator and denominator samples,
#' after library normalisation and a prior count scaled by library size. It needs no dispersion model, so it
#' works for one-vs-one contrasts; it produces no p-value.
#'
#' @param se `SummarizedExperiment` with a `counts` assay.
#' @param contrasts `list` of contrasts, each with `id`, `numerator`, `denominator` (see [contrast_samples()]).
#' @param rows `logical` or `NULL`. Genes to use (e.g. from [filter_genes()]); normalisation is computed on them.
#' @param prior_count `numeric(1)`. Prior count (counts units, scaled by library size).
#' @param norm_method `character(1)`. `"TMM"`, `"RLE"`, `"upperquartile"`, `"none"` or `"reference"`.
#' @param reference_genes `character`. Symbols used when `norm_method = "reference"`.
#' @param norm_factors `numeric` or `NULL`. Explicit normalisation factors (override `norm_method`).
#' @return Canonical contrast table with method `"logcpm_difference"` and extension column `ave_log2cpm`
#'   (mean log2 CPM over all samples). `df_residual` is the number of samples in the contrast minus 2.
#' @examples
#' set.seed(1)
#' m <- matrix(rpois(40, 50), 10, dimnames = list(paste0("g", 1:10), paste0("s", 1:4)))
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   colData = S4Vectors::DataFrame(sample_id = colnames(m), row.names = colnames(m))
#' )
#' head(estimate_lfc(se, list(list(id = "s1_vs_s2", numerator = "s1", denominator = "s2"))))
#' @export
estimate_lfc <- function(se, contrasts, rows = NULL, prior_count = 2,
                         norm_method = c("TMM", "RLE", "upperquartile", "none", "reference"),
                         reference_genes = character(), norm_factors = NULL) {
  norm_method <- match.arg(norm_method)
  counts <- SummarizedExperiment::assay(se, "counts")
  if (!is.null(rows)) counts <- counts[rows, , drop = FALSE]
  symbols <- SummarizedExperiment::rowData(se)$gene_symbol
  if (!is.null(rows)) symbols <- symbols[rows]
  if (is.null(norm_factors) && norm_method == "reference") {
    ref <- symbols %in% reference_genes
    if (!any(ref)) {
      abort("estimate_lfc(): norm_method 'reference' but none of the reference genes is among the selected genes")
    }
    norm_factors <- reference_norm_factors(counts, ref)
  }
  lc <- log_cpm(counts,
    norm_method = if (norm_method == "reference") "none" else norm_method,
    prior_count = prior_count, norm_factors = norm_factors
  )
  samples <- as.data.frame(SummarizedExperiment::colData(se))
  ave <- rowMeans(lc)
  out <- lapply(contrasts, function(cn) {
    num <- contrast_samples(cn$numerator, samples)
    den <- contrast_samples(cn$denominator, samples)
    lfc <- rowMeans(lc[, num, drop = FALSE]) - rowMeans(lc[, den, drop = FALSE])
    new_contrast_table(cn$id, rownames(lc), lfc,
      gene_symbol = symbols, significance_metric = "none",
      method = "logcpm_difference", df_residual = length(num) + length(den) - 2L,
      source = "computed", extra = data.frame(ave_log2cpm = unname(ave))
    )
  })
  ct <- do.call(rbind, out)
  rownames(ct) <- NULL
  attr(ct, "norm_factors") <- if (is.null(norm_factors)) {
    edgeR::calcNormFactors(edgeR::DGEList(counts), method = norm_method)$samples$norm.factors
  } else {
    norm_factors
  }
  ct
}

#' Shrunken fold changes for single-sample contrasts in a two-factor design
#'
#' Dispersions are estimated once on an additive design (which must have at least 1 residual degree of
#' freedom; with exactly 1 the result is exploratory), then transferred to a nested design
#' `~ other + other:factor` in which each one-vs-one contrast is a coefficient, and shrunk with apeglm.
#' No p-value is reported.
#'
#' @inheritParams estimate_lfc
#' @param design `formula`. Additive design used for size factors and dispersions, e.g. `~ batch + condition`.
#' @return Canonical contrast table with method `"deseq2_apeglm"`, `lfc_se` and `base_mean`; `df_residual` is the
#'   residual df of `design`, and extension column and attribute `exploratory` are `TRUE` when it equals 1.
#' @examples
#' \donttest{
#' if (requireNamespace("apeglm", quietly = TRUE)) {
#'   set.seed(1)
#'   m <- matrix(rnbinom(400, mu = 100, size = 10), 100,
#'     dimnames = list(paste0("g", 1:100), c("a1", "a2", "b1", "b2"))
#'   )
#'   cd <- S4Vectors::DataFrame(
#'     sample_id = colnames(m), cond = factor(c("c", "t", "c", "t")),
#'     batch = factor(c("x", "x", "y", "y")), row.names = colnames(m)
#'   )
#'   se <- SummarizedExperiment::SummarizedExperiment(list(counts = m), colData = cd)
#'   ct <- suppressWarnings(estimate_lfc_apeglm(se, list(list(
#'     id = "t_vs_c_x", numerator = "a2",
#'     denominator = "a1"
#'   )), design = ~ batch + cond))
#'   head(ct)
#' }
#' }
#' @export
estimate_lfc_apeglm <- function(se, contrasts, design, rows = NULL) {
  if (!requireNamespace("apeglm", quietly = TRUE)) abort("estimate_lfc_apeglm() needs the 'apeglm' package")
  counts <- SummarizedExperiment::assay(se, "counts")
  symbols <- SummarizedExperiment::rowData(se)$gene_symbol
  if (!is.null(rows)) {
    counts <- counts[rows, , drop = FALSE]
    symbols <- symbols[rows]
  }
  factors <- all.vars(design)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  assert_columns(cd, factors, "estimate_lfc_apeglm() sample metadata")
  cd <- cd[, factors, drop = FALSE]
  for (f in factors) if (!is.factor(cd[[f]])) cd[[f]] <- factor(cd[[f]])
  guard <- assert_df(design, cd, "DESeq2 dispersion estimate")
  dds_add <- DESeq2::DESeqDataSetFromMatrix(counts, cd, design = design)
  dds_add <- DESeq2::estimateSizeFactors(dds_add)
  dds_add <- DESeq2::estimateDispersions(dds_add, quiet = TRUE)
  samples <- data.frame(sample_id = colnames(counts), cd, check.names = FALSE)
  specs <- lapply(contrasts, function(cn) nested_coefficient(cn, samples, factors))
  fits <- list()
  out <- lapply(seq_along(contrasts), function(k) {
    sp <- specs[[k]]
    key <- paste(sp$other, sp$factor)
    if (is.null(fits[[key]])) {
      nested <- stats::as.formula(paste("~", sp$other, "+", sp$other, ":", sp$factor))
      dd <- DESeq2::DESeqDataSetFromMatrix(counts, cd, design = nested)
      DESeq2::sizeFactors(dd) <- DESeq2::sizeFactors(dds_add)
      DESeq2::dispersions(dd) <- DESeq2::dispersions(dds_add)
      fits[[key]] <<- DESeq2::nbinomWaldTest(dd, quiet = TRUE)
    }
    dd <- fits[[key]]
    if (!sp$coef %in% DESeq2::resultsNames(dd)) {
      abort("estimate_lfc_apeglm(): coefficient '", sp$coef, "' not in the nested model")
    }
    res <- DESeq2::lfcShrink(dd, coef = sp$coef, type = "apeglm", quiet = TRUE)
    new_contrast_table(contrasts[[k]]$id, rownames(counts), res$log2FoldChange,
      gene_symbol = symbols,
      lfc_se = res$lfcSE, base_mean = res$baseMean, significance_metric = "none",
      method = "deseq2_apeglm", df_residual = guard$df, source = "computed",
      extra = data.frame(exploratory = rep(guard$exploratory, nrow(counts)))
    )
  })
  ct <- do.call(rbind, out)
  rownames(ct) <- NULL
  attr(ct, "exploratory") <- guard$exploratory
  ct
}

# For a one-vs-one contrast in a two-factor design, the factor that differs between the two samples and the
# coefficient name of the nested model "~ other + other:factor".
nested_coefficient <- function(cn, samples, factors) {
  num <- contrast_samples(cn$numerator, samples)
  den <- contrast_samples(cn$denominator, samples)
  if (length(num) != 1 || length(den) != 1) {
    abort("contrast '", cn$id, "': deseq2_apeglm supports single-sample sides only")
  }
  if (length(factors) != 2) abort("contrast '", cn$id, "': deseq2_apeglm needs a design with exactly two factors")
  a <- samples[samples$sample_id == num, factors]
  b <- samples[samples$sample_id == den, factors]
  differs <- factors[vapply(factors, function(f) as.character(a[[f]]) != as.character(b[[f]]), logical(1))]
  if (length(differs) != 1) abort("contrast '", cn$id, "': the two samples must differ in exactly one factor")
  other <- setdiff(factors, differs)
  ref <- levels(samples[[differs]])[[1]]
  if (as.character(b[[differs]]) != ref) {
    abort("contrast '", cn$id, "': the denominator must be at the reference level '", ref, "' of '", differs, "'")
  }
  list(
    factor = differs, other = other,
    coef = paste0(other, as.character(a[[other]]), ".", differs, as.character(a[[differs]]))
  )
}

#' Spread of fold changes per contrast and gene subset
#'
#' @param ct Canonical contrast table (one or more methods).
#' @param genes `data.frame` of gene flags for [gene_subset()] (`feature_id`, `expression`, flags).
#' @param subsets `character`. Subset expressions.
#' @param threshold `numeric(1)`. |lfc| threshold counted in `n_abs_lfc_gt_threshold`.
#' @param organism `character(1)`.
#' @return `data.frame` with `method`, `subset`, `contrast_id`, `sd_lfc`, `n_abs_lfc_gt_threshold`.
#' @examples
#' ct <- new_contrast_table("c", c("a", "b", "c"), c(-1, 0, 2), method = "m", source = "computed")
#' contrast_spread(ct, data.frame(feature_id = c("a", "b", "c"), expression = 1:3), "all")
#' @export
contrast_spread <- function(ct, genes, subsets = "all", threshold = 1, organism = "human") {
  rows <- list()
  for (sub in subsets) {
    keep_ids <- genes$feature_id[gene_subset(genes, sub, organism)]
    for (m in unique(ct$method)) {
      for (cid in unique(ct$contrast_id)) {
        x <- ct$lfc[ct$method == m & ct$contrast_id == cid & ct$feature_id %in% keep_ids]
        rows[[length(rows) + 1]] <- data.frame(
          method = m, subset = sub, contrast_id = cid, sd_lfc = stats::sd(x),
          n_abs_lfc_gt_threshold = sum(abs(x) > threshold)
        )
      }
    }
  }
  do.call(rbind, rows)
}

#' Contrasts module
#'
#' Filters genes, then computes the configured fold-change estimators for every contrast in the config.
#' Estimators are descriptive (no p-values), so the module runs on one-vs-one contrasts.
#'
#' @param inputs `list` with `se` (SummarizedExperiment with counts) and optional `gene_flags`.
#' @param params `list`. `modules.contrasts` from the config.
#' @param context `list` with `cfg` (full config), `gene_lists`, `organism`.
#' @return `list(tables, figures, meta)`; tables `contrast_stats`, `contrast_spread`, `normalization_factors`.
#' @examples
#' set.seed(1)
#' m <- matrix(rpois(400, 50), 100, dimnames = list(paste0("g", 1:100), paste0("s", 1:4)))
#' rd <- S4Vectors::DataFrame(
#'   feature_id = rownames(m), gene_symbol = rownames(m), biotype = "protein_coding",
#'   chromosome = "1"
#' )
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   rowData = rd,
#'   colData = S4Vectors::DataFrame(sample_id = colnames(m), row.names = colnames(m))
#' )
#' cfg <- list(
#'   contrasts = list(list(id = "s1_vs_s2", numerator = "s1", denominator = "s2")),
#'   gene_filter = list(min_cpm = 1, min_samples = 1, sex_genes = "keep")
#' )
#' res <- run_contrasts(
#'   list(se = se), list(estimators = "logcpm_difference"),
#'   list(cfg = cfg, gene_lists = list(), organism = "human")
#' )
#' names(res$tables)
#' @export
run_contrasts <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "contrasts"))
  cfg <- context$cfg
  se <- inputs$se
  gf <- with_defaults(cfg$gene_filter, "gene_filter")
  keep <- config_filter_genes(se, cfg, context)
  contrasts <- cfg$contrasts
  if (length(contrasts) == 0) abort("contrasts module: no contrasts defined in the config")
  estimators <- unlist(params$estimators)
  pieces <- list()
  nf_rows <- list()
  if ("logcpm_difference" %in% estimators) {
    ct_a <- estimate_lfc(se, contrasts,
      rows = keep, prior_count = params$prior_count, norm_method = params$norm_method,
      reference_genes = resolve_gene_list(params$reference_genes, context$gene_lists)
    )
    nf_rows[[1]] <- data.frame(
      norm_variant = params$norm_method, sample_id = colnames(se),
      norm_factor = attr(ct_a, "norm_factors")
    )
    pieces$logcpm_difference <- ct_a
  }
  exploratory <- NA
  if ("deseq2_apeglm" %in% estimators) {
    d <- find_design(cfg, params$dispersion_design)
    ct_b <- estimate_lfc_apeglm(se, contrasts, design = stats::as.formula(d$formula), rows = keep)
    exploratory <- attr(ct_b, "exploratory")
    pieces$deseq2_apeglm <- ct_b
  }
  ct <- bind_contrast_tables(pieces)
  rd <- as.data.frame(SummarizedExperiment::rowData(se))[keep, , drop = FALSE]
  expr <- if (!is.null(pieces$logcpm_difference)) {
    pieces$logcpm_difference$ave_log2cpm[match(rd$feature_id, pieces$logcpm_difference$feature_id)]
  } else {
    NA_real_
  }
  genes <- merge_gene_flags(
    data.frame(feature_id = rd$feature_id, expression = expr, chromosome = rd$chromosome),
    inputs$gene_flags
  )
  spread <- contrast_spread(
    ct, genes, unlist(params$spread_subsets), params$spread_threshold,
    context$organism
  )
  list(
    tables = list(contrast_stats = ct, contrast_spread = spread, normalization_factors = do.call(rbind, nf_rows)),
    figures = list(),
    meta = list(
      contrast_stats = list(
        description = "Fold-change estimates per contrast and estimator (descriptive: no p-values).",
        units = list(
          lfc = "log2_ratio",
          ave_log2cpm = paste0("log2_cpm_", tolower(params$norm_method), "_prior", params$prior_count),
          base_mean = "normalized_counts"
        ),
        filters = list(gene_filter = gf, n_genes = sum(keep)), exploratory = exploratory
      ),
      contrast_spread = list(units = list(sd_lfc = "log2_ratio")),
      normalization_factors = list(units = list(norm_factor = "dimensionless"))
    )
  )
}

find_design <- function(cfg, id) {
  if (is.null(id)) abort("contrasts module: deseq2_apeglm needs modules.contrasts.dispersion_design")
  for (d in cfg$designs %||% list()) {
    if (identical(d$id, id)) {
      return(d)
    }
  }
  abort("contrasts module: design '", id, "' not found in designs")
}

merge_gene_flags <- function(genes, flags) {
  if (is.null(flags)) {
    return(genes)
  }
  extra <- setdiff(names(flags), names(genes))
  idx <- match(genes$feature_id, flags$feature_id)
  for (col in extra) genes[[col]] <- flags[[col]][idx]
  genes
}
