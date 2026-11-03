#' Reference axis module
#'
#' Relates the fold changes of each contrast (one estimator) to axes built from external reference profiles, and
#' summarises them on gene sets and contrast pairs. All measures are descriptive.
#' * `spearman_lfc_vs_axis`: per axis and contrast, Spearman correlation between the log2 fold change and the axis
#'   (`log2((mean numerator + pseudocount) / (mean denominator + pseudocount))` of HPA cell types) over the genes of
#'   `subset` that have an axis value.
#' * `mean_lfc_in_set`: per gene set and contrast, mean log2 fold change of the genes of `set_subset` in the set.
#' * `sign_concordance`: per contrast pair, among genes of `subset` with |log2 fold change| above
#'   `sign_test.threshold` in both contrasts, the fraction with the same sign and a two-sided binomial p-value
#'   (`inference_scope = "conditional_on_contrast"`: genes are not independent).
#' * `sample_axis_score`: per axis and sample, mean log CPM of the axis-high genes of `subset` (axis >= threshold)
#'   minus that of the axis-low genes (axis <= -threshold); `axis_score_shift`: per axis and configured contrast,
#'   mean score of the numerator samples minus that of the denominator samples. Both need counts.
#'
#' @param inputs `list` with `contrast_stats`, optional `se` (counts), `features` and `gene_flags`.
#' @param params `list`. `modules.reference_axis` from the config.
#' @param context `list` with `cfg`, `organism` and `resources` (paths by key).
#' @return `list(tables, figures, meta)`; tables `reference_axis` (one row per measure) and `axis_genes` (axis
#'   value and class of the genes of `subset`, per axis).
#' @examples
#' \donttest{
#' # called by the workflow on the project's contrast tables and HPA exports
#' }
#' @export
run_reference_axis <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "reference_axis"))
  ct <- inputs$contrast_stats
  if (is.null(ct) || nrow(ct) == 0) abort("reference_axis: no contrast tables")
  sel <- select_contrasts(ct, params$estimator, unlist(params$contrasts), context$cfg, "reference_axis")
  ct <- sel$ct
  contrasts <- sel$contrasts
  organism <- context$organism %||% "human"
  genes <- contrast_gene_table(ct, inputs$features, inputs$gene_flags)
  in_subset <- genes$feature_id[gene_subset(genes, params$subset, organism)]
  in_set_subset <- genes$feature_id[gene_subset(genes, params$set_subset, organism)]
  lfc_of <- function(k, ids) {
    x <- ct[ct$contrast_id == k, c("feature_id", "lfc")]
    x$lfc[match(ids, x$feature_id)]
  }
  axis_ids <- vapply(params$axes, `[[`, "", "id")
  if (anyDuplicated(axis_ids)) {
    abort("reference_axis: duplicated axis id(s): ", toString(unique(axis_ids[duplicated(axis_ids)])))
  }
  axes <- reference_axis_values(params$axes, genes$feature_id, context)
  rows <- list()
  for (a in params$axes) {
    v <- axes[match(in_subset, genes$feature_id), a$id]
    for (k in contrasts) {
      lfc <- lfc_of(k, in_subset)
      ok <- is.finite(v) & is.finite(lfc)
      rows[[length(rows) + 1]] <- reference_axis_row("spearman_lfc_vs_axis",
        axis_id = a$id, contrast_id = k, n = sum(ok),
        value = if (sum(ok) > 1) stats::cor(lfc[ok], v[ok], method = "spearman") else NA_real_
      )
    }
  }
  sets <- unlist(params$gene_sets$sets)
  if (length(sets) > 0) {
    key <- params$gene_sets$resource
    if (is.null(key)) abort("reference_axis: gene_sets.sets needs gene_sets.resource (an msigdb_export)")
    path <- context$resources[[key]]
    if (is.null(path)) abort("reference_axis: resource '", key, "' is not available in the cache")
    gs <- read_gene_sets(path, sets)
    by_symbol <- identical(context$cfg$annotation$id_type, "gene_symbol")
    member <- if (by_symbol) gs$gene_symbol else gs$ensembl_gene
    for (s in sets) {
      ids <- in_set_subset[in_set_subset %in% member[gs$gs_name == s]]
      for (k in contrasts) {
        lfc <- lfc_of(k, ids)
        lfc <- lfc[is.finite(lfc)]
        rows[[length(rows) + 1]] <- reference_axis_row("mean_lfc_in_set",
          set_id = s, contrast_id = k, n = length(lfc),
          value = if (length(lfc)) mean(lfc) else NA_real_
        )
      }
    }
  }
  thr <- params$sign_test$threshold
  pairs <- Filter(function(p) all(c(p$first, p$second) %in% contrasts), context$cfg$contrast_pairs %||% list())
  if (!is.null(thr)) {
    if (length(pairs) == 0) abort("reference_axis: sign_test needs contrast_pairs whose contrasts have estimates")
    for (p in pairs) {
      x <- lfc_of(p$first, in_subset)
      y <- lfc_of(p$second, in_subset)
      both <- is.finite(x) & is.finite(y) & abs(x) > thr & abs(y) > thr
      n <- sum(both)
      k <- sum(both & sign(x) == sign(y))
      rows[[length(rows) + 1]] <- reference_axis_row("sign_concordance",
        pair_id = p$id, n = n, n_concordant = k, value = if (n > 0) k / n else NA_real_,
        pvalue = if (n > 0) stats::binom.test(k, n)$p.value else NA_real_, inference_scope = "conditional_on_contrast"
      )
    }
  }
  axis_genes <- reference_axis_genes(params$axes, axes, genes, in_subset, ct)
  if (isTRUE(params$sample_scores) && length(params$axes) > 0 && !is.null(inputs$se)) {
    scores <- reference_axis_scores(inputs$se, params, axis_genes, context)
    rows <- c(rows, scores)
  }
  table <- do.call(rbind, c(list(reference_axis_row(NA_character_)[0, ]), rows))
  rownames(table) <- NULL
  list(
    tables = list(reference_axis = table, axis_genes = axis_genes),
    figures = list(),
    meta = list(
      reference_axis = list(
        description = paste0(
          "Descriptive measures of fold changes ('", sel$estimator, "') against reference axes, gene sets and ",
          "contrast pairs; see the module documentation for each measure. subset: ", params$subset,
          "; set_subset: ", params$set_subset, "."
        ),
        estimator = sel$estimator, subset = params$subset, set_subset = params$set_subset,
        axes = lapply(params$axes, function(a) {
          list(
            id = a$id, resource = a$resource, numerator = unlist(a$numerator), denominator = unlist(a$denominator),
            pseudocount = a$pseudocount, threshold = a$threshold
          )
        }),
        gene_sets = list(resource = params$gene_sets$resource %||% NA_character_, sets = sets),
        sign_test_threshold = thr %||% NA, norm_method = params$norm_method, prior_count = params$prior_count,
        units = list(value = "depends on measure: correlation, log2_ratio, fraction, log2_cpm difference"),
        caveats = c(
          "Reference axes come from other platforms and tissues: correlations with them are descriptive.",
          "The binomial p-value of sign_concordance treats genes as independent (conditional_on_contrast)."
        )
      ),
      axis_genes = list(
        description = paste0(
          "Axis value per gene of subset '", params$subset, "' (NA: gene absent from the reference); class high ",
          "(axis >= threshold), low (axis <= -threshold) or none."
        ),
        units = list(axis = "log2_ratio")
      )
    )
  )
}

reference_axis_row <- function(measure, axis_id = NA_character_, set_id = NA_character_, contrast_id = NA_character_,
                               pair_id = NA_character_, sample_id = NA_character_, value = NA_real_,
                               n = NA_integer_, n_high = NA_integer_, n_low = NA_integer_,
                               n_concordant = NA_integer_, pvalue = NA_real_, inference_scope = NA_character_) {
  data.frame(
    measure = measure, axis_id = axis_id, set_id = set_id, contrast_id = contrast_id, pair_id = pair_id,
    sample_id = sample_id, value = as.numeric(value), n = as.integer(n), n_high = as.integer(n_high),
    n_low = as.integer(n_low), n_concordant = as.integer(n_concordant), pvalue = as.numeric(pvalue),
    inference_scope = inference_scope, stringsAsFactors = FALSE
  )
}

# Axis values (one column per axis) for the given features, from HPA exports.
reference_axis_values <- function(axes, feature_ids, context) {
  out <- matrix(NA_real_, length(feature_ids), length(axes), dimnames = list(NULL, vapply(axes, `[[`, "", "id")))
  for (a in axes) {
    path <- context$resources[[a$resource]]
    if (is.null(path)) abort("reference_axis: resource '", a$resource, "' is not available in the cache")
    kind <- context$cfg$resources$items[[a$resource]]$kind %||% ""
    if (!kind %in% c("hpa_immune_cell", "hpa_single_cell")) {
      abort("reference_axis: resource '", a$resource, "' needs kind hpa_immune_cell or hpa_single_cell")
    }
    num <- unlist(a$numerator)
    den <- unlist(a$denominator)
    long <- read_reference_profiles(path, kind, unique(c(num, den)))
    side <- function(cells) {
      m <- vapply(cells, function(cl) {
        x <- long[long$cell == cl, , drop = FALSE]
        x$value[match(feature_ids, x$feature_id)]
      }, numeric(length(feature_ids)))
      rowMeans(matrix(m, nrow = length(feature_ids)))
    }
    out[, a$id] <- log2((side(num) + a$pseudocount) / (side(den) + a$pseudocount))
  }
  if (length(axes) > 0 && all(is.na(out))) {
    abort("reference_axis: no feature has a reference value (feature IDs must be the Ensembl IDs of the export)")
  }
  out
}

reference_axis_genes <- function(axes, values, genes, in_subset, ct) {
  empty <- data.frame(
    axis_id = character(), feature_id = character(), gene_symbol = character(), axis = numeric(),
    class = character(), stringsAsFactors = FALSE
  )
  i <- match(in_subset, genes$feature_id)
  out <- lapply(axes, function(a) {
    v <- values[i, a$id]
    data.frame(
      axis_id = a$id, feature_id = in_subset, gene_symbol = ct$gene_symbol[match(in_subset, ct$feature_id)],
      axis = v, class = ifelse(is.na(v), NA_character_, ifelse(v >= a$threshold, "high",
        ifelse(v <= -a$threshold, "low", "none")
      )), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(list(empty), out))
}

reference_axis_scores <- function(se, params, axis_genes, context) {
  keep <- config_filter_genes(se, context$cfg, context)
  counts <- SummarizedExperiment::assay(se, "counts")[keep, , drop = FALSE]
  lc <- log_cpm(counts, norm_method = params$norm_method, prior_count = params$prior_count)
  samples <- as.data.frame(SummarizedExperiment::colData(se))
  if (is.null(samples$sample_id)) samples$sample_id <- colnames(se)
  rows <- list()
  for (a in params$axes) {
    g <- axis_genes[axis_genes$axis_id == a$id, , drop = FALSE]
    high <- intersect(g$feature_id[g$class %in% "high"], rownames(lc))
    low <- intersect(g$feature_id[g$class %in% "low"], rownames(lc))
    score <- if (length(high) > 0 && length(low) > 0) {
      colMeans(lc[high, , drop = FALSE]) - colMeans(lc[low, , drop = FALSE])
    } else {
      stats::setNames(rep(NA_real_, ncol(lc)), colnames(lc))
    }
    for (s in colnames(lc)) {
      rows[[length(rows) + 1]] <- reference_axis_row("sample_axis_score",
        axis_id = a$id, sample_id = s, value = score[[s]], n_high = length(high), n_low = length(low)
      )
    }
    for (cn in context$cfg$contrasts %||% list()) {
      num <- contrast_samples(cn$numerator, samples)
      den <- contrast_samples(cn$denominator, samples)
      rows[[length(rows) + 1]] <- reference_axis_row("axis_score_shift",
        axis_id = a$id, contrast_id = cn$id, value = mean(score[num]) - mean(score[den]),
        n_high = length(high), n_low = length(low)
      )
    }
  }
  rows
}
