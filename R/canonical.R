#' Fixed columns of the canonical contrast table
#'
#' @return `character` vector of column names, in order.
#' @examples
#' contrast_table_columns()
#' @export
contrast_table_columns <- function() {
  c(
    "contrast_id", "feature_id", "gene_symbol", "lfc", "lfc_se", "stat", "pvalue", "padj", "base_mean",
    "significance_metric", "method", "df_residual", "source", "inference_scope"
  )
}

#' Build a canonical contrast table
#'
#' One row per contrast and feature. Columns that a producer cannot fill stay `NA`; extension columns
#' (e.g. `posterior_prob`, `ave_log2cpm`, `fold_change_signed`) go after the fixed ones.
#'
#' @param contrast_id `character(1)`. Contrast identifier.
#' @param feature_id `character`. Input gene identifiers.
#' @param lfc `numeric`. log2 fold change, numerator over denominator.
#' @param gene_symbol,lfc_se,stat,pvalue,padj,base_mean Optional per-feature vectors (`NA` when absent).
#'   `base_mean` is the mean of normalised counts.
#' @param significance_metric `character(1)`. `"padj"`, `"posterior_prob"` or `"none"`.
#' @param method `character(1)`. Estimator or test, e.g. `"logcpm_difference"`, `"imported_deseq2"`.
#' @param df_residual `integer(1)`. Residual degrees of freedom of the producing model, `NA` if unknown.
#' @param source `character(1)`. `"computed"` or `"imported:<file>"`.
#' @param inference_scope `character(1)`. `"sample"` when p-values come from a test whose units are samples,
#'   `"conditional_on_contrast"` for gene-level tests on a given contrast, `NA` when there is no p-value.
#' @param extra Optional `data.frame` of extension columns (same number of rows).
#' @return A `data.frame` that passes [validate_contrast_table()].
#' @examples
#' ct <- new_contrast_table("a_vs_b", c("g1", "g2"), c(1.2, -0.3),
#'   significance_metric = "none",
#'   method = "logcpm_difference", df_residual = 0L, source = "computed"
#' )
#' ct
#' @export
new_contrast_table <- function(contrast_id, feature_id, lfc, gene_symbol = NA_character_, lfc_se = NA_real_,
                               stat = NA_real_, pvalue = NA_real_, padj = NA_real_, base_mean = NA_real_,
                               significance_metric = "none", method, df_residual = NA_integer_, source,
                               inference_scope = NA_character_, extra = NULL) {
  n <- length(feature_id)
  ct <- data.frame(
    contrast_id = rep(contrast_id, n), feature_id = as.character(feature_id),
    gene_symbol = rep_len(as.character(gene_symbol), n), lfc = as.numeric(lfc),
    lfc_se = rep_len(as.numeric(lfc_se), n), stat = rep_len(as.numeric(stat), n),
    pvalue = rep_len(as.numeric(pvalue), n), padj = rep_len(as.numeric(padj), n),
    base_mean = rep_len(as.numeric(base_mean), n), significance_metric = rep(significance_metric, n),
    method = rep(method, n), df_residual = rep(as.integer(df_residual), n), source = rep(source, n),
    inference_scope = rep(as.character(inference_scope), n), stringsAsFactors = FALSE
  )
  if (!is.null(extra)) ct <- cbind(ct, extra)
  validate_contrast_table(ct)
}

#' Validate a canonical contrast table
#'
#' @param ct `data.frame`. Candidate contrast table (one or more contrasts).
#' @return `ct`, invisibly unchanged, or an error naming the column at fault.
#' @examples
#' ct <- new_contrast_table("c1", "g1", 1, method = "m", source = "computed")
#' validate_contrast_table(ct)
#' @export
validate_contrast_table <- function(ct) {
  what <- "contrast table"
  if (!is.data.frame(ct)) abort(what, ": expected a data.frame")
  fixed <- contrast_table_columns()
  assert_columns(ct, fixed, what)
  if (!identical(names(ct)[seq_along(fixed)], fixed)) {
    abort(what, ": fixed columns must come first, in the order ", paste(fixed, collapse = ", "))
  }
  for (col in c("lfc", "lfc_se", "stat", "pvalue", "padj", "base_mean")) {
    if (!is.numeric(ct[[col]])) abort(what, ": column '", col, "' must be numeric")
  }
  bad_metric <- setdiff(unique(ct$significance_metric), c("padj", "posterior_prob", "none"))
  if (length(bad_metric) > 0) abort(what, ": invalid significance_metric ", paste(bad_metric, collapse = ", "))
  if (any(ct$significance_metric == "posterior_prob") && !"posterior_prob" %in% names(ct)) {
    abort(what, ": significance_metric 'posterior_prob' needs an extension column 'posterior_prob'")
  }
  scope <- unique(ct$inference_scope[!is.na(ct$inference_scope)])
  bad_scope <- setdiff(scope, c("sample", "conditional_on_contrast"))
  if (length(bad_scope) > 0) abort(what, ": invalid inference_scope ", paste(bad_scope, collapse = ", "))
  has_p <- !is.na(ct$pvalue) | !is.na(ct$padj)
  if (any(has_p & is.na(ct$inference_scope))) {
    abort(what, ": rows with p-values need an inference_scope ('sample' or 'conditional_on_contrast')")
  }
  computed_sample_p <- has_p & ct$inference_scope %in% "sample" & ct$source == "computed"
  if (any(computed_sample_p & ct$df_residual %in% 0L)) {
    abort(
      what, ": sample-level p-values with 0 residual degrees of freedom (contrast ",
      paste(unique(ct$contrast_id[computed_sample_p & ct$df_residual %in% 0L]), collapse = ", "), ")"
    )
  }
  key <- paste(ct$contrast_id, ct$method, ct$feature_id, sep = "\r")
  if (anyDuplicated(key)) abort(what, ": duplicated contrast_id/method/feature_id rows")
  invisible(ct)
}

#' Significance values of a contrast table
#'
#' Returns adjusted p-values or posterior probabilities according to `significance_metric`, so that callers
#' never need to relabel one as the other.
#'
#' @param ct Canonical contrast table with a single significance metric.
#' @return A list with `metric` (`character(1)`) and `values` (`numeric`, `NA` when the metric is `"none"`).
#' @examples
#' ct <- new_contrast_table("c", c("g1", "g2"), c(1, 2),
#'   padj = c(0.01, 0.2), significance_metric = "padj",
#'   method = "m", source = "imported:x", inference_scope = "sample"
#' )
#' significance_values(ct)
#' @export
significance_values <- function(ct) {
  metric <- unique(ct$significance_metric)
  if (length(metric) != 1) {
    abort("significance_values(): one significance metric per call (got ", paste(metric, collapse = ", "), ")")
  }
  values <- switch(metric,
    padj = ct$padj,
    posterior_prob = ct$posterior_prob,
    none = rep(NA_real_, nrow(ct))
  )
  list(metric = metric, values = values)
}

#' Validate a sample-level SummarizedExperiment
#'
#' @param se A `SummarizedExperiment`.
#' @return `se` invisibly, or an error naming the assay or column at fault.
#' @examples
#' m <- matrix(1:4, 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = m),
#'   rowData = S4Vectors::DataFrame(feature_id = c("g1", "g2")),
#'   colData = S4Vectors::DataFrame(sample_id = c("s1", "s2"), row.names = c("s1", "s2"))
#' )
#' validate_se(se)
#' @export
validate_se <- function(se) {
  what <- "SummarizedExperiment"
  if (!methods::is(se, "SummarizedExperiment")) abort(what, ": expected a SummarizedExperiment")
  allowed <- c("counts", "tpm", "fpkm", "normalized")
  an <- SummarizedExperiment::assayNames(se)
  bad <- setdiff(an, allowed)
  if (length(bad) > 0) {
    abort(what, ": unknown assay(s) ", paste(bad, collapse = ", "), "; allowed: ", paste(allowed, collapse = ", "))
  }
  if ("counts" %in% an) {
    m <- SummarizedExperiment::assay(se, "counts")
    if (!is.integer(m) && !is_whole_number(m)) abort(what, ": assay 'counts' must contain integers")
    if (any(m < 0, na.rm = TRUE)) abort(what, ": assay 'counts' has negative values")
  }
  rd <- SummarizedExperiment::rowData(se)
  if (!"feature_id" %in% names(rd)) abort(what, ": rowData needs 'feature_id'")
  if (anyDuplicated(rd$feature_id)) abort(what, ": duplicated feature_id in rowData")
  if (!identical(as.character(rd$feature_id), rownames(se))) abort(what, ": rownames must equal rowData$feature_id")
  cd <- SummarizedExperiment::colData(se)
  if (!"sample_id" %in% names(cd) || !identical(as.character(cd$sample_id), colnames(se))) {
    abort(what, ": colnames must equal colData$sample_id")
  }
  invisible(se)
}

#' Validate a QC long table
#'
#' @param qc `data.frame` with `sample_id`, `metric`, `value` (numeric) and `source`.
#' @return `qc` invisibly.
#' @examples
#' qc <- data.frame(sample_id = "s1", metric = "n_unmapped", value = 3, source = "star_counts")
#' validate_qc_table(qc)
#' @export
validate_qc_table <- function(qc) {
  assert_columns(qc, c("sample_id", "metric", "value", "source"), "QC table")
  if (!is.numeric(qc$value)) abort("QC table: column 'value' must be numeric")
  if (anyDuplicated(paste(qc$sample_id, qc$metric, qc$source, sep = "\r"))) {
    abort("QC table: duplicated sample_id/metric/source rows")
  }
  invisible(qc)
}

#' Bind contrast tables with different extension columns
#'
#' @param tables `list` of canonical contrast tables (`NULL` entries are dropped).
#' @return One validated contrast table; extension columns missing from a table are filled with `NA`.
#' @examples
#' a <- new_contrast_table("a", "g", 1,
#'   method = "m", source = "computed", extra = data.frame(ave_log2cpm = 2)
#' )
#' b <- new_contrast_table("b", "g", 2, method = "m", source = "computed")
#' bind_contrast_tables(list(a, b))
#' @export
bind_contrast_tables <- function(tables) {
  tables <- Filter(Negate(is.null), tables)
  if (length(tables) == 0) {
    return(NULL)
  }
  ct <- as.data.frame(data.table::rbindlist(unname(tables), fill = TRUE, use.names = TRUE))
  fixed <- contrast_table_columns()
  validate_contrast_table(ct[, c(fixed, setdiff(names(ct), fixed)), drop = FALSE])
}
