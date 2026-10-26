#' Describe how a provider normalised counts
#'
#' Compares provider-normalised values with the raw counts of the same samples: whether each sample's values
#' are the counts divided by one constant (a size factor), whether DESeq2 median-of-ratios size factors
#' computed on all genes or on the provider's genes reproduce it, and which simple pre-filter
#' ("at least k reads in at least n samples") yields exactly the provider's gene set.
#'
#' @param counts Integer matrix, genes x samples.
#' @param normalized Numeric matrix with the same columns; rows are the provider's genes (a subset of
#'   `counts`), `NA` rows are ignored.
#' @param min_counts `numeric`. Values of k tried for the pre-filter rule.
#' @return `list(table, facts)`: `table` has one row per sample (`implied_size_factor`,
#'   `implied_ratio_rel_iqr`, `sf_deseq2_all_genes`, `sf_deseq2_provider_genes`); `facts` is a named list for
#'   the ingest report.
#' @examples
#' set.seed(1)
#' m <- matrix(rpois(600, 30), 100, dimnames = list(paste0("g", 1:100), paste0("s", 1:6)))
#' sf <- DESeq2::estimateSizeFactorsForMatrix(m)
#' check_normalization(m, t(t(m) / sf))$facts
#' @export
check_normalization <- function(counts, normalized, min_counts = c(1, 2, 3, 5, 10, 15, 20)) {
  normalized <- normalized[stats::complete.cases(normalized), , drop = FALSE]
  genes <- rownames(normalized)
  if (!all(genes %in% rownames(counts))) abort("check_normalization(): normalised genes absent from the counts")
  raw <- counts[genes, colnames(normalized), drop = FALSE]
  ratio <- raw / normalized
  ratio[!is.finite(ratio) | normalized == 0] <- NA
  implied <- apply(ratio, 2, stats::median, na.rm = TRUE)
  rel_iqr <- apply(ratio, 2, function(r) stats::IQR(r, na.rm = TRUE)) / implied
  sf_all <- DESeq2::estimateSizeFactorsForMatrix(counts[, colnames(normalized), drop = FALSE])
  sf_prov <- DESeq2::estimateSizeFactorsForMatrix(raw)
  table <- data.frame(
    sample_id = colnames(normalized), implied_size_factor = unname(implied), implied_ratio_rel_iqr = unname(rel_iqr),
    sf_deseq2_all_genes = unname(sf_all), sf_deseq2_provider_genes = unname(sf_prov), stringsAsFactors = FALSE
  )
  in_set <- rownames(counts) %in% genes
  rules <- list()
  for (k in min_counts) {
    for (n in seq_len(ncol(counts))) {
      pass <- rowSums(counts >= k) >= n
      rules[[length(rules) + 1]] <- data.frame(min_count = k, min_samples = n, agreement = mean(pass == in_set))
    }
  }
  rules <- do.call(rbind, rules)
  exact <- rules[rules$agreement == 1, , drop = FALSE]
  best <- rules[which.max(rules$agreement), , drop = FALSE]
  max_dev <- function(sf) max(abs(sf - implied) / implied)
  facts <- list(
    normalization_genes = sprintf("%d of %d count genes", length(genes), nrow(counts)),
    normalization_constant_per_sample = sprintf("max relative IQR of counts/normalised = %.2g", max(rel_iqr)),
    normalization_size_factors = sprintf(
      "implied vs DESeq2 median-of-ratios: max relative difference %.2g on provider genes, %.2g on all genes",
      max_dev(sf_prov), max_dev(sf_all)
    ),
    normalization_filter_rule = if (nrow(exact) > 0) {
      paste0("exact: >= ", exact$min_count, " reads in >= ", exact$min_samples, " samples", collapse = "; ")
    } else {
      sprintf(
        "no exact rule among k in {%s}; best: >= %g reads in >= %d samples (agreement %.4f)",
        paste(min_counts, collapse = ", "), best$min_count, best$min_samples, best$agreement
      )
    }
  )
  list(table = table, facts = facts)
}
