#' Select genes with a subset expression
#'
#' Subsets are written as tokens joined by `+`, all of which must hold:
#' `all`; `no_<flag>` (exclude genes where logical column `<flag>` is `TRUE`); `in_<flag>` (keep genes where it is
#' `TRUE`); `no_xy` (exclude X and Y chromosomes); `ge_median`, `ge_q25`, `ge_q75` (expression at least the
#' median / quartile of **all** genes in `genes`, not of the subset, so strata are comparable across subsets).
#'
#' @param genes `data.frame` with one row per gene: `expression` (numeric) and any logical flag columns;
#'   `chromosome` is needed for `no_xy`.
#' @param subset `character(1)`. Subset expression, e.g. `"no_masked+ge_median"`.
#' @param organism `character(1)`. Organism key (sex chromosome names).
#' @return `logical`, one value per row of `genes`.
#' @examples
#' g <- data.frame(expression = 1:4, masked = c(TRUE, FALSE, FALSE, FALSE))
#' gene_subset(g, "no_masked+ge_median")
#' @export
gene_subset <- function(genes, subset, organism = "human") {
  keep <- rep(TRUE, nrow(genes))
  org <- organism_info(organism)
  for (tok in strsplit(subset, "+", fixed = TRUE)[[1]]) {
    keep <- keep & switch(tok,
      all = TRUE,
      no_xy = !genes$chromosome %in% c(org$chromosome_x, org$chromosome_y),
      ge_median = genes$expression >= stats::median(genes$expression),
      ge_q25 = genes$expression >= stats::quantile(genes$expression, 0.25),
      ge_q75 = genes$expression >= stats::quantile(genes$expression, 0.75),
      {
        flag <- sub("^(no|in)_", "", tok)
        if (!grepl("^(no|in)_", tok) || !flag %in% names(genes)) {
          abort(
            "gene subset '", subset, "': unknown token '", tok, "' (flags available: ",
            paste(setdiff(names(genes), c("feature_id", "expression")), collapse = ", "), ")"
          )
        }
        v <- as.logical(genes[[flag]])
        v[is.na(v)] <- FALSE
        if (startsWith(tok, "no_")) !v else v
      }
    )
  }
  keep
}
