#' Residual degrees of freedom of a design
#'
#' Computed as the number of samples minus the rank of the design matrix, so that confounded or saturated
#' designs are detected whatever the number of groups.
#'
#' @param design A model matrix, or a formula evaluated on `data`.
#' @param data `data.frame` of sample metadata, used when `design` is a formula.
#' @return `integer(1)`, residual degrees of freedom.
#' @examples
#' d <- data.frame(g = factor(c("a", "a", "b", "b")), b = factor(c("x", "y", "x", "y")))
#' residual_df(~ g + b, d) # 4 - 3 = 1
#' @export
residual_df <- function(design, data = NULL) {
  if (inherits(design, "formula")) {
    if (is.null(data)) abort("residual_df(): a formula design needs `data`")
    design <- stats::model.matrix(design, data = data)
  }
  as.integer(nrow(design) - qr(design)$rank)
}

#' Guard sample-level tests against missing residual degrees of freedom
#'
#' Must be called before any test whose units are samples (DESeq2/edgeR/limma tests, dispersion-based
#' p-values). With 0 residual degrees of freedom there is no estimate of between-sample variance, so only
#' descriptive estimators are allowed; with 1 the test runs but is flagged exploratory.
#'
#' @param design A model matrix or a formula (see [residual_df()]).
#' @param data `data.frame` of sample metadata when `design` is a formula.
#' @param what `character(1)`. Name of the analysis, used in messages.
#' @return A list with `df` (`integer`) and `exploratory` (`logical`), invisibly. Errors when `df == 0`.
#' @examples
#' d <- data.frame(g = factor(c("a", "a", "b", "b")))
#' assert_df(~g, d, "example test")$df
#' try(assert_df(~g, d[c(1, 3), , drop = FALSE], "one vs one"))
#' @export
assert_df <- function(design, data = NULL, what = "test") {
  df <- residual_df(design, data)
  if (df <= 0) {
    abort(
      what, ": 0 residual degrees of freedom (n samples - rank of the design matrix); ",
      "sample-level p-values are not allowed, only descriptive estimators"
    )
  }
  exploratory <- df == 1
  if (exploratory) {
    warning(what, ": 1 residual degree of freedom; results are flagged exploratory", call. = FALSE)
  }
  invisible(list(df = df, exploratory = exploratory))
}
