#' Read a differential-expression table
#'
#' Imports a provider's DE results without recomputing them. DESeq2 and NOISeq column names are recognised;
#' any other layout uses an explicit column mapping. A NOISeq posterior probability is kept as
#' `posterior_prob` with `significance_metric = "posterior_prob"`: it is never relabelled as an adjusted
#' p-value. A `FoldChange` column is checked against `log2FoldChange` and recorded as signed (`-1/FC` for
#' negative log fold changes) or unsigned.
#'
#' @param file `character(1)`. TSV, CSV or Excel file.
#' @param contrast_id `character(1)`. Identifier of the contrast.
#' @param adapter `character(1)`. `"deseq2"`, `"noiseq"` or `"generic"`.
#' @param sheet `character(1)` or `NULL`. Excel sheet.
#' @param columns Named `list` or `character`. Canonical column -> input column; overrides defaults and is
#'   required for `"generic"` (at least `feature_id` and `lfc`; optionally `lfc_se`, `stat`, `pvalue`,
#'   `padj`, `base_mean`, `posterior_prob`).
#' @param df_residual `integer(1)` or `NULL`. Residual degrees of freedom if known.
#' @return A canonical contrast table (see [new_contrast_table()]) with attribute `checks` (a list
#'   describing the fold-change sign convention and missing values).
#' @examples
#' f <- tempfile(fileext = ".tsv")
#' writeLines(c(
#'   "Gene_ID\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj",
#'   "A\t10\t1\t0.5\t2\t0.04\t0.08"
#' ), f)
#' read_de_table(f, "a_vs_b", "deseq2")
#' @export
read_de_table <- function(file, contrast_id, adapter = c("deseq2", "noiseq", "generic"), sheet = NULL,
                          columns = NULL, df_residual = NULL) {
  adapter <- match.arg(adapter)
  x <- read_any_table(file, sheet)
  columns <- as.list(columns %||% list())
  defaults <- switch(adapter,
    deseq2 = list(
      lfc = "log2FoldChange", lfc_se = "lfcSE", stat = "stat", pvalue = "pvalue", padj = "padj",
      base_mean = "baseMean"
    ),
    noiseq = list(lfc = "M", posterior_prob = "prob"),
    generic = list()
  )
  map <- utils::modifyList(defaults, columns)
  if (is.null(map$feature_id)) {
    candidates <- intersect(c("Gene_ID", "gene_id", "GeneID", "Geneid", "feature_id", "id"), names(x))
    map$feature_id <- if (length(candidates) > 0) candidates[[1]] else names(x)[[1]]
  }
  if (adapter == "generic" && is.null(columns$lfc)) abort("DE table ", file, ": the generic adapter needs columns$lfc")
  assert_columns(x, unlist(map), paste0("DE table ", file, if (!is.null(sheet)) paste0(" [", sheet, "]")))
  get <- function(k) if (is.null(map[[k]])) NA_real_ else as.numeric(x[[map[[k]]]])
  metric <- if (!is.null(map$padj)) "padj" else if (!is.null(map$posterior_prob)) "posterior_prob" else "none"
  extra <- data.frame(row.names = seq_len(nrow(x)))
  if (!is.null(map$posterior_prob)) extra$posterior_prob <- get("posterior_prob")
  checks <- list(
    n_rows = nrow(x), n_missing_lfc = sum(is.na(get("lfc"))),
    n_missing_padj = if (metric == "padj") sum(is.na(get("padj"))) else NA_integer_
  )
  if ("FoldChange" %in% names(x)) {
    lfc <- get("lfc")
    fc <- as.numeric(x$FoldChange)
    signed <- ifelse(lfc >= 0, 2^lfc, -2^(-lfc))
    unsigned <- 2^lfc
    ok <- !is.na(lfc) & !is.na(fc)
    checks$fold_change_convention <- if (isTRUE(all.equal(fc[ok], signed[ok], tolerance = 1e-6))) {
      "signed (-1/FC for lfc < 0)"
    } else if (isTRUE(all.equal(fc[ok], unsigned[ok], tolerance = 1e-6))) {
      "unsigned (2^lfc)"
    } else {
      "inconsistent with log2FoldChange"
    }
    extra$fold_change_provider <- fc
  }
  if (ncol(extra) == 0) extra <- NULL
  ct <- new_contrast_table(
    contrast_id = contrast_id, feature_id = as.character(x[[map$feature_id]]), lfc = get("lfc"),
    lfc_se = get("lfc_se"), stat = get("stat"), pvalue = get("pvalue"), padj = get("padj"),
    base_mean = get("base_mean"), significance_metric = metric, method = paste0("imported_", adapter),
    df_residual = df_residual %||% NA_integer_, source = paste0("imported:", basename(file)),
    inference_scope = if (metric == "none") NA_character_ else "sample", extra = extra
  )
  attr(ct, "checks") <- checks
  attr(ct, "provenance") <- file_provenance(
    file, paste0("de_table_", adapter),
    list(sheet = sheet %||% NA, columns = map)
  )
  ct
}
