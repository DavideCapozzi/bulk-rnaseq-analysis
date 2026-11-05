de_heatmap_fixture <- function() {
  se <- nb_se(n_genes = 200, groups = c("a", "a", "a", "b", "b", "b"), seed = 5, effect = 6)
  SummarizedExperiment::colData(se)$group <- c("a", "a", "a", "b", "b", "b")
  ids <- rownames(se)
  padj <- c(seq(1e-6, 1e-3, length.out = 20), rep(0.5, 180))
  lfc <- c(rep(c(2, -2), 10), rep(0.1, 180))
  ct <- new_contrast_table("b_vs_a", ids, lfc,
    gene_symbol = paste0("S", seq_along(ids)), padj = padj, significance_metric = "padj", method = "deseq2",
    source = "computed", inference_scope = "sample"
  )
  cfg <- list(samples = list(factors = list(group = c("a", "b"))), thresholds = list(padj = 0.05, abs_fold_change = 1.5))
  list(se = se, ct = ct, ctx = list(cfg = cfg))
}

test_that("de_heatmap selects the most significant UP/DOWN genes and z-scores their normalised counts", {
  x <- de_heatmap_fixture()
  res <- run_de_heatmap(list(se = x$se, contrast_stats = x$ct), list(n_top = 8), x$ctx)
  expect_registry_outputs(res, "de_heatmap")
  g <- res$tables$de_heatmap_genes
  expect_equal(g$feature_id, rownames(x$se)[1:8])
  expect_equal(g$rank, 1:8)
  expect_setequal(g$row_order, 1:8)
  expect_equal(unique(g$class), c("UP", "DOWN"))
  v <- res$tables$de_heatmap_values
  expect_equal(nrow(v), 8 * 6)
  counts <- SummarizedExperiment::assay(x$se, "counts")
  sf <- DESeq2::estimateSizeFactorsForMatrix(counts)
  g1 <- v[v$feature_id == rownames(x$se)[1], ]
  expect_equal(g1$log2_value, unname(log2(counts[1, ] / sf + 1)))
  expect_equal(g1$z, as.vector(scale(g1$log2_value)))
  expect_equal(res$meta$de_heatmap_genes$values[[1]], "counts / median-of-ratios size factors")
  expect_named(res$figures, "de_heatmap_b_vs_a")
})

test_that("de_heatmap can take half of the genes from each direction", {
  x <- de_heatmap_fixture()
  ct <- x$ct
  ct$lfc[1:20] <- c(rep(2, 14), rep(-2, 6))
  one <- run_de_heatmap(list(se = x$se, contrast_stats = ct), list(n_top = 8), x$ctx)$tables$de_heatmap_genes
  expect_equal(unique(one$class), "UP")
  both <- run_de_heatmap(list(se = x$se, contrast_stats = ct), list(n_top = 8, per_direction = TRUE), x$ctx)
  g <- both$tables$de_heatmap_genes
  expect_equal(as.vector(table(g$class)[c("UP", "DOWN")]), c(4L, 4L))
  expect_equal(g$feature_id, rownames(x$se)[c(1:4, 15:18)])
})

test_that("de_heatmap prefers the normalized input and skips contrasts without significant genes", {
  x <- de_heatmap_fixture()
  SummarizedExperiment::assay(x$se, "normalized") <- SummarizedExperiment::assay(x$se, "counts") * 1
  SummarizedExperiment::assay(x$se, "normalized")[200, ] <- NA
  res <- run_de_heatmap(list(se = x$se, contrast_stats = x$ct), list(n_top = 3), x$ctx)
  v <- res$tables$de_heatmap_values
  expect_equal(v$log2_value[1], log2(SummarizedExperiment::assay(x$se, "counts")[1, 1] + 1))
  expect_equal(res$meta$de_heatmap_genes$values[[1]], "normalized input")
  ns <- x$ct
  ns$padj <- 0.9
  res <- run_de_heatmap(list(se = x$se, contrast_stats = ns), list(), x$ctx)
  expect_equal(nrow(res$tables$de_heatmap_genes), 0)
  expect_length(res$figures, 0)
})

test_that("de_heatmap fails loudly on bad configuration", {
  x <- de_heatmap_fixture()
  run <- function(p, ct = x$ct) run_de_heatmap(list(se = x$se, contrast_stats = ct), p, x$ctx)
  expect_error(run(list(contrasts = list("zz"))), "without a significance metric: zz")
  expect_error(run(list(values = "normalized")), "needs the normalized input")
  ct <- x$ct
  ct$significance_metric <- "none"
  expect_error(run(list(), ct), "no contrast with a significance metric")
})
