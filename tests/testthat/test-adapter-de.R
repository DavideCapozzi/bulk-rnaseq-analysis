test_that("DESeq2 workbook: every row read despite <dimension ref=\"A1\"/>, sheet names with '>' work", {
  f <- fixture("scenario_c", "comparison.xlsx")
  stat <- as.data.frame(readxl::read_excel(f, sheet = "Statistic"))
  all <- read_de_table(f, "c", "deseq2", "DESeq2_results", df_residual = 4)
  up <- read_de_table(f, "c", "deseq2", "padj<=0.05+FoldChange>=1.5", df_residual = 4)
  expect_gt(nrow(all), 100)
  expect_equal(nrow(up), stat$UP)
  expect_equal(attr(all, "checks")$fold_change_convention, "signed (-1/FC for lfc < 0)")
  expect_equal(unique(all$significance_metric), "padj")
  expect_equal(unique(all$inference_scope), "sample")
  expect_equal(unique(all$source), "imported:comparison.xlsx")
  expect_true(all(up$fold_change_provider >= 1.5))
})

test_that("NOISeq posterior probabilities are not relabelled as adjusted p-values", {
  x <- read_de_table(fixture("scenario_c", "noiseq_results.tsv"), "n", "noiseq", df_residual = 0)
  expect_equal(unique(x$significance_metric), "posterior_prob")
  expect_true(all(is.na(x$padj)))
  raw <- utils::read.delim(fixture("scenario_c", "noiseq_results.tsv"))
  expect_equal(x$posterior_prob, raw$prob)
  expect_equal(x$lfc, raw$M)
  # Without replicates (df 0) the posterior probabilities are conditional on the contrast, not sample-level.
  expect_equal(unique(x$inference_scope), "conditional_on_contrast")
  expect_error(read_de_table(fixture("scenario_c", "noiseq_results.tsv"), "n", "noiseq"), "need the residual degrees")
})

test_that("generic DE tables need an explicit column mapping", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("gene,effect,q", "A,1.5,0.01", "B,-2,0.2"), f)
  expect_error(read_de_table(f, "g", "generic"), "needs columns\\$lfc")
  x <- read_de_table(f, "g", "generic",
    columns = list(feature_id = "gene", lfc = "effect", padj = "q"), df_residual = 2
  )
  expect_equal(x$padj, c(0.01, 0.2))
  expect_error(read_de_table(f, "g", "generic", columns = list(lfc = "nope")), "missing column")
})
