test_that("contrast tables are built and validated", {
  ct <- new_contrast_table("c", c("g1", "g2"), c(1, -1), method = "logcpm_difference", source = "computed")
  expect_equal(names(ct), contrast_table_columns())
  expect_error(validate_contrast_table(ct[, -1]), "missing column\\(s\\) 'contrast_id'")
  expect_error(validate_contrast_table(ct[, rev(names(ct))]), "fixed columns must come first")
  bad <- ct
  bad$significance_metric <- "qvalue"
  expect_error(validate_contrast_table(bad), "invalid significance_metric")
})

test_that("p-values need an inference scope and never come from a 0-df sample-level model", {
  expect_error(
    new_contrast_table("c", "g", 1, pvalue = 0.01, method = "m", source = "computed"),
    "need an inference_scope"
  )
  expect_error(
    new_contrast_table("c", "g", 1,
      pvalue = 0.01, method = "m", source = "computed", df_residual = 0L, inference_scope = "sample"
    ),
    "0 residual degrees of freedom"
  )
  ok <- new_contrast_table("c", "g", 1,
    pvalue = 0.01, method = "camera", source = "computed", df_residual = 0L,
    inference_scope = "conditional_on_contrast"
  )
  expect_equal(ok$inference_scope, "conditional_on_contrast")
})

test_that("posterior probabilities keep their name", {
  expect_error(
    new_contrast_table("c", "g", 1, significance_metric = "posterior_prob", method = "m", source = "x"),
    "needs an extension column 'posterior_prob'"
  )
  ct <- new_contrast_table("c", "g", 1,
    significance_metric = "posterior_prob", method = "m", source = "x",
    extra = data.frame(posterior_prob = 0.9)
  )
  sig <- significance_values(ct)
  expect_equal(sig$metric, "posterior_prob")
  expect_equal(sig$values, 0.9)
  expect_true(is.na(ct$padj))
})

test_that("SummarizedExperiment and QC validators reject malformed objects", {
  m <- matrix(c(1.5, 2, 3, 4), 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  expect_error(validate_se(toy_se(m)), "must contain integers")
  m2 <- matrix(1:4, 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  expect_silent(validate_se(toy_se(m2)))
  expect_error(validate_qc_table(data.frame(sample_id = "a", metric = "m", value = "x", source = "s")), "numeric")
})
