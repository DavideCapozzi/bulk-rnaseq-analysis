make_toy <- function() {
  set.seed(3)
  m <- matrix(rpois(60 * 4, 200), 60, dimnames = list(paste0("g", 1:60), c("t1", "t2", "c1", "c2")))
  m[1:5, 1:2] <- m[1:5, 1:2] * 4L
  rd <- S4Vectors::DataFrame(
    feature_id = rownames(m), gene_symbol = c("XIST", "TSIX", paste0("G", 3:60)),
    biotype = c(rep("lncRNA", 2), rep("protein_coding", 57), "rRNA"),
    chromosome = c("X", "X", "Y", "MT", rep("1", 56))
  )
  cd <- S4Vectors::DataFrame(
    sample_id = colnames(m), cond = factor(c("trt", "trt", "ctrl", "ctrl"), c("ctrl", "trt")),
    batch = factor(c("b1", "b2", "b1", "b2")), row.names = colnames(m)
  )
  toy_se(m, rd, cd)
}

test_that("filter_genes applies biotype, mitochondrial, expression and sex-chromosome rules", {
  se <- make_toy()
  k <- filter_genes(
    se,
    min_cpm = 0, min_samples = 1, sex_genes = "drop_xist_y", x_inactivation_genes = c("XIST", "TSIX")
  )
  expect_false(any(k[c(1, 2, 3, 4, 60)]))
  expect_true(all(k[5:59]))
  expect_true(filter_genes(se, min_cpm = 0, min_samples = 1, sex_genes = "keep")[3])
  expect_false(filter_genes(se, min_cpm = 0, min_samples = 1, sex_genes = "drop_xy")[1])
})

test_that("log-CPM difference equals the difference of edgeR log-CPM and has 0 df for one-vs-one", {
  se <- make_toy()
  ct <- estimate_lfc(se, list(list(id = "t1_vs_c1", numerator = "t1", denominator = "c1")), prior_count = 2)
  lc <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(SummarizedExperiment::assay(se))), log = TRUE, prior.count = 2)
  expect_equal(ct$lfc, unname(lc[, "t1"] - lc[, "c1"]))
  expect_equal(unique(ct$df_residual), 0L)
  expect_true(all(is.na(ct$pvalue)))
  grp <- estimate_lfc(se, list(list(id = "trt_vs_ctrl", numerator = "cond=trt", denominator = "cond=ctrl")))
  expect_equal(unique(grp$df_residual), 2L)
  expect_equal(grp$lfc, unname(rowMeans(lc[, 1:2]) - rowMeans(lc[, 3:4])))
})

test_that("reference normalisation uses the reference genes' share of each library", {
  m <- matrix(c(10L, 30L, 60L, 20L, 30L, 50L), 3)
  nf <- reference_norm_factors(m, 1)
  expect_equal(prod(nf), 1)
  expect_equal(nf[2] / nf[1], (20 / 100) / (10 / 100))
  se <- toy_se(matrix(rpois(40, 100), 10, dimnames = list(paste0("g", 1:10), paste0("s", 1:4))))
  ab <- list(list(id = "a", numerator = "s1", denominator = "s2"))
  expect_error(estimate_lfc(se, ab, norm_method = "reference", reference_genes = "NOPE"), "none of the reference genes")
})

test_that("apeglm estimator runs on single-sample contrasts and is flagged exploratory with 1 df", {
  skip_if_not_installed("apeglm")
  se <- make_toy()
  # DESeq2 also warns about NaNs in its outlier cutoff with so few samples; only the guard's warning is checked.
  w <- character()
  ct <- withCallingHandlers(
    estimate_lfc_apeglm(se, list(list(id = "b1", numerator = "t1", denominator = "c1")), design = ~ batch + cond),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("exploratory", w)))
  expect_true(attr(ct, "exploratory"))
  expect_equal(unique(ct$df_residual), 1L)
  expect_true(all(is.finite(ct$lfc)))
  xt <- list(list(id = "x", numerator = "c1", denominator = "t1"))
  expect_error(suppressWarnings(estimate_lfc_apeglm(se, xt, design = ~ batch + cond)), "reference level")
})

test_that("contrast_samples resolves sample ids and factor terms", {
  ss <- data.frame(sample_id = c("a", "b", "c"), g = c("x", "y", "x"))
  expect_equal(contrast_samples(c("g=x"), ss), c("a", "c"))
  expect_error(contrast_samples("zz", ss), "unknown sample_id")
  expect_error(contrast_samples("h=x", ss), "no sample sheet column 'h'")
})

test_that("contrasts module returns descriptive fold changes and the declared tables", {
  set.seed(1)
  m <- matrix(rpois(400, 50), 100, dimnames = list(sprintf("g%03d", 1:100), paste0("s", 1:4)))
  cfg <- list(
    contrasts = list(list(id = "s1_vs_s2", numerator = "s1", denominator = "s2")),
    gene_filter = list(min_cpm = 1, min_samples = 1, sex_genes = "keep")
  )
  res <- run_contrasts(list(se = toy_se(m)), list(), list(cfg = cfg, gene_lists = list(), organism = "human"))
  expect_registry_outputs(res, "contrasts")
  ct <- res$tables$contrast_stats
  expect_equal(unique(ct$method), "logcpm_difference")
  expect_true(all(is.na(ct$pvalue)))
  expect_equal(unique(ct$df_residual), 0L)
  expect_equal(res$tables$contrast_spread$subset, "all")
})
