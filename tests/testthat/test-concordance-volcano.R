test_that("concordance metrics count concordant and discordant genes above a threshold", {
  x <- c(2, -2, 3, 0.5, -3)
  y <- c(1.5, 2, 2, 3, -2)
  m <- concordance_metrics(x, y, 1)
  expect_equal(unname(m[c("n_both", "n_concordant", "n_discordant", "n_conc_up", "n_conc_down")]), c(4, 3, 1, 2, 1))
  expect_equal(unname(m["pct_concordant"]), 75)
  expect_equal(unname(concordance_metrics(1, 1, 5)["pct_concordant"]), 0)
})

test_that("gene subsets combine tokens; strata use quantiles of all genes", {
  g <- data.frame(
    feature_id = letters[1:6], expression = 1:6, masked = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE),
    chromosome = c("1", "X", "1", "Y", "1", "1")
  )
  expect_equal(which(gene_subset(g, "no_masked+no_xy")), c(3, 5))
  expect_equal(which(gene_subset(g, "no_masked+ge_median")), c(4, 5))
  expect_equal(which(gene_subset(g, "in_masked")), c(1, 6))
  expect_error(gene_subset(g, "no_unknownflag"), "unknown token 'no_unknownflag'")
})

test_that("concordance module reports per method, subset, threshold and pair; missing contrasts are errors", {
  ids <- paste0("g", 1:30)
  set.seed(2)
  ct <- rbind(
    new_contrast_table("a", ids, rnorm(30), method = "m", source = "computed"),
    new_contrast_table("b", ids, rnorm(30), method = "m", source = "computed")
  )
  cfg <- list(contrast_pairs = list(list(id = "ab", first = "a", second = "b")))
  params <- list(thresholds = c(0.5, 1), subsets = "all")
  res <- run_concordance(list(contrast_stats = ct), params, list(cfg = cfg, organism = "human"))
  expect_equal(nrow(res$tables$concordance_summary), 2)
  expect_s3_class(res$figures$concordance_scatter, "ggplot")
  expect_registry_outputs(res, "concordance")
  cfg$contrast_pairs[[1]]$second <- "zz"
  expect_error(
    run_concordance(list(contrast_stats = ct), list(thresholds = 1), list(cfg = cfg, organism = "human")),
    "without estimates from a common method"
  )
  cfg$contrast_pairs[[1]]$second <- "b"
  extra <- new_contrast_table("de", ids, rnorm(30), method = "deseq2_wald", source = "computed")
  res <- run_concordance(
    list(contrast_stats = rbind(ct, extra)), list(thresholds = 1), list(cfg = cfg, organism = "human")
  )
  expect_equal(unique(res$tables$concordance_summary$method), "m")
})

test_that("volcano classes use the contrast's own significance metric and linear fold-change threshold", {
  x <- read_de_table(fixture("scenario_c", "comparison.xlsx"), "c", "deseq2", "DESeq2_results", df_residual = 4)
  cls <- classify_de(x, list(padj = 0.05, abs_fold_change = 1.5))
  up <- read_de_table(fixture("scenario_c", "comparison.xlsx"), "c", "deseq2", "padj<=0.05+FoldChange>=1.5", df_residual = 4)
  down <- read_de_table(fixture("scenario_c", "comparison.xlsx"), "c", "deseq2", "padj<=0.05+FoldChange<=-1.5", df_residual = 4)
  expect_equal(sum(cls == "UP"), nrow(up))
  expect_equal(sum(cls == "DOWN"), nrow(down))
  n <- read_de_table(fixture("scenario_c", "noiseq_results.tsv"), "n", "noiseq", df_residual = 0)
  ctx <- list(cfg = list(thresholds = list(posterior_prob = 0.8, abs_fold_change = 1.5)))
  res <- run_volcano(list(contrast_stats = n), list(label_top = 3), ctx)
  expect_equal(unique(res$tables$volcano_points$significance_metric), "posterior_prob")
  expect_registry_outputs(res, "volcano")
  # A provider model with 1 residual df: the volcano points say the call is exploratory.
  one_df <- read_de_table(fixture("scenario_c", "comparison.xlsx"), "c", "deseq2", "DESeq2_results", df_residual = 1)
  pts <- run_volcano(list(contrast_stats = one_df), list(), ctx)$tables$volcano_points
  expect_true(all(pts$exploratory))
  expect_false(any(res$tables$volcano_points$exploratory)) # NOISeq read with df_residual = 0
  expect_equal(res$figures$volcano_n$labels$y, "posterior probability")
  desc <- new_contrast_table("d", "g", 1, method = "m", source = "computed")
  expect_error(
    run_volcano(list(contrast_stats = desc), list(), list(cfg = list())),
    "no contrast with a significance metric"
  )
})

test_that("concordance reports how many samples the two contrasts of a pair share", {
  set.seed(3)
  m <- matrix(rpois(400, 50), 100, dimnames = list(sprintf("g%03d", 1:100), paste0("s", 1:4)))
  se <- toy_se(m)
  cfg <- list(
    contrasts = list(
      list(id = "a", numerator = "s1", denominator = "s2"), list(id = "b", numerator = "s3", denominator = "s4"),
      list(id = "c", numerator = "s1", denominator = "s3")
    ),
    contrast_pairs = list(list(id = "ab", first = "a", second = "b"), list(id = "ac", first = "a", second = "c"))
  )
  ct <- estimate_lfc(se, cfg$contrasts)
  imported <- transform(ct[ct$contrast_id == "b", ], contrast_id = "x", source = "imported:x")
  s <- run_concordance(list(contrast_stats = ct, se = se), list(), list(cfg = cfg, organism = "human"))$tables$concordance_summary
  expect_equal(unique(s$shared_samples[s$pair_id == "ab"]), 0L)
  expect_equal(unique(s$shared_samples[s$pair_id == "ac"]), 1L)
  cfg$contrast_pairs <- list(list(id = "ax", first = "a", second = "x"))
  s <- run_concordance(list(contrast_stats = rbind(ct, imported), se = se), list(), list(cfg = cfg, organism = "human"))
  expect_true(all(is.na(s$tables$concordance_summary$shared_samples)))
})
