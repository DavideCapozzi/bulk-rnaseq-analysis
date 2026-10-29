test_that("expression QC counts library shares on all genes; ribosomal proteins must be protein-coding", {
  m <- matrix(c(50L, 30L, 10L, 10L, 0L, 0L, 70L, 20L, 5L, 5L), 5,
    dimnames = list(paste0("g", 1:5), c("s1", "s2"))
  )
  rd <- S4Vectors::DataFrame(
    feature_id = rownames(m), gene_symbol = c("MT-CO1", "RPL3", "RPS4X", "HBB", "A"),
    biotype = c("protein_coding", "protein_coding", "lncRNA", "protein_coding", "rRNA"), chromosome = c("MT", "1", "X", "1", "1")
  )
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = m, tpm = m * 1e4), rowData = rd)
  qc <- expression_qc(se, hemoglobin_genes = "HBB", top_n_share = 2)
  v <- function(metric, s) qc$value[qc$metric == metric & qc$sample_id == s]
  expect_equal(v("library_size", "s1"), 100)
  expect_equal(v("pct_counts_mito", "s1"), 50)
  expect_equal(v("pct_counts_ribo_protein", "s1"), 30) # RPS4X is lncRNA here, so it does not count
  expect_equal(v("pct_counts_hemoglobin", "s2"), 5)
  expect_equal(v("pct_counts_top2_genes", "s2"), 90)
  expect_equal(v("pct_tpm_mito", "s1"), 50 * 1e4 / 1e4)
  expect_equal(v("genes_count_ge10", "s1"), 4)
})

test_that("qc_wide keeps metric order and disambiguates metrics reported by two sources", {
  qc <- data.frame(
    sample_id = c("s1", "s2", "s1", "s2", "s1"), metric = c("b", "b", "a", "a", "b"), value = 1:5,
    source = c("x", "x", "x", "x", "y")
  )
  w <- qc_wide(qc, c("s2", "s1"))
  expect_equal(names(w), c("sample_id", "b_x", "a", "b_y"))
  expect_equal(w$a, c(4, 3))
  expect_true(is.na(w$b_y[[1]]))
})

test_that("qc module runs on imported QC alone and fails when nothing is available", {
  qc <- data.frame(sample_id = c("s1", "s2"), metric = "n_unmapped", value = c(3, 4), source = "star_counts")
  res <- run_qc(list(qc = qc), list(), list(cfg = list(), gene_lists = list(), organism = "human"))
  expect_equal(names(res$tables), "qc_summary")
  expect_registry_outputs(res, "qc")
  expect_length(res$figures, 0)
  expect_error(run_qc(list(), list(), list(cfg = list())), "no QC metrics available")
})

test_that("pc_separation labels factors and the interaction; ties go to the first factor", {
  ss <- data.frame(f = c("a", "a", "b", "b"), r = c("x", "y", "x", "y"))
  expect_equal(pc_separation(c(-2, -2, 2, 2), ss, c("f", "r")), "f")
  expect_equal(pc_separation(c(-2, 2, -2, 2), ss, c("f", "r")), "r")
  expect_equal(pc_separation(c(1, -1, -1, 1), ss, c("f", "r")), "interaction")
  expect_equal(pc_separation(c(1, 1, 1, 1), ss, c("f", "r")), "f")
  expect_true(is.na(pc_separation(1:4, ss, character())))
})

test_that("pca_top_variance keeps the most variable genes and reports percent variance", {
  set.seed(3)
  x <- matrix(rnorm(200), 50)
  x[1:5, ] <- x[1:5, ] * 20
  p <- pca_top_variance(x, 5)
  expect_equal(p$n_genes, 5)
  expect_setequal(as.integer(sub("V", "", rownames(p$rotation) %||% as.character(1:5))), 1:5)
  expect_equal(sum(p$ve), 100)
  expect_equal(pca_top_variance(x, Inf)$n_genes, 50)
})

test_that("pca module runs the variant grid, orients PC1 and reproduces seeded thinning", {
  se <- nb_se()
  se$label <- c("A1", "A2", "B1", "B2")
  cfg <- list(
    samples = list(factors = list(g = c("a", "b")), label_column = "label"),
    gene_filter = list(min_cpm = 1, min_samples = 2, sex_genes = "keep")
  )
  params <- list(
    ntop = list(100, "all"), scores_ntop = 100, sex_gene_variants = c("keep", "drop_xy"), symbol_exclusion_variants = "^g00[0-4]",
    thinning = TRUE, orient_pc1 = list(factor = "g", level = "b"), top_loadings = 5
  )
  ctx <- list(cfg = cfg, gene_lists = list(), organism = "human", seed = 7)
  res <- run_pca(list(se = se), params, ctx)
  s <- res$tables$pca_summary
  expect_equal(nrow(s), 2 * 2 * 2 * 2 * 2) # sex x exclusion x data x transform x ntop
  expect_true(all(s$pc1_separates == "g"))
  expect_true(all(grepl("B[12] < B[12]$", s$pc1_order)))
  expect_lt(max(s$n_genes[s$symbol_exclusion != "none"]), min(s$n_genes[s$symbol_exclusion == "none"]))
  expect_equal(res$tables$pc1_top_loadings$direction[[1]], "towards_b")
  expect_equal(sort(unique(res$tables$sample_correlations$method)), c("pearson", "spearman"))
  expect_equal(nrow(res$tables$sample_distances), 16)
  expect_s3_class(res$figures$pca, "ggplot")
  expect_registry_outputs(res, "pca")
  again <- run_pca(list(se = se), params, ctx)
  expect_identical(again$tables$pca_summary, s)
  # The schema default (an empty list, as in a resolved config) also means every factor of samples.factors.
  resolved <- run_pca(list(se = se), c(params, list(separation_factors = list())), ctx)
  expect_identical(resolved$tables$pca_summary, s)
  other_seed <- run_pca(list(se = se), params, modifyList(ctx, list(seed = 8)))
  thinned <- s$data == "thinned"
  expect_false(identical(other_seed$tables$pca_summary$pc1_var[thinned], s$pc1_var[thinned]))
})

test_that("deseq2_wald guards residual degrees of freedom and labels sample-level inference", {
  se <- nb_se(groups = rep(c("a", "b"), each = 3))
  d <- list(id = "b_vs_a", formula = "~ g", contrast = c("g", "b", "a"))
  ct <- deseq2_wald(se, d, min_count = 5, min_samples = 3)
  expect_true(all(ct$inference_scope == "sample"))
  expect_equal(unique(ct$df_residual), 4L)
  expect_false(any(ct$exploratory))
  expect_false(attr(ct, "exploratory"))
  expect_gt(stats::median(ct$lfc[ct$feature_id <= "g0100"]), 1)
  one_each <- se[, c(1, 4)]
  expect_error(deseq2_wald(one_each, d), "0 residual degrees of freedom")
  expect_error(deseq2_wald(se, modifyList(d, list(contrast = c("g", "z", "a")))), "unknown level")
})

test_that("de_test compares computed and imported contrasts and rejects unknown engines", {
  se <- nb_se(groups = rep(c("a", "b"), each = 3))
  cfg <- list(designs = list(list(id = "b_vs_a", formula = "~ g", contrast = c("g", "b", "a"))))
  base <- run_de_test(list(se = se), list(), list(cfg = cfg))
  imported <- base$tables$contrast_stats
  imported$contrast_id <- "provider"
  imported$source <- "imported:x"
  imported$method <- "imported_deseq2"
  params <- list(compare_with = list(list(design = "b_vs_a", contrast_table = "provider")))
  res <- run_de_test(list(se = se, contrast_stats = imported), params, list(cfg = cfg))
  pc <- res$tables$provider_comparison
  expect_equal(pc$n_same_class, nrow(imported))
  expect_equal(pc$lfc_max_abs_diff, 0)
  expect_equal(res$tables$design_summary$df_residual, 4L)
  expect_registry_outputs(res, "de_test")
  params$compare_with[[1]]$contrast_table <- "absent"
  expect_error(run_de_test(list(se = se, contrast_stats = imported), params, list(cfg = cfg)), "'absent' not found")
  cfg$designs[[1]]$engine <- "edger"
  expect_error(run_de_test(list(se = se), list(), list(cfg = cfg)), "not implemented")
})

test_that("check_normalization recovers size factors and an exact pre-filter rule", {
  se <- nb_se(groups = rep(c("a", "b"), each = 3))
  m <- SummarizedExperiment::assay(se, "counts")
  m[1:50, ] <- 0L
  m[1:50, 1:2] <- 15L # 2 samples with >= 10 reads: excluded by ">= 10 in >= 3"
  keep <- rowSums(m >= 10) >= 3
  sf <- DESeq2::estimateSizeFactorsForMatrix(m[keep, ])
  nc <- check_normalization(m, t(t(m[keep, ]) / sf), min_counts = c(5, 10))
  expect_equal(nc$table$implied_size_factor, unname(sf))
  expect_equal(nc$table$sf_deseq2_provider_genes, unname(sf))
  expect_match(nc$facts$normalization_filter_rule, ">= 10 reads in >= 3 samples")
  expect_true(all(nc$table$implied_ratio_rel_iqr < 1e-12))
})
