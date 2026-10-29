# Two donors (d1, d2), each sorted into two fractions (neg, pos). Gene values are chosen so that CPM and the
# expected distances can be computed by hand: every library has 1e6 reads.
pairing_se <- function() {
  genes <- c("XIST", "RPS4Y1", "DDX3Y", "RPS4X", "HLA-A", "GSTM1", "FILLER")
  m <- cbind(
    s1 = c(100, 300, 100, 1000, 63, 7, 0),
    s2 = c(100, 300, 100, 1000, 31, 3, 0),
    s3 = c(0, 800, 200, 1000, 255, 1, 0),
    s4 = c(0, 800, 200, 1000, 127, 0, 0)
  )
  m[length(genes), ] <- 1e6 - colSums(m)
  rownames(m) <- paste0("g", seq_along(genes))
  storage.mode(m) <- "integer"
  rd <- S4Vectors::DataFrame(
    feature_id = rownames(m), gene_symbol = genes, biotype = "protein_coding",
    chromosome = c("X", "Y", "Y", "X", "6", "1", "1"), row.names = rownames(m)
  )
  cd <- S4Vectors::DataFrame(
    sample_id = colnames(m), donor = c("d1", "d1", "d2", "d2"), fraction = c("neg", "pos", "neg", "pos"),
    row.names = colnames(m)
  )
  SummarizedExperiment::SummarizedExperiment(list(counts = m), rowData = rd, colData = cd)
}
pairing_context <- function() {
  list(cfg = list(), organism = "human", gene_lists = read_gene_lists())
}

test_that("pairing reports sex genes, fingerprint distances by pair type and the chance level of the check", {
  params <- list(
    group_column = "donor", across_factor = "fraction",
    fingerprint_variants = list(list(id = "fp", genes = c("HLA-A", "GSTM1"))),
    male_reference_samples = c("s3", "s4"),
    sex_denominators = list(list(id = "chr_y"), list(id = "rps4", numerator = "RPS4Y1", denominator = "RPS4X"))
  )
  res <- run_pairing(list(se = pairing_se()), params, pairing_context())
  expect_registry_outputs(res, "pairing")
  s <- res$tables$sex_summary
  expect_equal(s$xist_cpm, c(100, 100, 0, 0))
  expect_equal(s$chr_y_cpm_sum, c(400, 400, 1000, 1000))
  d <- res$tables$fingerprint_distances
  # log2(CPM + 1): HLA-A 6, 5, 8, 7; GSTM1 3, 2, 1, 0.
  expect_equal(d$distance[d$sample_a == "s1" & d$sample_b == "s2"], 1)
  expect_equal(d$distance[d$sample_a == "s1" & d$sample_b == "s4"], (1 + 3) / 2)
  expect_equal(d$pair_type[d$sample_a == "s1" & d$sample_b == "s3"], "different_group_same_level")
  expect_equal(d$pair_type[d$sample_a == "s1" & d$sample_b == "s4"], "different_group_different_level")
  ps <- res$tables$pairing_summary
  expect_true(ps$same_group_pairs_nearest_across)
  expect_equal(ps$chance_at_most, 1 / 2)
  mf <- res$tables$male_fraction
  expect_equal(mf$male_fraction_chr_y, c(0.4, 0.4, 1, 1))
  expect_equal(mf$male_fraction_rps4, c(300 / 800, 300 / 800, 1, 1))
  expect_equal(mf$xist_over_rps4, c(0.1, 0.1, 0, 0))
  g <- res$tables$male_fraction_genes
  expect_equal(g$male_fraction[g$gene_symbol == "DDX3Y" & g$sample_id == "s1"], 0.5)
  expect_match(res$meta$pairing_summary$interpretation, "does\\s+not resolve it")
})

test_that("without male reference samples the male fraction is not estimated, and the meta says so", {
  res <- run_pairing(list(se = pairing_se()), list(), pairing_context())
  expect_null(res$tables$male_fraction)
  expect_null(res$tables$pairing_summary)
  expect_match(res$meta$sex_summary$notes, "male fraction not estimated", all = FALSE)
  expect_equal(unique(res$tables$fingerprint_distances$variant_id), "fingerprint")
  expect_equal(unique(res$tables$fingerprint_distances$n_genes), 2) # the default list minus absent genes
})

test_that("pairing skips variants without genes, says so, and stops only when nothing can be computed", {
  params <- list(fingerprint_variants = list(list(id = "none", genes = "ABSENT1"), list(id = "fp", genes = "GSTM1")))
  res <- run_pairing(list(se = pairing_se()), params, pairing_context())
  expect_equal(unique(res$tables$fingerprint_distances$variant_id), "fp")
  expect_match(res$meta$fingerprint_distances$notes, "without genes among the features: none", all = FALSE)
  params$fingerprint_variants <- list(list(id = "none", genes = "ABSENT1"))
  params$sex_genes <- "ABSENT2"
  expect_error(run_pairing(list(se = pairing_se()), params, pairing_context()), "neither the sex genes")
})

test_that("Hart sample sizes match the reference values", {
  expect_equal(n_hart(0.2, 20, 1.5, 0.001), 19)
  expect_equal(n_hart(0.2, 20, 1.5, 0.01), 13)
})

test_that("power reports the BCV, its per-gene interval for the residual df, and the sample size grid", {
  se <- nb_se(groups = c("a", "b", "a", "b"))
  SummarizedExperiment::colData(se)$h <- c("x", "x", "y", "y")
  cfg <- list(designs = list(list(id = "add", formula = "~ g + h")), gene_filter = list(min_cpm = 1, min_samples = 2))
  params <- list(bootstrap_replicates = 5, grid = list(cv_labels = c("w", "x", "y", "z")))
  res <- run_power(list(se = se), params, list(cfg = cfg, organism = "human", gene_lists = list(), seed = 3))
  expect_registry_outputs(res, "power")
  b <- res$tables$bcv
  expect_equal(b$residual_df, 1)
  expect_equal(round(c(b$gene_bcv_ci95_factor_low, b$gene_bcv_ci95_factor_high), c(3, 1)), c(0.446, 31.9))
  expect_true(b$common_bcv_boot_q025 <= b$common_bcv_boot_q975)
  again <- run_power(list(se = se), params, list(cfg = cfg, organism = "human", gene_lists = list(), seed = 3))
  expect_identical(again$tables$bcv, b)
  grid <- res$tables$sample_size_grid
  expect_equal(nrow(grid), 4 * 2 * 2 * 2)
  expect_equal(grid$scenario[grid$cv == 0.6][[1]], "z")
  expect_false(is.unsorted(grid$cv))
  expect_error(run_power(list(se = se), list(design = "zz"), list(cfg = cfg)), "design 'zz' not found")
})
