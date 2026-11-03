axis_fixture <- function() {
  set.seed(2)
  n <- 120
  ids <- sprintf("ENSG%011d", seq_len(n))
  # Cell "A" high on the first 40 genes, "B1"/"B2" high on genes 41-80; the rest similar.
  a <- c(rep(100, 40), rep(5, 40), rep(20, 40))
  b1 <- c(rep(5, 40), rep(100, 40), rep(20, 40))
  b2 <- c(rep(3, 40), rep(60, 40), rep(20, 40))
  hpa <- data.frame(
    Gene = rep(ids, 3), `Gene name` = "", `Immune cell` = rep(c("A", "B1", "B2"), each = n),
    TPM = 0, pTPM = 0, nTPM = c(a, b1, b2), check.names = FALSE
  )
  hf <- tempfile(fileext = ".tsv")
  data.table::fwrite(hpa, hf, sep = "\t")
  gs <- data.frame(
    gs_collection = "C7", gs_subcollection = "IMMUNESIGDB", gs_name = rep(c("UP_SET", "DN_SET"), each = 20),
    ensembl_gene = ids[c(1:20, 41:60)], gene_symbol = "", db_version = "v1"
  )
  gf <- tempfile(fileext = ".tsv")
  data.table::fwrite(gs, gf, sep = "\t")
  lfc1 <- c(rep(1, 40), rep(-1, 40), rnorm(40, sd = 0.2)) + rnorm(n, sd = 0.1)
  lfc2 <- c(rep(1, 40), rep(1, 40), rnorm(40, sd = 0.2)) + rnorm(n, sd = 0.1)
  ct <- rbind(
    new_contrast_table("c1", ids, lfc1, method = "m", source = "computed", extra = data.frame(ave_log2cpm = 1:n)),
    new_contrast_table("c2", ids, lfc2, method = "m", source = "computed", extra = data.frame(ave_log2cpm = 1:n))
  )
  counts <- matrix(rpois(n * 4, 200), n, 4, dimnames = list(ids, c("s1", "s2", "s3", "s4")))
  counts[1:40, c("s1", "s2")] <- counts[1:40, c("s1", "s2")] * 4
  se <- toy_se(counts, rd = data.frame(
    gene_symbol = paste0("G", seq_len(n)), biotype = "protein_coding", chromosome = "1", row.names = ids
  ), cd = data.frame(sample_id = colnames(counts), group = c("x", "x", "y", "y"), row.names = colnames(counts)))
  cfg <- list(
    resources = list(items = list(hpa = list(kind = "hpa_immune_cell"), sets = list(kind = "msigdb_export"))),
    gene_filter = list(min_cpm = 0, min_samples = 1, biotypes = "protein_coding", sex_genes = "keep"),
    contrasts = list(
      list(id = "c1", numerator = "group=x", denominator = "group=y"),
      list(id = "c2", numerator = "s1", denominator = "s3")
    ),
    contrast_pairs = list(list(id = "p", first = "c1", second = "c2"))
  )
  ctx <- list(cfg = cfg, organism = "human", resources = list(hpa = hf, sets = gf))
  list(ct = ct, se = se, ctx = ctx, ids = ids)
}

axis_params <- list(
  axes = list(list(id = "a_vs_b", resource = "hpa", numerator = "A", denominator = c("B1", "B2"), threshold = 2)),
  gene_sets = list(resource = "sets", sets = c("UP_SET", "DN_SET")),
  sign_test = list(threshold = 0.5)
)

test_that("reference_axis reports correlations, set means, sign concordance and sample scores", {
  x <- axis_fixture()
  res <- run_reference_axis(list(contrast_stats = x$ct, se = x$se), axis_params, x$ctx)
  expect_registry_outputs(res, "reference_axis")
  t <- res$tables$reference_axis
  sp <- t[t$measure == "spearman_lfc_vs_axis", ]
  expect_equal(sp$contrast_id, c("c1", "c2"))
  expect_gt(sp$value[sp$contrast_id == "c1"], 0.5)
  expect_equal(sp$n, c(120L, 120L))
  m <- t[t$measure == "mean_lfc_in_set", ]
  expect_equal(nrow(m), 4)
  lfc1 <- x$ct$lfc[x$ct$contrast_id == "c1"]
  expect_equal(m$value[m$set_id == "UP_SET" & m$contrast_id == "c1"], mean(lfc1[1:20]))
  sg <- t[t$measure == "sign_concordance", ]
  lfc2 <- x$ct$lfc[x$ct$contrast_id == "c2"]
  both <- abs(lfc1) > 0.5 & abs(lfc2) > 0.5
  k <- sum(both & sign(lfc1) == sign(lfc2))
  expect_equal(sg$n, sum(both))
  expect_equal(sg$n_concordant, k)
  expect_equal(sg$pvalue, stats::binom.test(k, sum(both))$p.value)
  expect_equal(sg$inference_scope, "conditional_on_contrast")
  ax <- res$tables$axis_genes
  expect_equal(ax$axis[[1]], log2(101 / (mean(c(5, 3)) + 1)))
  expect_equal(sum(ax$class == "high"), 40)
  expect_equal(sum(ax$class == "low"), 40)
  sc <- t[t$measure == "sample_axis_score", ]
  expect_equal(sc$sample_id, c("s1", "s2", "s3", "s4"))
  expect_equal(unique(sc$n_high), 40L)
  lc <- log_cpm(SummarizedExperiment::assay(x$se, "counts"))
  expect_equal(sc$value, unname(colMeans(lc[1:40, ]) - colMeans(lc[41:80, ])))
  sh <- t[t$measure == "axis_score_shift", ]
  expect_equal(sh$value[sh$contrast_id == "c1"], mean(sc$value[1:2]) - mean(sc$value[3:4]))
  expect_equal(sh$value[sh$contrast_id == "c2"], sc$value[1] - sc$value[3])
  write_module_outputs(res, tempfile(), list())
})

test_that("reference_axis honours subsets and runs without counts", {
  x <- axis_fixture()
  flags <- data.frame(feature_id = x$ids, masked = seq_along(x$ids) <= 10)
  p <- utils::modifyList(axis_params, list(subset = "no_masked", set_subset = "no_masked+ge_median"))
  res <- run_reference_axis(list(contrast_stats = x$ct, gene_flags = flags), p, x$ctx)
  t <- res$tables$reference_axis
  expect_equal(unique(t$n[t$measure == "spearman_lfc_vs_axis"]), 110L)
  # ge_median keeps ave_log2cpm >= 60.5, so neither set has genes.
  expect_true(all(t$n[t$measure == "mean_lfc_in_set"] == 0))
  expect_false(any(t$measure %in% c("sample_axis_score", "axis_score_shift")))
  expect_equal(nrow(res$tables$axis_genes), 110)
})

test_that("reference_axis fails loudly on bad configuration", {
  x <- axis_fixture()
  run <- function(p, ctx = x$ctx) run_reference_axis(list(contrast_stats = x$ct), p, ctx)
  expect_error(run(list(gene_sets = list(sets = "UP_SET"))), "needs gene_sets.resource")
  expect_error(run(list(gene_sets = list(resource = "sets", sets = "NOPE"))), "not found: NOPE")
  expect_error(run(list(axes = list(list(id = "a", resource = "zz", numerator = "A", denominator = "B1")))), "'zz'")
  expect_error(run(list(axes = list(list(id = "a", resource = "hpa", numerator = "Z", denominator = "B1")))), "not found: Z")
  expect_error(run(list(axes = rep(axis_params$axes, 2))), "duplicated axis")
  ctx <- x$ctx
  ctx$cfg$contrast_pairs <- NULL
  expect_error(run(list(sign_test = list(threshold = 1)), ctx), "needs contrast_pairs")
  ctx <- x$ctx
  ctx$cfg$resources$items$hpa$kind <- "tsv"
  expect_error(run(axis_params, ctx), "needs kind hpa_immune_cell")
})
