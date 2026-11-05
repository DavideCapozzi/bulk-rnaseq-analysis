markers_se <- function() {
  m <- matrix(
    c(
      100L, 5L, 895L, 0L,
      10L, 50L, 930L, 10L,
      300L, 1L, 699L, 0L
    ), 4,
    dimnames = list(c("E1", "E2", "E3", "E4"), c("s1", "s2", "s3"))
  )
  rd <- S4Vectors::DataFrame(
    feature_id = rownames(m), gene_symbol = c("FOXP3", "IL7R", "ACTB", "IL7R"),
    biotype = c("protein_coding", "protein_coding", "protein_coding", "lncRNA"), chromosome = "1"
  )
  cd <- S4Vectors::DataFrame(
    sample_id = colnames(m), fraction = c("pos", "neg", "pos"), label = c("S 1", "S 2", "S 3"),
    row.names = colnames(m)
  )
  se <- toy_se(m, rd, cd)
  SummarizedExperiment::assay(se, "tpm") <- m * 2
  se
}

markers_ctx <- list(
  cfg = list(samples = list(factors = list(fraction = c("neg", "pos")), label_column = "label")),
  gene_lists = list(treg_core = c("FOXP3", "IL7R", "MISSING1"))
)

test_that("markers reports CPM, TPM, panel scores and ratios", {
  se <- markers_se()
  p <- list(
    panels = list(list(id = "treg", genes = "treg_core"), list(id = "house", genes = list("ACTB"))),
    biotypes = list("protein_coding"),
    ratios = list(list(
      id = "foxp3_pos_neg", numerator = list(gene = "FOXP3", sample = "s1"),
      denominator = list(gene = "FOXP3", sample = "s2")
    ))
  )
  res <- run_markers(list(se = se), p, markers_ctx)
  expect_registry_outputs(res, "markers")
  e <- res$tables$marker_expression
  expect_equal(nrow(e), 3 * 3)
  expect_equal(unique(e$gene_symbol), c("FOXP3", "IL7R", "ACTB"))
  expect_equal(e$cpm[e$gene_symbol == "FOXP3"], c(100, 10, 300) / 1000 * 1e6)
  expect_equal(e$tpm, e$count * 2)
  expect_false("normalized" %in% names(e))
  expect_equal(res$meta$marker_expression$genes_not_found, "treg:MISSING1")
  s <- res$tables$marker_panel_scores
  expect_equal(s$cpm_sum[s$panel_id == "treg"], c(105, 60, 301) / 1000 * 1e6)
  expect_equal(s$n_genes, rep(c(2L, 1L), each = 3))
  r <- res$tables$marker_ratios
  expect_equal(r$log2_ratio, log2((1e5 + 1) / (1e4 + 1)))
  expect_equal(r$label, "foxp3_pos_neg")
  expect_s3_class(res$figures$markers_heatmap, "ggplot")
  expect_s3_class(res$figures$marker_panels, "ggplot")
})

test_that("markers keeps every feature of a symbol without a biotype filter", {
  res <- run_markers(list(se = markers_se()), list(panels = list(list(id = "p", genes = "IL7R"))), markers_ctx)
  expect_equal(unique(res$tables$marker_expression$feature_id), c("E2", "E4"))
  expect_equal(res$tables$marker_panel_scores$n_genes, rep(2L, 3))
})

test_that("markers fails loudly on bad configuration", {
  se <- markers_se()
  run <- function(p) run_markers(list(se = se), p, markers_ctx)
  expect_error(run(list()), "no panels configured")
  expect_error(run(list(panels = list(list(id = "a", genes = "ACTB"), list(id = "a", genes = "ACTB")))), "duplicated")
  expect_error(run(list(panels = list(list(id = "a", genes = "NOPE")))), "no panel gene matches")
  ratio <- function(gene, sample) {
    list(
      panels = list(list(id = "a", genes = "ACTB")),
      ratios = list(list(id = "r", numerator = list(gene = gene, sample = sample), denominator = list(
        gene = "ACTB", sample = "s1"
      )))
    )
  }
  expect_error(run(ratio("IL7R", "s1")), "matches 2 features")
  expect_error(run(ratio("ACTB", "zz")), "unknown sample_id 'zz'")
})

test_that("markers runs only with panels, and reads the normalized assay when present", {
  cfg <- make_cfg(list(
    project = list(name = "x", output_dir = "o"), inputs = list(counts = list(adapter = "matrix", file = "c.tsv")),
    resources = list(items = list())
  ))
  caps <- resolve_capabilities(cfg)
  expect_match(caps$reason[caps$module == "markers"], "needs modules.markers.panels")
  cfg$modules$markers$panels <- list(list(id = "a", genes = "ACTB"))
  caps <- resolve_capabilities(cfg)
  expect_equal(caps$status[caps$module == "markers"], "run")
  se <- markers_se()
  SummarizedExperiment::assay(se, "normalized") <- SummarizedExperiment::assay(se, "counts") / 2
  res <- run_markers(list(se = se), list(panels = list(list(id = "a", genes = "ACTB"))), markers_ctx)
  expect_equal(res$tables$marker_expression$normalized, c(895, 930, 699) / 2)
})
