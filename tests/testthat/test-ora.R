test_that("ORA tests GO and a cached KEGG table with gene-level inference scope", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed(organism_info("human")$orgdb)
  x <- read_de_table(fixture("scenario_c", "comparison.xlsx"), "c", "deseq2", "DESeq2_results")
  ctx <- list(
    cfg = list(annotation = list(id_type = "gene_symbol"), thresholds = list(padj = 0.05, abs_fold_change = 1.5)),
    organism = "human", resources = list(kegg_pathways = fixture("resources", "kegg_pathways.tsv.gz"))
  )
  params <- list(databases = c("GO_BP", "KEGG"), pvalue_cutoff = 1)
  res <- suppressMessages(run_ora(list(contrast_stats = x), params, ctx))
  expect_named(res$tables, c("ora_input_genes", "ora_go_bp", "ora_kegg"))
  expect_equal(nrow(res$tables$ora_input_genes), nrow(x))
  expect_gt(nrow(res$tables$ora_kegg), 0)
  expect_true(all(startsWith(res$tables$ora_kegg$term_id, "hsa990")))
  expect_equal(unique(res$tables$ora_kegg$inference_scope), "conditional_on_contrast")
  expect_error(
    run_ora(list(contrast_stats = x), list(databases = "KEGG"), list(cfg = ctx$cfg, organism = "human")),
    "kegg_pathways"
  )
})
