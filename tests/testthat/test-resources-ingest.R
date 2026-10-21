test_that("resources are fetched once, recorded with md5 and verified before use", {
  src <- tempfile(fileext = ".tsv")
  writeLines("a\tb", src)
  cfg <- make_cfg(list(
    project = list(name = "x", output_dir = "o"), inputs = list(),
    resources = list(cache_dir = "cache", items = list(demo = list(source = paste0("file://", src))))
  ))
  expect_error(verify_resources(cfg), "run the explicit fetch step")
  fetch_resources(cfg)
  p <- verify_resources(cfg)
  expect_true(file.exists(p[["demo"]]))
  writeLines("changed", p[["demo"]])
  expect_error(verify_resources(cfg), "md5 mismatch")
  cfg$resources$items$demo$md5 <- strrep("0", 32)
  expect_error(fetch_resources(cfg), "differs from the declared")
})

test_that("ingest of the DE-only scenario produces canonical contrast tables with mapping facts", {
  dir <- tempfile("proj_")
  dir.create(dir)
  cfg <- make_cfg(
    list(
      project = list(name = "x", output_dir = "out"),
      inputs = list(
        base_dir = normalizePath(fixture("scenario_c")),
        contrast_tables = list(
          list(id = "c1", adapter = "deseq2", file = "comparison.xlsx", sheet = "DESeq2_results"),
          list(id = "n1", adapter = "noiseq", file = "noiseq_results.tsv")
        )
      ),
      annotation = list(id_type = "gene_symbol"),
      resources = list(
        cache_dir = "cache",
        items = list(gtf = list(source = paste0("file://", normalizePath(fixture("resources", "annotation.gtf.gz")))))
      )
    ),
    dir
  )
  fetch_resources(cfg)
  ing <- run_ingest(cfg)
  expect_equal(sort(unique(ing$contrast_stats$contrast_id)), c("c1", "n1"))
  expect_equal(ing$contrast_stats$gene_symbol, ing$contrast_stats$feature_id)
  expect_gt(ing$report$contrast_symbols_in_annotation, 0)
  written <- write_ingest_outputs(ing, file.path(dir, "out"))
  expect_true(file.exists(file.path(dir, "out", "01_ingest", "contrasts.tsv.gz")))
  expect_true(file.exists(file.path(dir, "out", "01_ingest", "ingest_report.md")))
})
