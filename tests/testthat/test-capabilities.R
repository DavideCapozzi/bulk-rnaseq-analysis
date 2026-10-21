cfg_c <- function(n_tables = 1) {
  tabs <- list(
    list(id = "t1", adapter = "deseq2", file = "comparison.xlsx"),
    list(id = "t2", adapter = "noiseq", file = "noiseq_results.tsv")
  )[seq_len(n_tables)]
  make_cfg(list(
    project = list(name = "x", output_dir = "o"), inputs = list(contrast_tables = tabs),
    resources = list(items = list(kegg_pathways = list(source = "kegg_rest"))),
    contrast_pairs = if (n_tables == 2) list(list(id = "p", first = "t1", second = "t2")) else list()
  ))
}

test_that("DE-only input runs volcano and ORA, skips sample-level modules with the reason", {
  caps <- resolve_capabilities(cfg_c(1))
  st <- stats::setNames(caps$status, caps$module)
  expect_equal(st[["volcano"]], "run")
  expect_equal(st[["ora"]], "run")
  expect_equal(st[["contrasts"]], "skipped")
  expect_match(
    caps$reason[caps$module == "pca"],
    "pca skipped: needs data level\\(s\\) counts; input provides contrast_stats, ranked"
  )
  expect_match(
    caps$reason[caps$module == "concordance"],
    "concordance skipped: needs >= 2 contrasts with per-gene statistics; input provides 1 \\(comparison.xlsx\\)"
  )
  expect_match(caps$reason[caps$module == "qc"], "needs one of counts, qc")
  expect_equal(caps$depends_on[caps$module == "volcano"], "ingest")
})

test_that("two contrast tables and a pair enable concordance", {
  caps <- resolve_capabilities(cfg_c(2))
  expect_equal(caps$status[caps$module == "concordance"], "run")
})

test_that("descriptive contrasts do not enable significance-based modules; disabled modules say so", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("sample_id\tinput_column\tcondition", "a\tA\tc", "b\tB\tt", "c\tC\tc", "d\tD\tt"), f)
  cfg <- make_cfg(list(
    project = list(name = "x", output_dir = "o"), samples = list(sheet = f),
    inputs = list(counts = list(adapter = "matrix", file = "m.tsv", id_column = "id")),
    resources = list(items = list(gtf = list(source = "x"))),
    contrasts = list(list(id = "a_vs_b", numerator = list("a"), denominator = list("b"))),
    modules = list(markers = list(enabled = FALSE))
  ))
  caps <- resolve_capabilities(cfg, read_sample_sheet(f))
  expect_equal(caps$status[caps$module == "contrasts"], "run")
  expect_match(caps$reason[caps$module == "volcano"], "descriptive only")
  expect_equal(caps$status[caps$module == "markers"], "disabled")
  expect_match(caps$depends_on[caps$module == "concordance"] %||% "", "")
})
