test_that("methods templates substitute values and keep conditional blocks only when set", {
  r <- function(t, p = list(), cfg = list(), x = list()) render_methods_template(t, cfg, p, x)
  expect_equal(r("n = {{p.n}}", list(n = 3)), "n = 3")
  expect_equal(r("{{p.v}}", list(v = list("a", "b"))), "a, b")
  expect_equal(r("a{{#p.flag}} flagged{{/p.flag}}.", list(flag = FALSE)), "a.")
  expect_equal(r("a{{#p.flag}} flagged{{/p.flag}}.", list(flag = TRUE)), "a flagged.")
  expect_equal(r("{{#x.k}}k={{x.k}}{{/x.k}}", x = list(k = "v")), "k=v")
  expect_equal(r("seed {{project.seed}}", cfg = list(project = list(seed = 7))), "seed 7")
  expect_equal(r("{{p.a.b}}", list(a = list(b = 0.05))), "0.05")
  expect_match(r("{{pkg.edgeR}}"), "^[0-9.]+$")
  expect_error(r("{{p.missing}}"), "no value for 'p.missing'")
  expect_error(r("{{#p.a}} open"), "is not closed")
})

test_that("every implemented module has a methods template with known references", {
  reg <- read_registry()
  implemented <- names(reg)[vapply(reg, module_implemented, TRUE)]
  for (m in c("inputs", implemented)) {
    tpl <- methods_template(m)
    expect_false(is.null(tpl), info = m)
    expect_silent(methods_references(tpl$refs))
  }
  refs <- yaml::read_yaml(package_file("templates", "methods", "references.yaml"))
  expect_true(all(nzchar(unlist(refs))))
})

test_that("methods phrases describe inputs, filter, designs and thresholds", {
  cfg <- make_cfg(list(
    project = list(name = "x", output_dir = "o"),
    inputs = list(
      counts = list(adapter = "matrix", file = "c.xlsx"),
      contrast_tables = list(list(id = "t", adapter = "deseq2", file = "t.xlsx", df_residual = 4))
    ),
    annotation = list(id_type = "gene_symbol"),
    designs = list(list(id = "d", formula = "~ group", contrast = list("group", "b", "a")))
  ))
  x <- methods_phrases(cfg)
  expect_equal(x$id_type, "gene symbols")
  expect_equal(x$contrast_tables, "t (deseq2 table, 4 residual df)")
  expect_equal(x$designs, "d (~ group; group b vs a)")
  expect_match(x$gene_filter, "at least 2 counts per million in at least 2 samples; chromosome Y genes and XIST")
  expect_match(x$de_thresholds, "\\|log2 fold change\\| >= 0.585")
})

# A results folder with ingest objects, a DE table and markers outputs, as the workflow writes them.
workbook_run <- function() {
  out <- tempfile("run_")
  dir.create(file.path(out, "00_provenance"), recursive = TRUE)
  se <- markers_like_se()
  dir.create(file.path(out, "01_ingest", "objects"), recursive = TRUE)
  saveRDS(se, file.path(out, "01_ingest", "objects", "se.rds"))
  ct <- new_contrast_table("b_vs_a", rownames(se), c(2, -2, 0.1),
    gene_symbol = c("FOXP3", "IL7R", "ACTB"), lfc_se = 0.5, stat = c(4, -4, 0.2), pvalue = c(1e-5, 1e-4, 0.8),
    padj = c(1e-4, 1e-3, 0.9), base_mean = 100, significance_metric = "padj", method = "imported_deseq2",
    source = "imported:t.xlsx", inference_scope = "sample"
  )
  saveRDS(ct, file.path(out, "01_ingest", "objects", "contrast_stats.rds"))
  caps <- data.frame(
    module = c("markers", "volcano"), number = c(4L, 14L), dir = c("04_markers", "14_volcano"),
    status = c("run", "skipped"), reason = c("", "x"), depends_on = c("ingest", ""), stringsAsFactors = FALSE
  )
  write_json(caps, file.path(out, "00_provenance", "capabilities.json"))
  panels <- list(list(id = "treg", genes = c("FOXP3", "IL7R")))
  res <- run_markers(list(se = se), list(panels = panels), list())
  write_module_outputs(res, file.path(out, "04_markers"))
  sheet <- tempfile(fileext = ".tsv")
  data.table::fwrite(data.frame(
    sample_id = colnames(se), input_column = colnames(se), group = c("a", "a", "b"), label = c("A1", "A2", "B1")
  ), sheet, sep = "\t")
  cfg <- make_cfg(list(
    project = list(name = "demo", output_dir = out), inputs = list(),
    samples = list(sheet = sheet, factors = list(group = c("a", "b")), label_column = "label"),
    modules = list(markers = list(panels = panels)), export = list(researcher_workbook = list(enabled = TRUE))
  ))
  list(out = out, cfg = cfg)
}

markers_like_se <- function() {
  m <- matrix(c(100L, 5L, 895L, 10L, 50L, 940L, 300L, 1L, 699L), 3,
    dimnames = list(c("E1", "E2", "E3"), c("s1", "s2", "s3"))
  )
  rd <- S4Vectors::DataFrame(feature_id = rownames(m), gene_symbol = c("FOXP3", "IL7R", "ACTB"), chromosome = "1")
  cd <- S4Vectors::DataFrame(
    sample_id = colnames(m), group = c("a", "a", "b"), label = c("A1", "A2", "B1"),
    row.names = colnames(m)
  )
  toy_se(m, rd, cd)
}

test_that("the researcher workbook collects DE and module tables with a README and a glossary", {
  x <- workbook_run()
  write_methods(x$cfg, x$out)
  path <- write_researcher_workbook(x$cfg, x$out)
  expect_equal(basename(path), "results_workbook.xlsx")
  sheets <- openxlsx::getSheetNames(path)
  expect_equal(sheets, c(
    "README", "Methods", "DE b_vs_a", "UP b_vs_a", "DOWN b_vs_a", "Markers CPM", "Marker panels"
  ))
  de <- openxlsx::read.xlsx(path, "DE b_vs_a")
  expect_equal(de$gene_symbol, c("FOXP3", "IL7R", "ACTB"))
  expect_equal(de$fold_change, c(4, -4, 2^0.1))
  expect_equal(de$class, c("UP", "DOWN", "NS"))
  mk <- openxlsx::read.xlsx(path, "Markers CPM", check.names = FALSE)
  expect_equal(names(mk), c("panel_id", "gene_symbol", "feature_id", "A1", "A2", "B1"))
  readme <- openxlsx::read.xlsx(path, "README")
  expect_true(all(c("fold_change", "panel_id") %in% readme$item))
  expect_true("UP b_vs_a" %in% readme$item)
  first <- tools::md5sum(path)
  write_researcher_workbook(x$cfg, x$out)
  expect_equal(tools::md5sum(path), first)
  methods <- readLines(file.path(x$out, "report", "methods.md"))
  expect_true("## 04_markers" %in% methods)
  expect_false(any(grepl("14_volcano", methods)))
})

test_that("the researcher workbook refuses columns without a glossary entry", {
  x <- workbook_run()
  f <- file.path(x$out, "04_markers", "tables", "marker_ratios.tsv")
  data.table::fwrite(data.frame(ratio_id = "r", undocumented_column = 1), f, sep = "\t")
  expect_error(write_researcher_workbook(x$cfg, x$out), "without a glossary entry: undocumented_column")
})
