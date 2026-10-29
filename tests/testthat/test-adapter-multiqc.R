skip_if_no_python <- function() testthat::skip_if(!nzchar(Sys.which("python3")), "python3 not available")

samples_a <- function() {
  read_sample_sheet(fixture("scenario_a", "samples.tsv"), list(condition = c("control", "treated")))
}
report_a <- function(read = "r1") {
  list(list(id = "run", file = fixture("scenario_a", "multiqc", "multiqc_report.html"), read = read))
}
aliases_a <- list(A1 = "lib_A1", A2 = "lib_A2", B1 = "lib_B1", B2 = "lib_B2")
metric <- function(qc, m, s) qc$value[qc$metric == m & qc$sample_id == s]

# Expected values follow from the definitions in tests/fixtures/make_multiqc_fixture.py.
test_that("MultiQC metrics are derived from the STAR, RSeQC and raw FastQC sections", {
  skip_if_no_python()
  mq <- read_multiqc(report_a(), samples_a(), aliases_a)
  qc <- mq$qc
  expect_true(all(qc$source == "multiqc:run"))
  star_in <- 6e5 + 9e4 + 5e3 + 2.5e5 + 4e3
  expect_equal(metric(qc, "star_input_reads", "s1"), star_in)
  expect_equal(metric(qc, "pct_unique", "s1"), 100 * 6e5 / star_in)
  expect_equal(metric(qc, "raw_read_pairs", "s1"), 1e6)
  expect_equal(metric(qc, "pct_raw_to_star", "s1"), 100 * star_in / 1e6)
  expect_equal(metric(qc, "assigned_reverse", "s1"), 3.9e5)
  expect_equal(metric(qc, "strand_ratio_star", "s1"), 9e3 / 3.9e5)
  expect_equal(metric(qc, "pct_assigned_unique", "s1"), 100 * 3.9e5 / 6e5)
  expect_equal(metric(qc, "rseqc_antisense_pct", "s1"), 93)
  expect_equal(metric(qc, "pct_utr", "s1"), 100 * 8e4 / (3e5 + 2e4 + 6e4 + 8e4 + 2e3 + 3e3 + 1e4))
  expect_equal(metric(qc, "pct_intergenic_or_flank", "s1"), 100 * 1.5e4 / (3e5 + 2e4 + 6e4 + 8e4 + 1.5e4))
  expect_equal(metric(qc, "genebody_3to5", "s1"), (0.4 + 0.006 * 80.5) / (0.4 + 0.006 * 20.5))
  expect_equal(metric(qc, "adapter_pct_last_pos_r1", "s2"), 0.02 * 40 * 2)
  expect_equal(mq$adapter_positions, list(run = 40))
  q <- 20:41
  expect_equal(metric(qc, "mean_q_r1", "s1"), sum(q * (q - 20)) / sum(q - 20))
  expect_equal(metric(qc, "dup_lvl1_pct_r1", "s2"), 55)
  expect_equal(metric(qc, "dup_ge100_pct_r1", "s2"), 5 * 2)
  # The QC-only library: raw FastQC only, low quality from position 37 (38 - 0.5 x < 20).
  expect_equal(metric(qc, "raw_read_pairs", "s5"), 4e4)
  expect_equal(metric(qc, "first_pos_meanq_lt20_r1", "s5"), 37)
  expect_true(is.na(metric(qc, "first_pos_meanq_lt20_r1", "s1")))
  expect_length(metric(qc, "star_input_reads", "s5"), 0)
  expect_setequal(basename(mq$files), c("run_plotdata.json", "run_general_stats.tsv", "plot_index.tsv"))
})

test_that("FastQC sections are read only from reports marked with `read`", {
  skip_if_no_python()
  qc <- read_multiqc(report_a(read = NULL), samples_a(), aliases_a)$qc
  expect_false(any(grepl("_r1$|^raw_read_pairs$", qc$metric)))
  expect_true("pct_unique" %in% qc$metric)
})

test_that("MultiQC sample names must match exactly one sample", {
  skip_if_no_python()
  expect_error(read_multiqc(report_a(), samples_a()), "MultiQC sample 'A1' matches no input_column")
  ss <- samples_a()
  ss$input_column[ss$sample_id == "s5"] <- "lib_"
  expect_error(read_multiqc(report_a(), ss, aliases_a), "matches lib_A1, lib_")
})

test_that("the qc input level comes from MultiQC reports, and inputs.qc is accepted", {
  cfg <- list(inputs = list(qc = list(multiqc = list(list(id = "run", file = "r.html")))))
  expect_equal(input_levels(cfg)$qc, "MultiQC (run)")
  expect_true(validate_config(c(cfg, list(samples = list(sheet = "s.tsv")))))
})

test_that("an empty factor is allowed only for samples without expression data", {
  ss <- samples_a()
  expect_true(is.na(ss$condition[ss$sample_id == "s5"]))
  f <- tempfile(fileext = ".tsv")
  writeLines(c("sample_id\tinput_column\tcondition", "a\tA\tctrl", "b\tB\t"), f)
  expect_error(read_sample_sheet(f, list(condition = "ctrl")), "'condition' is empty for sample\\(s\\) with expression data: b")
})

test_that("ingest writes the decoded MultiQC files and passes QC-only samples to qc_summary", {
  skip_if_no_python()
  dir <- withr::local_tempdir()
  cfg <- make_cfg(list(
    project = list(name = "x", output_dir = "out"),
    samples = list(sheet = normalizePath(fixture("scenario_a", "samples.tsv"))),
    inputs = list(
      base_dir = normalizePath(fixture("scenario_a")),
      qc = list(multiqc = report_a()[1], sample_aliases = aliases_a)
    )
  ), dir)
  cfg$inputs$qc$multiqc[[1]]$file <- "multiqc/multiqc_report.html"
  ing <- run_ingest(cfg)
  expect_equal(ing$report$multiqc_adapter_position, list(run = 40))
  write_ingest_outputs(ing, file.path(dir, "out"))
  mq <- file.path(dir, "out", "01_ingest", "multiqc")
  expect_setequal(list.files(mq), c("run_plotdata.json", "run_general_stats.tsv", "plot_index.tsv"))
  res <- run_qc(list(qc = ing$qc), list(), list(cfg = cfg, gene_lists = list(), organism = "human"))
  expect_equal(res$tables$qc_summary$sample_id, paste0("s", 1:5))
  expect_registry_outputs(res, "qc")
})
