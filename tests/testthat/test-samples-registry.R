test_that("sample sheet sets factor levels and parses inputs_available", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("sample_id\tinput_column\tcondition\tinputs_available", "a\tA\tctrl\t", "b\tB\ttrt\tqc,counts"), f)
  ss <- read_sample_sheet(f, list(condition = c("ctrl", "trt")))
  expect_equal(levels(ss$condition), c("ctrl", "trt"))
  expect_true(is.na(ss$inputs_available[[1]]))
  expect_equal(ss$inputs_available[[2]], c("qc", "counts"))
  expect_equal(nrow(bulkrnaseq:::samples_with_level(ss, "abundance")), 1)
})

test_that("sample sheet errors name the problem", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("sample_id\tinput_column\tcondition", "a\tA\tctrl", "a\tB\ttrt"), f)
  expect_error(read_sample_sheet(f), "duplicated sample_id")
  writeLines(c("sample_id\tinput_column\tcondition", "a\tA\tctrl", "b\tB\tother"), f)
  expect_error(read_sample_sheet(f, list(condition = c("ctrl", "trt"))), "values not in samples.factors: other")
  writeLines(c("sample_id\tcondition", "a\tctrl"), f)
  expect_error(read_sample_sheet(f), "missing column\\(s\\) 'input_column'")
})

test_that("registry and organism table are consistent", {
  reg <- read_registry()
  expect_false(anyDuplicated(vapply(reg, function(m) m$number, 1L)) > 0)
  expect_equal(bulkrnaseq:::module_dir_name("contrasts", reg), "09_contrasts")
  schema <- yaml::read_yaml(config_schema_path())
  expect_setequal(names(reg), names(schema$properties$modules$properties))
  expect_equal(organism_info("human")$chromosome_mito, "MT")
  expect_error(organism_info("mouse"), "not in organisms.yaml")
})

test_that("organism-specific identifiers appear only in the organism table", {
  roots <- c(testthat::test_path("..", "..", "R"), testthat::test_path("..", "..", "workflow"))
  roots <- roots[dir.exists(roots)]
  skip_if(length(roots) == 0, "source tree not available (installed package)")
  files <- list.files(roots, pattern = "(\\.(R|py|smk|yaml|tsv)|Snakefile)$", recursive = TRUE, full.names = TRUE)
  hits <- unlist(lapply(files, function(f) {
    x <- readLines(f, warn = FALSE)
    f[any(grepl("org\\.Hs\\.eg\\.db|\"hsa\"|\"MT-\"", x))]
  }))
  expect_length(hits, 0)
})
