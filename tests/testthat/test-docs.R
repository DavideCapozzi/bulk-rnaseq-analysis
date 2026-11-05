# Generated and hand-written documentation that must follow the code. docs/, README.md and config/ are not in the
# built package, so these checks run from the source tree (devtools::test(), CI) and are skipped by R CMD check.
repo_file <- function(...) test_path("..", "..", ...)

test_that("docs/MODULES.md and docs/CONFIG.md are up to date with the registry and the schema", {
  path <- repo_file("docs", "MODULES.md")
  skip_if_not(file.exists(path), "docs/ is not in the built package")
  expect_identical(readLines(path), modules_doc(), info = "regenerate with: Rscript tools/make_docs.R")
  expect_identical(readLines(repo_file("docs", "CONFIG.md")), config_doc(), info = "regenerate with: Rscript tools/make_docs.R")
})

test_that("every parameter of the config schema has a description", {
  expect_identical(schema_params_without_description(), character())
})

test_that("schema_param_rows() renders nullable objects, arrays of objects and maps of objects", {
  s <- list(properties = list(
    a = list(type = list("object", "null"), properties = list(b = list(type = "string", description = "B."))),
    l = list(type = "array", items = list(type = "object", properties = list(c = list(type = "integer")))),
    m = list(type = "object", additionalProperties = list(type = "object", properties = list(d = list(type = "string"))))
  ))
  rows <- schema_param_rows(s, s)
  expect_identical(sub("^\\| `([^`]+)`.*$", "\\1", rows), c("a", "a.b", "l", "l[].c", "m", "m.<key>.d"))
  expect_identical(schema_params_without_description(s), c("a", "l", "l[].c", "m", "m.<key>.d"))
})

test_that("the README lists the modules each example scenario runs", {
  readme <- repo_file("README.md")
  skip_if_not(file.exists(readme), "README.md is not in the built package")
  rows <- grep("^\\| `scenario_[a-z_]+\\.yaml` \\|", readLines(readme), value = TRUE)
  examples <- list.files(repo_file("config", "examples"), pattern = "^scenario_.*\\.yaml$")
  configs <- gsub("`", "", vapply(strsplit(rows, "|", fixed = TRUE), function(cells) trimws(cells[[2]]), ""))
  expect_setequal(configs, examples)
  for (i in seq_along(rows)) {
    cells <- strsplit(rows[[i]], "|", fixed = TRUE)[[1]]
    listed <- regmatches(cells[[length(cells)]], gregexpr("(?<=`)[a-z_]+(?=`)", cells[[length(cells)]], perl = TRUE))
    cfg <- read_config(repo_file("config", "examples", configs[[i]]))
    samples <- if (!is.null(cfg$samples$sheet)) {
      read_sample_sheet(config_path(cfg, cfg$samples$sheet), cfg$samples$factors %||% list())
    }
    caps <- resolve_capabilities(cfg, samples, read_registry(), project_gene_lists(cfg))
    expect_identical(listed[[1]], caps$module[caps$status == "run"], label = configs[[i]])
  }
})

test_that("modules_doc() renders schema types, defaults and conditions", {
  expect_identical(doc_default(NULL), "null")
  expect_identical(doc_default(list()), "[]")
  expect_identical(doc_default(stats::setNames(list(), character())), "{}")
  expect_identical(doc_default(c(0.2, 0.3)), "[0.2,0.3]")
  expect_identical(doc_default(list(list(id = "chr_y"))), "[{\"id\":\"chr_y\"}]")
  expect_identical(schema_type(list(type = c("string", "null"))), "string or null")
  expect_identical(
    schema_type(list(type = "array", items = list(type = "string", enum = c("vst", "logcpm_tmm")))),
    "array of string (vst, logcpm_tmm)"
  )
  expect_identical(doc_condition(list(configured = c("profiles", "mask.include"))), "configured: profiles, mask.include")
  expect_identical(doc_condition(list(contrasts_defined = TRUE)), "contrasts_defined: true")
  expect_identical(doc_text("a |b\n c"), "a \\|b c")
  expect_true(all(c("## project", "## inputs", "## modules") %in% config_doc()))
  doc <- modules_doc()
  expect_true(all(sprintf(
    "## %02d %s", vapply(read_registry(), function(m) as.integer(m$number), 1L),
    names(read_registry())
  ) %in% doc))
})

test_that("the schema has no empty keys (e.g. a comma inside an unquoted flow-mapping description)", {
  empty_keys <- function(x, path = "") {
    if (!is.list(x)) {
      return(character())
    }
    nms <- names(x) %||% rep("", length(x))
    unlist(lapply(seq_along(x), function(i) {
      p <- paste0(path, "/", if (nzchar(nms[[i]])) nms[[i]] else i)
      if (is.null(x[[i]]) && nms[[i]] != "default") p else empty_keys(x[[i]], p)
    }))
  }
  expect_length(empty_keys(yaml::read_yaml(config_schema_path())), 0)
})
