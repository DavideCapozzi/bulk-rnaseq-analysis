test_that("schema defaults are filled recursively, including module parameters", {
  cfg <- make_cfg(list(project = list(name = "x", output_dir = "o"), inputs = list()))
  expect_equal(cfg$project$seed, 1)
  expect_true(cfg$modules$volcano$enabled)
  expect_equal(cfg$modules$concordance$thresholds, list(0.58, 1))
  expect_equal(cfg$gene_filter$min_cpm, 2)
})

test_that("module functions take their defaults from the schema, keeping the values given", {
  p <- with_defaults(list(top_loadings = 5), c("modules", "pca"))
  expect_equal(p$top_loadings, 5)
  expect_equal(p$ntop, list(500L, 1000L, 2000L, "all"))
  expect_equal(with_defaults(NULL, "thresholds")$padj, 0.05)
  expect_error(with_defaults(list(), c("modules", "nope")), "no schema section modules.nope")
})

test_that("read_config reads a YAML file and records its folder", {
  f <- tempfile(fileext = ".yaml")
  yaml::write_yaml(list(project = list(name = "demo", output_dir = "out"), inputs = list()), f)
  cfg <- read_config(f)
  expect_equal(attr(cfg, "config_dir"), dirname(normalizePath(f)))
})

test_that("cross-field rules give errors naming the field", {
  expect_error(validate_config(list(inputs = list(counts = list(adapter = "matrix")))), "samples.sheet")
  cfg <- list(contrasts = list(list(id = "a"), list(id = "a")), inputs = list())
  expect_error(validate_config(cfg), "duplicated contrast/design id")
  pairs <- list(list(id = "p", first = "a", second = "zz"))
  cfg <- list(contrasts = list(list(id = "a")), contrast_pairs = pairs, inputs = list())
  expect_error(validate_config(cfg), "undefined contrast\\(s\\): zz")
  cfg <- list(inputs = list(checksums = list(file = "md5sum.txt")))
  expect_error(validate_config(cfg), "inputs.checksums not implemented")
  tab <- list(id = "t", adapter = "deseq2", file = "x.xlsx")
  expect_error(validate_config(list(inputs = list(contrast_tables = list(tab)))), "'t' needs df_residual")
  tab$df_residual <- 0
  expect_true(validate_config(list(inputs = list(contrast_tables = list(tab)))))
})

test_that("contrast pairs use first/second because YAML 1.1 reads the key 'y' as a boolean", {
  expect_named(yaml::yaml.load("{id: a, x: C1, y: C2}"), c("id", "x", "TRUE"))
  pair <- yaml::yaml.load("{id: a, first: C1, second: C2}")
  expect_equal(pair$second, "C2")
})
