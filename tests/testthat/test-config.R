test_that("schema defaults are filled recursively, including module parameters", {
  cfg <- make_cfg(list(project = list(name = "x", output_dir = "o"), inputs = list()))
  expect_equal(cfg$project$seed, 1)
  expect_true(cfg$modules$volcano$enabled)
  expect_equal(cfg$modules$concordance$thresholds, list(0.58, 1))
  expect_equal(cfg$gene_filter$min_cpm, 2)
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
})

test_that("contrast pairs use first/second because YAML 1.1 reads the key 'y' as a boolean", {
  expect_named(yaml::yaml.load("{id: a, x: C1, y: C2}"), c("id", "x", "TRUE"))
  pair <- yaml::yaml.load("{id: a, first: C1, second: C2}")
  expect_equal(pair$second, "C2")
})
