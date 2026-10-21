test_that("the run manifest records the toolkit identity checked by the workflow", {
  out <- withr::local_tempdir()
  prov <- file.path(out, "00_provenance")
  dir.create(prov)
  caps <- data.frame(module = "volcano", number = 14L, dir = "14_volcano", status = "run", reason = "")
  jsonlite::write_json(caps, file.path(prov, "capabilities.json"))
  toolkit <- list(
    version = "0.1.0", pinned_version = "0.1.0", git_commit = strrep("a", 40), git_describe = "v0.1.0",
    git_dirty = FALSE, code_md5 = strrep("b", 32), dev_mode = FALSE
  )
  jsonlite::write_json(toolkit, file.path(prov, "toolkit.json"), auto_unbox = TRUE)
  cfg <- structure(list(project = list(name = "demo", seed = 1), resources = list(cache_dir = out)), config_dir = out)
  m <- jsonlite::fromJSON(write_run_manifest(cfg, out, repo = "/nonexistent"))
  expect_equal(m$toolkit, toolkit)
})
