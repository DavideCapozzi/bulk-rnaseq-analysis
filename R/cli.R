#' Run one workflow step
#'
#' Single entry point for the Snakemake wrappers, so that workflow scripts contain no logic. Steps:
#' `capabilities`, `fetch_resources`, `ingest`, `module` (with `module` set), `provenance`, `report`.
#'
#' @param step `character(1)`. Step name.
#' @param config `character(1)`. Path to `00_provenance/config_resolved.yaml` (or a project config).
#' @param module `character(1)` or `NULL`. Module name for `step = "module"`.
#' @return Paths written, invisibly.
#' @examples
#' \donttest{
#' # run_step_cli("capabilities", "results/00_provenance/config_resolved.yaml")
#' }
#' @export
run_step_cli <- function(step, config, module = NULL) {
  cfg <- read_config(config)
  out <- output_dir(cfg)
  switch(step,
    capabilities = write_capabilities(cfg, out),
    fetch_resources = fetch_resources(cfg),
    ingest = {
      written <- write_ingest_outputs(run_ingest(cfg), out)
      done_file(file.path(out, "01_ingest"), written)
    },
    module = run_module_step(cfg, out, module),
    provenance = write_run_manifest(cfg, out),
    report = write_report(cfg, out),
    abort("run_step_cli(): unknown step '", step, "'")
  )
}

project_gene_lists <- function(cfg) {
  read_gene_lists(vapply(unlist(cfg$resources$gene_lists %||% list()), function(f) config_path(cfg, f), ""))
}

#' Write the capability decision of a run
#'
#' @param cfg Config from [read_config()].
#' @param out `character(1)`. Project output folder.
#' @return Path of `00_provenance/capabilities.json`, invisibly.
#' @examples
#' cfg <- list(
#'   project = list(name = "x", output_dir = tempfile()), inputs = list(), resources = list()
#' )
#' write_capabilities(cfg, tempfile())
#' @export
write_capabilities <- function(cfg, out) {
  samples <- if (!is.null(cfg$samples$sheet)) {
    read_sample_sheet(config_path(cfg, cfg$samples$sheet), cfg$samples$factors %||% list())
  }
  missing <- Filter(Negate(file.exists), configured_input_files(cfg, samples))
  if (length(missing) > 0) {
    abort(
      length(missing), " input file(s) not found (paths from inputs.base_dir and the sample sheet): ",
      paste(utils::head(missing, 5), collapse = ", "), if (length(missing) > 5) ", ..."
    )
  }
  caps <- resolve_capabilities(cfg, samples, read_registry(), project_gene_lists(cfg))
  dir <- file.path(out, "00_provenance")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, "capabilities.json")
  write_json(caps, path)
  invisible(path)
}

read_capabilities <- function(out) {
  as.data.frame(jsonlite::fromJSON(file.path(out, "00_provenance", "capabilities.json")))
}

done_file <- function(dir, written) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  rel <- sort(sub(paste0("^", normalizePath(dir), "/?"), "", normalizePath(written)))
  writeLines(rel, file.path(dir, "step.done"))
  invisible(file.path(dir, "step.done"))
}

run_module_step <- function(cfg, out, module) {
  registry <- read_registry()
  caps <- read_capabilities(out)
  row <- caps[caps$module == module, , drop = FALSE]
  if (nrow(row) != 1 || row$status != "run") abort("module '", module, "' is not scheduled to run: ", row$reason)
  m <- registry[[module]]
  dir <- file.path(out, row$dir)
  inputs <- module_inputs(cfg, out, caps, strsplit(row$depends_on, ",")[[1]])
  resource_key <- function(r) if (r == "annotation") cfg$annotation$gtf %||% "gtf" else r
  needed <- vapply(unlist(m[["resources"]]), resource_key, "")
  needed <- needed[needed %in% names(cfg$resources$items %||% list())]
  context <- list(
    cfg = cfg, organism = cfg$project$organism %||% "human", gene_lists = project_gene_lists(cfg),
    resources = as.list(verify_resources(cfg, needed)), seed = cfg$project$seed %||% 1
  )
  fn <- getExportedValue("bulkrnaseq", m$`function`)
  params <- module_params(cfg, module)
  result <- fn(inputs, params, context)
  xlsx <- if (isTRUE(cfg$export$xlsx)) file.path(out, "report", paste0(row$dir, ".xlsx"))
  written <- write_module_outputs(result, dir, params,
    xlsx_path = xlsx, gzip_tables = "contrast_stats",
    figure_sizes = result$figure_sizes %||% list()
  )
  if ("contrast_stats" %in% unlist(m[["provides"]])) {
    dir.create(file.path(dir, "objects"), showWarnings = FALSE)
    saveRDS(result$tables$contrast_stats, file.path(dir, "objects", "contrast_stats.rds"), version = 3)
    written <- c(written, file.path(dir, "objects", "contrast_stats.rds"))
  }
  done_file(dir, written)
}

module_inputs <- function(cfg, out, caps, depends_on) {
  ing <- file.path(out, "01_ingest", "objects")
  inputs <- list()
  if (file.exists(file.path(ing, "se.rds"))) {
    inputs$se <- readRDS(file.path(ing, "se.rds"))
    inputs$features <- as.data.frame(SummarizedExperiment::rowData(inputs$se))
  }
  if (file.exists(file.path(ing, "qc.rds"))) inputs$qc <- readRDS(file.path(ing, "qc.rds"))
  pieces <- list()
  if (file.exists(file.path(ing, "contrast_stats.rds"))) pieces$ingest <- readRDS(file.path(ing, "contrast_stats.rds"))
  for (dep in setdiff(depends_on, "ingest")) {
    f <- file.path(out, caps$dir[caps$module == dep], "objects", "contrast_stats.rds")
    if (file.exists(f)) pieces[[dep]] <- readRDS(f)
  }
  inputs$contrast_stats <- bind_contrast_tables(pieces)
  inputs
}
