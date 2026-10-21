#' Write the run manifest
#'
#' Records what is needed to reproduce a run: md5 of the resolved config, toolkit identity (version pinned by
#' the config, git commit, dirty flag, code md5, development mode), versions of R and of the packages used, md5
#' of every input file, resources, and every module's status with its reason. The toolkit identity is the one
#' the workflow checked before the run (`00_provenance/toolkit.json`); without it, the git state of `repo` is
#' read now. Wall-clock times go to a separate `run_timestamps.json`, so the manifest is identical across runs
#' of the same inputs.
#'
#' @param cfg Config from [read_config()].
#' @param out `character(1)`. Project output folder.
#' @param repo `character(1)` or `NULL`. Toolkit repository, read only when `toolkit.json` is missing.
#' @return Path of `run_manifest.json`, invisibly.
#' @examples
#' \donttest{
#' # called by the workflow after all modules
#' }
#' @export
write_run_manifest <- function(cfg, out, repo = NULL) {
  dir <- file.path(out, "00_provenance")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  resolved <- file.path(dir, "config_resolved.yaml")
  caps <- read_capabilities(out)
  pkgs <- c("bulkrnaseq", sort(unique(c(desc_packages("Imports"), desc_packages("Suggests")))))
  version_of <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  versions <- vapply(pkgs, version_of, "")
  inputs <- file.path(out, "01_ingest", "input_files.tsv")
  toolkit_json <- file.path(dir, "toolkit.json")
  toolkit <- if (file.exists(toolkit_json)) {
    jsonlite::fromJSON(toolkit_json)
  } else if (!is.null(repo)) {
    git_state(repo)
  } else {
    list(git_commit = NA_character_, git_dirty = NA)
  }
  cache_manifest <- file.path(config_path(cfg, cfg$resources$cache_dir %||% "resources"), "resources_manifest.tsv")
  if (file.exists(cache_manifest)) file.copy(cache_manifest, file.path(dir, "resources_manifest.tsv"), overwrite = TRUE)
  manifest <- list(
    project = cfg$project$name,
    config_md5 = if (file.exists(resolved)) md5_file(resolved) else NA_character_,
    toolkit = toolkit,
    r_version = R.version.string,
    packages = as.list(versions[!is.na(versions)]),
    seed = cfg$project$seed,
    inputs = if (file.exists(inputs)) as.data.frame(data.table::fread(inputs, sep = "\t")) else list(),
    modules = caps[, c("module", "number", "dir", "status", "reason")]
  )
  path <- file.path(dir, "run_manifest.json")
  write_json(manifest, path)
  stamp <- list(provenance_written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  write_json(stamp, file.path(dir, "run_timestamps.json"))
  invisible(path)
}

desc_packages <- function(field) {
  d <- utils::packageDescription("bulkrnaseq")[[field]]
  if (is.null(d)) {
    return(character())
  }
  p <- trimws(sub("\\(.*\\)", "", strsplit(d, ",")[[1]]))
  setdiff(p[nzchar(p)], c("R", "methods", "stats", "tools", "utils"))
}

git_state <- function(repo) {
  commit <- tryCatch(system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_, warning = function(w) NA_character_
  )
  status <- tryCatch(system2("git", c("-C", shQuote(repo), "status", "--porcelain"), stdout = TRUE, stderr = FALSE),
    error = function(e) character(), warning = function(w) character()
  )
  list(git_commit = if (length(commit)) commit[[1]] else NA_character_, git_dirty = length(status) > 0)
}

#' Write the run summary
#'
#' @param cfg Config from [read_config()].
#' @param out `character(1)`. Project output folder.
#' @return Path of `report/summary.md`, invisibly.
#' @examples
#' \donttest{
#' # called by the workflow as the last step
#' }
#' @export
write_report <- function(cfg, out) {
  caps <- read_capabilities(out)
  lines <- c(
    paste0("# Run summary: ", cfg$project$name), "", "| Module | Status | Reason |", "|---|---|---|",
    sprintf("| %s | %s | %s |", caps$dir, caps$status, caps$reason), "", "## Tables", ""
  )
  for (d in caps$dir[caps$status == "run"]) {
    tabs <- sort(list.files(file.path(out, d, "tables"), pattern = "\\.tsv(\\.gz)?$"))
    lines <- c(lines, paste0("- ", d, ": ", if (length(tabs)) paste(tabs, collapse = ", ") else "(none)"))
  }
  dir.create(file.path(out, "report"), showWarnings = FALSE)
  path <- file.path(out, "report", "summary.md")
  writeLines(lines, path)
  invisible(path)
}
