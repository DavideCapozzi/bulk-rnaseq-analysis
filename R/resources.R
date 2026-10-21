#' Fetch external resources into the project cache
#'
#' The only function that downloads or copies external resources. Each item of `resources.items` is retrieved
#' into `resources.cache_dir`, checked against its declared md5 (if any), and recorded in
#' `resources_manifest.tsv` with md5, size, version, retrieval date and source. Analysis steps never download:
#' they read the cache through [verify_resources()].
#'
#' Sources: `file://<path>` (copy; relative paths resolve against the config folder), `http(s)://...`
#' (download), `kegg_rest` (KEGG pathway-gene table for the organism, via KEGGREST).
#'
#' @param cfg Config from [read_config()].
#' @param items `character` or `NULL`. Subset of resource keys to fetch (default all).
#' @return Path of the manifest, invisibly.
#' @examples
#' src <- tempfile(fileext = ".tsv")
#' writeLines("a\tb", src)
#' items <- list(demo = list(source = paste0("file://", src)))
#' cfg <- list(resources = list(cache_dir = tempfile(), items = items))
#' attr(cfg, "config_dir") <- tempdir()
#' fetch_resources(cfg)
#' @export
fetch_resources <- function(cfg, items = NULL) {
  cache <- config_path(cfg, cfg$resources$cache_dir %||% "resources")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(cache, "resources_manifest.tsv")
  manifest <- if (file.exists(manifest_path)) read_resources_manifest(manifest_path) else NULL
  todo <- cfg$resources$items %||% list()
  if (!is.null(items)) todo <- todo[items]
  for (key in names(todo)) {
    it <- todo[[key]]
    src <- it$source
    if (src == "kegg_rest") {
      dest <- file.path(cache, paste0(key, ".tsv.gz"))
      version <- fetch_kegg(cfg$project$organism %||% "human", dest)
    } else if (startsWith(src, "file://")) {
      from <- config_path(cfg, sub("^file://", "", src))
      assert_file(from, paste0("resource '", key, "'"))
      dest <- file.path(cache, basename(from))
      file.copy(from, dest, overwrite = TRUE)
      version <- it$version %||% NA_character_
    } else if (grepl("^https?://", src)) {
      dest <- file.path(cache, basename(src))
      utils::download.file(src, dest, mode = "wb", quiet = TRUE)
      version <- it$version %||% NA_character_
    } else {
      abort("resource '", key, "': unsupported source '", src, "'")
    }
    md5 <- md5_file(dest)
    if (!is.null(it$md5) && !identical(md5, it$md5)) {
      abort("resource '", key, "': md5 ", md5, " differs from the declared ", it$md5, " (", src, ")")
    }
    row <- data.frame(
      resource = key, file = basename(dest), md5 = md5, bytes = file.size(dest),
      version = version, kind = it$kind %||% NA_character_,
      retrieved = format(Sys.Date()), source = src, stringsAsFactors = FALSE
    )
    manifest <- rbind(manifest[manifest$resource != key, , drop = FALSE], row)
  }
  manifest <- manifest[order(manifest$resource), , drop = FALSE]
  data.table::fwrite(manifest, manifest_path, sep = "\t")
  invisible(manifest_path)
}

#' Read a resources manifest
#'
#' @param path `character(1)`. `resources_manifest.tsv`.
#' @return `data.frame` with `resource`, `file`, `md5`, `bytes`, `version`, `kind`, `retrieved`, `source`.
#' @examples
#' f <- tempfile()
#' writeLines(c(
#'   "resource\tfile\tmd5\tbytes\tversion\tkind\tretrieved\tsource",
#'   "x\tx.tsv\tabc\t1\tNA\ttsv\t2026-01-01\tfile://x.tsv"
#' ), f)
#' read_resources_manifest(f)
#' @export
read_resources_manifest <- function(path) {
  as.data.frame(data.table::fread(path, sep = "\t", colClasses = "character", na.strings = "NA"))
}

#' Verify cached resources before an analysis step
#'
#' Fails if a needed resource is missing from the cache manifest or its md5 differs, so that analyses never run
#' on a changed or partially downloaded file.
#'
#' @param cfg Config from [read_config()].
#' @param needed `character`. Resource keys the step needs (default: all configured).
#' @return Named `character` of absolute file paths in the cache.
#' @examples
#' src <- tempfile(fileext = ".tsv")
#' writeLines("a\tb", src)
#' items <- list(demo = list(source = paste0("file://", src)))
#' cfg <- list(resources = list(cache_dir = tempfile(), items = items))
#' attr(cfg, "config_dir") <- tempdir()
#' fetch_resources(cfg)
#' verify_resources(cfg)
#' @export
verify_resources <- function(cfg, needed = names(cfg$resources$items)) {
  if (length(needed) == 0) {
    return(stats::setNames(character(), character()))
  }
  cache <- config_path(cfg, cfg$resources$cache_dir %||% "resources")
  manifest_path <- file.path(cache, "resources_manifest.tsv")
  hint <- "; run the explicit fetch step: snakemake fetch_resources (or bulkrnaseq::fetch_resources())"
  if (!file.exists(manifest_path)) abort("resource cache has no manifest (", manifest_path, ")", hint)
  man <- read_resources_manifest(manifest_path)
  out <- character()
  for (key in needed) {
    row <- man[man$resource == key, , drop = FALSE]
    if (nrow(row) != 1) abort("resource '", key, "' is not in the cache manifest", hint)
    path <- file.path(cache, row$file)
    if (!file.exists(path)) abort("resource '", key, "': cached file missing (", path, ")", hint)
    if (!identical(md5_file(path), row$md5)) {
      abort("resource '", key, "': cached file changed since it was fetched (md5 mismatch)")
    }
    out[[key]] <- path
  }
  out
}

fetch_kegg <- function(organism, dest) {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) abort("KEGG resource needs the 'KEGGREST' package")
  code <- organism_info(organism)$kegg_code
  links <- KEGGREST::keggLink("pathway", code)
  names_vec <- KEGGREST::keggList("pathway", code)
  tab <- data.frame(
    pathway_id = sub("^path:", "", unname(links)), entrez_id = sub("^[a-z]+:", "", names(links)),
    stringsAsFactors = FALSE
  )
  pn <- data.frame(pathway_id = sub("^path:", "", names(names_vec)), pathway_name = unname(names_vec))
  tab$pathway_name <- pn$pathway_name[match(tab$pathway_id, pn$pathway_id)]
  tab <- tab[order(tab$pathway_id, as.integer(tab$entrez_id)), c("pathway_id", "pathway_name", "entrez_id")]
  tmp <- tempfile(fileext = ".tsv")
  data.table::fwrite(tab, tmp, sep = "\t")
  gzip_deterministic(tmp, dest)
  # KEGG publishes no release number for the REST API; the "Last update" date identifies the content.
  info <- KEGGREST::keggInfo("pathway")
  updated <- regmatches(info, regexpr("Last update [0-9/]+", info))
  paste0(
    "KEGG pathway REST, ", if (length(updated)) tolower(updated) else "last update unknown", ", ",
    length(unique(tab$pathway_id)), " pathways"
  )
}
