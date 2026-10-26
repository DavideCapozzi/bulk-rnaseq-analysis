#' Read a project configuration
#'
#' Reads a project YAML, fills defaults from the JSON Schema shipped with the package and checks the
#' cross-field rules that a JSON Schema cannot express. Inside the workflow the config has already been
#' validated by `snakemake.utils.validate`; filling defaults again is idempotent, so the same function serves
#' interactive use.
#'
#' @param path `character(1)`. Path to the project YAML, or to `00_provenance/config_resolved.yaml` written by
#'   the workflow (which carries the original config folder in `_config_dir`).
#' @param schema `character(1)`. Path to the schema; defaults to the copy installed with the package.
#' @return A named list. Attribute `config_dir` holds the folder that relative paths are resolved against.
#' @examples
#' cfg <- tempfile(fileext = ".yaml")
#' yaml::write_yaml(list(project = list(name = "demo", output_dir = "out"), inputs = list()), cfg)
#' x <- read_config(cfg)
#' x$project$seed
#' @export
read_config <- function(path, schema = config_schema_path()) {
  assert_file(path, "config")
  cfg <- yaml::read_yaml(path)
  config_dir <- cfg[["_config_dir"]] %||% dirname(normalizePath(path))
  cfg[["_config_dir"]] <- NULL
  cfg <- apply_schema_defaults(cfg, yaml::read_yaml(schema))
  validate_config(cfg)
  attr(cfg, "config_dir") <- config_dir
  cfg
}

#' Path of the configuration schema
#'
#' @return `character(1)`, path to `config/schema.yaml` as installed in `inst/schema/`.
#' @examples
#' file.exists(config_schema_path())
#' @export
config_schema_path <- function() package_file("schema", "schema.yaml")

#' Fill defaults from a JSON Schema
#'
#' Walks the schema and sets every missing property that declares a `default`, recursing into objects and
#' into arrays of objects. This mirrors what `snakemake.utils.validate(set_default = TRUE)` does, so that
#' package functions behave the same with or without Snakemake.
#'
#' @param x `list`. Parsed config (or a sub-tree).
#' @param schema `list`. Parsed schema (or the matching sub-schema).
#' @param root `list`. Whole schema, used to resolve `$ref`.
#' @return `x` with defaults filled in.
#' @examples
#' s <- list(type = "object", properties = list(a = list(type = "number", default = 2)))
#' apply_schema_defaults(list(), s)
#' @export
apply_schema_defaults <- function(x, schema, root = schema) {
  schema <- deref_schema(schema, root)
  if (is.null(x)) x <- list()
  props <- schema$properties
  if (is.list(x) && !is.null(props)) {
    for (nm in names(props)) {
      sub <- deref_schema(props[[nm]], root)
      if (is.null(x[[nm]]) && !is.null(sub$default)) {
        x[[nm]] <- sub$default
      }
      if (!is.null(x[[nm]])) {
        if (identical(sub$type, "object")) {
          x[[nm]] <- apply_schema_defaults(x[[nm]], sub, root)
        } else if (identical(sub$type, "array") && !is.null(sub$items$properties)) {
          x[[nm]] <- lapply(x[[nm]], apply_schema_defaults, schema = sub$items, root = root)
        }
      }
    }
  }
  x
}

deref_schema <- function(s, root) {
  ref <- s[["$ref"]]
  if (is.null(ref)) {
    return(s)
  }
  parts <- strsplit(sub("^#/", "", ref), "/", fixed = TRUE)[[1]]
  target <- root
  for (p in parts) target <- target[[p]]
  utils::modifyList(target, s[setdiff(names(s), "$ref")])
}

#' Check cross-field rules of a config
#'
#' Rules that JSON Schema cannot express: sample-level inputs need a sample sheet; inputs whose adapter is not
#' implemented yet are refused instead of being ignored; contrast, pair and design identifiers are unique;
#' contrast pairs refer to defined contrasts.
#'
#' @param cfg `list`. Parsed config with defaults.
#' @return `TRUE` invisibly; otherwise an error naming the field.
#' @examples
#' validate_config(list(project = list(name = "x", output_dir = "o"), inputs = list()))
#' @export
validate_config <- function(cfg) {
  sample_level <- c("counts", "abundance", "normalized", "coverage", "qc")
  if (any(sample_level %in% names(cfg$inputs)) && is.null(cfg$samples$sheet)) {
    abort(
      "config: inputs ", paste(intersect(sample_level, names(cfg$inputs)), collapse = ", "),
      " need samples.sheet (sample metadata come only from the sample sheet)"
    )
  }
  pending <- intersect(c("qc", "coverage", "checksums"), names(cfg$inputs))
  if (length(pending) > 0) {
    abort(
      "config: inputs.", paste(pending, collapse = ", inputs."), " not implemented in bulkrnaseq ",
      utils::packageVersion("bulkrnaseq"), "; remove it from the config (it would be ignored)"
    )
  }
  for (t in cfg$inputs$contrast_tables %||% list()) {
    has_metric <- t$adapter %in% c("deseq2", "noiseq") || !is.null(t$columns$padj) || !is.null(t$columns$posterior_prob)
    if (has_metric && is.null(t$df_residual)) {
      abort(
        "config: inputs.contrast_tables '", t$id, "' needs df_residual, the residual degrees of freedom of the ",
        "provider's model (n samples - rank of the design; 0 without replicates): it decides whether its ",
        "p-values are sample-level"
      )
    }
  }
  check_unique <- function(ids, what) {
    dup <- unique(ids[duplicated(ids)])
    if (length(dup) > 0) abort("config: duplicated ", what, " id(s): ", paste(dup, collapse = ", "))
  }
  contrast_ids <- c(
    vapply(cfg$contrasts %||% list(), `[[`, "", "id"),
    vapply(cfg$inputs$contrast_tables %||% list(), `[[`, "", "id"),
    vapply(cfg$designs %||% list(), `[[`, "", "id")
  )
  check_unique(contrast_ids, "contrast/design")
  check_unique(vapply(cfg$contrast_pairs %||% list(), `[[`, "", "id"), "contrast pair")
  for (p in cfg$contrast_pairs %||% list()) {
    bad <- setdiff(c(p$first, p$second), contrast_ids)
    if (length(bad) > 0) {
      abort("config: contrast_pairs '", p$id, "' refers to undefined contrast(s): ", paste(bad, collapse = ", "))
    }
  }
  invisible(TRUE)
}

# Absolute path of a config entry, resolved against the config folder.
config_path <- function(cfg, path) resolve_path(path, attr(cfg, "config_dir") %||% getwd())

input_path <- function(cfg, path) {
  base <- config_path(cfg, cfg$inputs$base_dir %||% ".")
  resolve_path(path, base)
}

output_dir <- function(cfg) config_path(cfg, cfg$project$output_dir)

module_params <- function(cfg, module) cfg$modules[[module]] %||% list()

# One config section with the schema defaults filled in, e.g. with_defaults(params, c("modules", "pca")) or
# with_defaults(cfg$thresholds, "thresholds"). The schema is the only place where defaults live: inside the
# workflow the config is already complete, and module functions called directly get the same values.
with_defaults <- function(x, section) {
  root <- yaml::read_yaml(config_schema_path())
  sub <- root
  for (key in section) sub <- deref_schema(sub$properties[[key]], root)
  if (is.null(sub)) abort("with_defaults(): no schema section ", paste(section, collapse = "."))
  apply_schema_defaults(x %||% list(), sub, root)
}
