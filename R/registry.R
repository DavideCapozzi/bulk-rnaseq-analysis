#' Read the module registry
#'
#' The registry (`inst/modules.yaml`) is the single place where a module declares the data levels it needs,
#' the levels it provides, its resources, conditions and outputs. Capability resolution and documentation are
#' generated from it.
#'
#' @param path `character(1)`. Registry file; defaults to the installed copy.
#' @return A named list of module entries; each entry gains a `name` field.
#' @examples
#' reg <- read_registry()
#' reg$volcano$requires
#' @export
read_registry <- function(path = package_file("modules.yaml")) {
  reg <- yaml::read_yaml(path)$modules
  for (nm in names(reg)) reg[[nm]]$name <- nm
  numbers <- vapply(reg, function(m) as.integer(m$number), 1L)
  if (anyDuplicated(numbers)) abort("modules.yaml: duplicated module numbers")
  reg[order(numbers)]
}

#' Organism information
#'
#' All organism-specific identifiers (OrgDb package, KEGG code, chromosome names, mitochondrial prefix) come
#' from `inst/organisms.yaml`; package code never hard-codes them.
#'
#' @param organism `character(1)`. Key of the organism table (only `"human"` is implemented).
#' @return A named list with fields `species`, `orgdb`, `kegg_code`, `msigdb_species`, `ensembl_gene_prefix`,
#'   `chromosome_x`, `chromosome_y`, `chromosome_mito`, `mito_symbol_prefix`.
#' @examples
#' organism_info("human")$chromosome_mito
#' @export
organism_info <- function(organism = "human") {
  tab <- yaml::read_yaml(package_file("organisms.yaml"))
  if (!organism %in% names(tab)) {
    abort("organism '", organism, "' is not in organisms.yaml (available: ", paste(names(tab), collapse = ", "), ")")
  }
  tab[[organism]]
}

# Folder name of a module output, e.g. "09_contrasts".
module_dir_name <- function(module, registry = read_registry()) {
  if (module %in% c("provenance", "ingest")) {
    return(c(provenance = "00_provenance", ingest = "01_ingest")[[module]])
  }
  m <- registry[[module]]
  if (is.null(m)) abort("unknown module '", module, "'")
  sprintf("%02d_%s", as.integer(m$number), module)
}

# TRUE when the module's function exists in the package (otherwise the module is skipped as "not implemented").
module_implemented <- function(m) exists(m[["function"]], envir = asNamespace("bulkrnaseq"), inherits = FALSE)

# Markdown reference of the modules (docs/MODULES.md), built from the registry and the config schema only.
# tools/make_modules_doc.R writes the file; a test fails when the committed file differs from this output.
modules_doc <- function(registry = read_registry(), schema = yaml::read_yaml(config_schema_path())) {
  implemented <- vapply(registry, module_implemented, TRUE)
  module_schema <- schema$properties$modules$properties
  code <- function(x) if (length(x) == 0) "none" else paste0("`", x, "`", collapse = ", ")
  status <- function(i) if (i) "implemented" else "not implemented"
  index <- vapply(registry, function(m) {
    sprintf(
      "| %02d | [`%s`](#%02d-%s) | %s | %s | %s |", as.integer(m$number), m$name, as.integer(m$number), m$name,
      status(implemented[[m$name]]), doc_levels(m, code), code(unlist(m[["provides"]]))
    )
  }, "")
  sections <- unlist(lapply(registry, function(m) {
    dir <- sprintf("%02d_%s", as.integer(m$number), m$name)
    params <- if (!is.null(module_schema[[m$name]])) schema_param_rows(module_schema[[m$name]], schema)
    c(
      "", sprintf("## %02d %s", as.integer(m$number), m$name), "",
      sprintf(
        "Status: %s (`%s`). Output folder: `%s/`.", status(implemented[[m$name]]), m[["function"]], dir
      ),
      "",
      paste0("- Data levels: ", doc_levels(m, code)),
      paste0("- Optional levels: ", code(unlist(m[["optional"]]))),
      paste0("- Resources: ", code(unlist(m[["resources"]]))),
      paste0("- Conditions: ", code(vapply(m[["conditions"]] %||% list(), doc_condition, ""))),
      paste0("- Provides: ", code(unlist(m[["provides"]]))),
      paste0("- Tables: ", code(unlist(m[["outputs"]]))),
      paste0("- Figures: ", code(unlist(m[["figures"]]))),
      if (length(m[["caveats"]])) c("", "Caveats:", "", paste0("- ", unlist(m[["caveats"]]))),
      "",
      if (is.null(params)) {
        "No parameters in the schema."
      } else {
        c(
          sprintf(
            "Parameters (`modules.%s`%s):", m$name,
            if (implemented[[m$name]]) "" else "; draft, not read by any code yet"
          ),
          "", "| Parameter | Type | Default | Description |", "|---|---|---|---|", params
        )
      }
    )
  }), use.names = FALSE)
  defs <- vapply(names(schema$`$defs`), function(nm) {
    paste0("- `", nm, "`: ", doc_text(schema$`$defs`[[nm]]$description %||% schema_type(schema$`$defs`[[nm]])))
  }, "")
  unname(c(
    "# Modules", "",
    paste(
      "<!-- Generated from inst/modules.yaml and inst/schema/schema.yaml by `Rscript tools/make_modules_doc.R`;",
      "do not edit by hand (a test fails when this file is out of date). -->"
    ),
    "",
    paste(
      "Every module declares the data levels it needs. A module whose needs are not met, or that is not",
      "implemented yet, is skipped with the reason in `00_provenance/capabilities.json` and `report/summary.md`.",
      "Defaults come from the config schema; `bulkrnaseq plan <config>` shows what a project runs."
    ),
    "", sprintf(
      "%d modules: %d implemented, %d not implemented.", length(registry), sum(implemented), sum(!implemented)
    ),
    "", "| NN | Module | Status | Data levels | Provides |", "|---|---|---|---|---|", index,
    "", "Schema types used below:", "", defs,
    sections
  ))
}

doc_levels <- function(m, code) {
  all_of <- unlist(m[["requires"]])
  any_of <- unlist(m[["requires_any"]])
  parts <- c(if (length(all_of)) code(all_of), if (length(any_of)) paste("one of", code(any_of)))
  if (length(parts)) paste(parts, collapse = "; ") else "none"
}

doc_condition <- function(cond) {
  value <- unlist(cond)
  if (is.logical(value)) value <- tolower(value)
  paste0(names(cond), ": ", paste(value, collapse = ", "))
}

doc_text <- function(x) gsub("|", "\\|", trimws(gsub("\\s+", " ", paste(x, collapse = " "))), fixed = TRUE)

schema_type <- function(p) {
  type <- paste(unlist(p$type), collapse = " or ")
  if (identical(p$type, "array") && !is.null(p$items)) type <- paste("array of", schema_type(p$items))
  if (!is.null(p$enum)) type <- paste0(type, " (", paste(unlist(p$enum), collapse = ", "), ")")
  if (!nzchar(type) && !is.null(p$oneOf)) type <- paste(vapply(p$oneOf, schema_type, ""), collapse = " or ")
  type
}

# One Markdown table row per parameter, recursing into objects (a.b) and arrays of objects (a[].b).
schema_param_rows <- function(s, root, prefix = "") {
  s <- deref_schema(s, root)
  unlist(lapply(names(s$properties), function(nm) {
    raw <- s$properties[[nm]]
    p <- deref_schema(raw, root)
    path <- paste0(prefix, nm)
    type <- if (!is.null(raw$`$ref`)) sub(".*/", "", raw$`$ref`) else schema_type(p)
    default <- if ("default" %in% names(p)) paste0("`", doc_default(p$default), "`") else ""
    # a $ref inherits the description of its definition, which the "Schema types" list already gives
    description <- if (!is.null(raw$`$ref`)) raw$description else p$description
    row <- sprintf("| `%s` | %s | %s | %s |", path, doc_text(type), doc_text(default), doc_text(description %||% ""))
    items <- if (identical(p$type, "array") && !is.null(p$items)) deref_schema(p$items, root)
    c(
      row,
      if (identical(p$type, "object") && !is.null(p$properties)) schema_param_rows(p, root, paste0(path, ".")),
      if (!is.null(items$properties)) schema_param_rows(items, root, paste0(path, "[]."))
    )
  }))
}

doc_default <- function(x) {
  if (is.null(x)) {
    return("null")
  }
  if (length(x) == 0) {
    return(if (is.null(names(x))) "[]" else "{}")
  }
  as.character(jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA))
}
