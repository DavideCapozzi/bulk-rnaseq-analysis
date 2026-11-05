# Generated reference documentation: docs/MODULES.md (registry + module parameters) and docs/CONFIG.md (every
# other section of the config schema). tools/make_docs.R writes both files; a test fails when a committed file
# differs from this output, and another when a schema parameter has no description.

# Markdown reference of the modules (docs/MODULES.md), built from the registry and the config schema only.
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
          sprintf("Parameters (`modules.%s`):", m$name),
          "", "| Parameter | Type | Default | Description |", "|---|---|---|---|", params
        )
      }
    )
  }), use.names = FALSE)
  unname(c(
    "# Modules", "",
    doc_generated_note(),
    "",
    paste(
      "Every module declares the data levels it needs. A module whose needs are not met, or that is not",
      "implemented yet, is skipped with the reason in `00_provenance/capabilities.json` and `report/summary.md`.",
      "Defaults come from the config schema; `bulkrnaseq plan <config>` shows what a project runs.",
      "The other sections of the config are in [CONFIG.md](CONFIG.md)."
    ),
    "", sprintf(
      "%d modules: %d implemented, %d not implemented.", length(registry), sum(implemented), sum(!implemented)
    ),
    "", "| NN | Module | Status | Data levels | Provides |", "|---|---|---|---|---|", index,
    "", "Schema types used below:", "", doc_schema_types(schema),
    sections
  ))
}

# Markdown reference of the config sections other than `modules` (docs/CONFIG.md), from the schema only.
config_doc <- function(schema = yaml::read_yaml(config_schema_path())) {
  required <- unlist(schema$required)
  sections <- unlist(lapply(setdiff(names(schema$properties), "modules"), function(nm) {
    s <- deref_schema(schema$properties[[nm]], schema)
    # the first row is the section itself, already described by the paragraph above the table
    rows <- schema_param_rows(list(properties = stats::setNames(list(s), nm)), schema)[-1]
    c(
      "", paste0("## ", nm), "",
      paste0(doc_text(s$description %||% ""), if (nm %in% required) " Required." else ""),
      if (length(s$required)) {
        c("", paste0("Required fields: ", paste0("`", unlist(s$required), "`", collapse = ", "), "."))
      },
      "", "| Parameter | Type | Default | Description |", "|---|---|---|---|", rows
    )
  }), use.names = FALSE)
  modules <- names(schema$properties$modules$properties)
  numbers <- vapply(read_registry()[modules], function(m) as.integer(m$number), 1L)
  unname(c(
    "# Configuration", "",
    doc_generated_note(),
    "",
    paste(
      "A project is one YAML file validated against `inst/schema/schema.yaml` (`config/schema.yaml` links to it).",
      "Defaults are filled from the schema and recorded in `00_provenance/config_resolved.yaml`; relative paths",
      "resolve against the folder of the config file. `bulkrnaseq init` writes a skeleton and",
      "`config/examples/` holds one complete config per input scenario."
    ),
    "",
    paste0(
      "Sections: ", paste0("[`", names(schema$properties), "`](#", names(schema$properties), ")", collapse = ", "),
      ". Parameters of the analysis modules are in [MODULES.md](MODULES.md)."
    ),
    "", "Schema types used below:", "", doc_schema_types(schema),
    sections,
    "", "## modules", "",
    paste0(
      "Parameters of the analysis modules: one section per module, each with `enabled` (default `true`), ",
      "documented in [MODULES.md](MODULES.md): ",
      paste0("[`", modules, "`](MODULES.md#", sprintf("%02d", numbers), "-", modules, ")", collapse = ", "), "."
    )
  ))
}

doc_generated_note <- function() {
  paste(
    "<!-- Generated from inst/modules.yaml and inst/schema/schema.yaml by `Rscript tools/make_docs.R`;",
    "do not edit by hand (a test fails when this file is out of date). -->"
  )
}

doc_schema_types <- function(schema) {
  vapply(names(schema$`$defs`), function(nm) {
    paste0("- `", nm, "`: ", doc_text(schema$`$defs`[[nm]]$description %||% schema_type(schema$`$defs`[[nm]])))
  }, "", USE.NAMES = FALSE)
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

# One Markdown table row per parameter, recursing into objects (a.b, also when the type allows null), arrays of
# objects (a[].b) and maps whose values are objects (a.<key>.b).
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
    items <- if ("array" %in% unlist(p$type) && is.list(p$items)) deref_schema(p$items, root)
    values <- if (is.list(p$additionalProperties)) deref_schema(p$additionalProperties, root)
    c(
      row,
      if ("object" %in% unlist(p$type) && !is.null(p$properties)) schema_param_rows(p, root, paste0(path, ".")),
      if (!is.null(items$properties)) schema_param_rows(items, root, paste0(path, "[].")),
      if (!is.null(values$properties)) schema_param_rows(values, root, paste0(path, ".<key>."))
    )
  }))
}

# Parameters (rows of schema_param_rows() over the whole schema) whose description is empty.
schema_params_without_description <- function(schema = yaml::read_yaml(config_schema_path())) {
  rows <- schema_param_rows(schema, schema)
  sub("^\\| `([^`]+)`.*$", "\\1", rows[endsWith(rows, "|  |")])
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
