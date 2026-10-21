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
