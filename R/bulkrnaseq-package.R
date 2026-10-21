#' bulkrnaseq: reproducible analysis of bulk RNA-seq provider outputs
#'
#' The package has four layers. Adapters (`read_*()`) turn one file format into a canonical object;
#' canonical objects (a `SummarizedExperiment`, a long contrast table, a long QC table) are checked by
#' `validate_*()`; analysis modules (`run_<module>()`) are pure functions of canonical objects and a parameter
#' list; writers serialise tables with a `.meta.json` sidecar. The Snakemake workflow in `workflow/` calls
#' these functions through thin wrappers.
#'
#' @keywords internal
#' @importFrom ggplot2 .data
"_PACKAGE"

# data.table is used for fast, exact TSV reading and writing only; package code does not rely on its
# non-standard evaluation, so no global variable declarations are needed.
.datatable.aware <- TRUE # nolint: object_name_linter.
