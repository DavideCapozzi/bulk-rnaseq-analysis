fixture <- function(...) testthat::test_path("..", "fixtures", ...)

# Minimal config object with config_dir set, as read_config() would return.
make_cfg <- function(x, dir = tempfile("cfg_")) {
  dir.create(dir, showWarnings = FALSE)
  x <- apply_schema_defaults(x, yaml::read_yaml(config_schema_path()))
  attr(x, "config_dir") <- dir
  x
}

toy_se <- function(m, rd = NULL, cd = NULL) {
  rd <- rd %||% S4Vectors::DataFrame(
    feature_id = rownames(m), gene_symbol = rownames(m), biotype = "protein_coding",
    chromosome = "1"
  )
  cd <- cd %||% S4Vectors::DataFrame(sample_id = colnames(m), row.names = colnames(m))
  SummarizedExperiment::SummarizedExperiment(list(counts = m), rowData = rd, colData = cd)
}
