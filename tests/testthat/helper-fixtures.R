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

# The tables and figures a module returned are declared in inst/modules.yaml (figures may add a _<suffix>).
expect_registry_outputs <- function(res, module) {
  m <- read_registry()[[module]]
  testthat::expect_true(all(names(res$tables) %in% unlist(m$outputs)), info = paste(module, "tables"))
  declared <- unlist(m$figures)
  known <- vapply(names(res$figures), function(f) f %in% declared || any(startsWith(f, paste0(declared, "_"))), TRUE)
  testthat::expect_true(all(known), info = paste(module, "figures"))
}

# Negative-binomial counts with a 10% block of genes changed by `effect` in group "b".
nb_se <- function(n_genes = 1000, groups = c("a", "a", "b", "b"), seed = 1, effect = 3) {
  set.seed(seed)
  n <- length(groups)
  mu <- matrix(rep(stats::rexp(n_genes, 1 / 200), n), n_genes)
  de <- seq_len(n_genes) <= n_genes / 10
  mu[de, groups == "b"] <- mu[de, groups == "b"] * effect
  m <- matrix(stats::rnbinom(n_genes * n, mu = mu, size = 20), n_genes,
    dimnames = list(sprintf("g%04d", seq_len(n_genes)), paste0("s", seq_len(n)))
  )
  storage.mode(m) <- "integer"
  rd <- S4Vectors::DataFrame(
    feature_id = rownames(m), gene_symbol = rownames(m), biotype = rep(c("protein_coding", "lncRNA"), length.out = n_genes),
    chromosome = "1", start = 1L, end = as.integer(seq(1000, 2e5, length.out = n_genes))
  )
  cd <- S4Vectors::DataFrame(sample_id = colnames(m), g = groups, row.names = colnames(m))
  SummarizedExperiment::SummarizedExperiment(list(counts = m), rowData = rd, colData = cd)
}
