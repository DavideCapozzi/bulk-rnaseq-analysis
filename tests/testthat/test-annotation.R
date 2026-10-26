test_that("GTF genes are read by streaming", {
  g <- read_gtf_genes(fixture("resources", "annotation.gtf.gz"))
  expect_named(g, c("ensembl_gene_id", "gene_symbol", "biotype", "chromosome", "start", "end", "strand"))
  expect_equal(sum(is.na(g$gene_symbol)), 1)
  expect_type(g$chromosome, "character")
  expect_true(all(c("MT", "X", "Y") %in% g$chromosome))
})

test_that("build_se maps input columns to sample_id and keeps provider-normalised values as NA where absent", {
  ss <- read_sample_sheet(fixture("scenario_b", "samples.tsv"))
  counts <- read_matrix_table(fixture("scenario_b", "raw_counts.xlsx"), "Geneid", "Sheet1", "counts")$values
  f <- fixture("scenario_b", "normalized_counts.xlsx")
  norm <- read_matrix_table(f, "Gene_ID", "NormalizedCount", "normalized")$values
  ann <- read_gtf_genes(fixture("resources", "annotation.gtf.gz"))
  b <- build_se(ss, counts = counts, normalized = norm, annotation = ann, id_type = "gene_symbol")
  se <- b$se
  expect_equal(colnames(se), ss$sample_id)
  expect_equal(b$report$normalized_genes, nrow(norm))
  expect_equal(sum(is.na(SummarizedExperiment::assay(se, "normalized")[, 1])), nrow(counts) - nrow(norm))
  expect_gt(b$report$symbols_mapped, 200)
  bad <- ss
  bad$input_column[1] <- "missing_column"
  expect_error(build_se(bad, counts = counts), "input_column value\\(s\\) not found.*missing_column")
})

test_that("symbols map to Ensembl only when unique in the annotation", {
  ann <- data.frame(
    ensembl_gene_id = c("E1", "E2", "E3"), gene_symbol = c("A", "B", "B"), biotype = "pc",
    chromosome = "1", start = 1L, end = 2L, strand = "+"
  )
  rd <- annotate_features(c("A", "B", "C"), ann, "gene_symbol")
  expect_equal(rd$ensembl_gene_id, c("E1", NA, NA))
  expect_equal(attr(rd, "report")$symbols_ambiguous, 1)
  expect_equal(attr(rd, "report")$symbols_unmapped, 1)
})
