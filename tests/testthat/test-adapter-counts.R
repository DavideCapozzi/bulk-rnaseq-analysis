star_files <- function() {
  libs <- c("lib_A1", "lib_A2", "lib_B1", "lib_B2")
  stats::setNames(fixture("scenario_a", "counts", paste0(libs, "_ReadsPerGene.out.tab")), libs)
}

test_that("STAR 4-column files: chosen strand column, N_* rows moved to QC", {
  r <- read_star_counts(star_files(), strand_column = "reverse")
  raw <- utils::read.delim(star_files()[[1]], header = FALSE)
  genes <- raw[!startsWith(raw$V1, "N_"), ]
  expect_equal(unname(r$counts[, "lib_A1"]), genes$V4)
  expect_type(r$counts, "integer")
  expect_false(any(startsWith(rownames(r$counts), "N_")))
  expect_setequal(unique(r$qc$metric), c("n_unmapped", "n_multimapping", "n_nofeature", "n_ambiguous"))
  expect_equal(r$qc$value[r$qc$sample_id == "lib_A1" & r$qc$metric == "n_unmapped"], raw$V4[raw$V1 == "N_unmapped"])
  expect_equal(nrow(r$provenance), 4)
})

test_that("STAR strand setting must match the file layout", {
  expect_error(read_star_counts(c(s = fixture("formats", "star_2col.txt")), "reverse"), "strand_column: single")
  expect_error(read_star_counts(star_files()[1], "single"), "4 columns")
  two <- read_star_counts(c(s = fixture("formats", "star_2col.txt")), "single")
  expect_equal(nrow(two$counts), 20)
  expect_equal(nrow(two$qc), 0)
})

test_that("duplicated gene IDs follow on_duplicate", {
  f <- c(s = fixture("formats", "star_duplicated.txt"))
  expect_error(read_star_counts(f, "reverse"), "duplicated gene ID")
  summed <- read_star_counts(f, "reverse", on_duplicate = "sum")$counts
  first <- read_star_counts(f, "reverse", on_duplicate = "first")$counts
  expect_equal(nrow(summed), 5)
  expect_equal(summed[2, 1], 2L + 6L)
  expect_equal(first[2, 1], 2L)
})

test_that("StringTie: duplicated Gene ID handled explicitly, '-' names become NA", {
  f <- stats::setNames(
    fixture("scenario_a", "abundance", c("lib_A1_gene_abundance.tab", "lib_A2_gene_abundance.tab")),
    c("lib_A1", "lib_A2")
  )
  expect_error(read_stringtie(f), "duplicated Gene ID")
  s <- read_stringtie(f, on_duplicate = "sum")
  raw <- utils::read.delim(f[[1]], check.names = FALSE)
  id <- raw[["Gene ID"]][1]
  expect_equal(unname(s$tpm[id, "lib_A1"]), sum(raw$TPM[raw[["Gene ID"]] == id]))
  expect_equal(nrow(s$tpm), length(unique(raw[["Gene ID"]])))
  expect_true(anyNA(s$features$gene_name))
  expect_equal(s$name_mismatches, 0)
  first <- read_stringtie(f, on_duplicate = "first")
  expect_equal(unname(first$tpm[id, "lib_A1"]), raw$TPM[1])
})

test_that("matrix tables from Excel keep counts as integers and never accept decimals as counts", {
  r <- read_matrix_table(fixture("scenario_b", "raw_counts.xlsx"), "Geneid", "Sheet1", "counts")
  expect_type(r$values, "integer")
  expect_equal(ncol(r$values), 6)
  expect_error(
    read_matrix_table(fixture("scenario_b", "normalized_counts.xlsx"), "Gene_ID", "NormalizedCount", "counts"),
    "never converted back to counts"
  )
  n <- read_matrix_table(fixture("scenario_b", "normalized_counts.xlsx"), "Gene_ID", "NormalizedCount", "normalized")
  expect_type(n$values, "double")
  expect_error(read_matrix_table(fixture("scenario_b", "raw_counts.xlsx"), "Geneid", "Nope"), "sheet 'Nope' not found")
})
