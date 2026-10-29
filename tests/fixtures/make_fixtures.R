#!/usr/bin/env Rscript
# Generates every synthetic fixture used by the tests and the quickstart (seeded, deterministic).
# Run from the repository root, after python3 -I tests/fixtures/make_multiqc_fixture.py:
#   Rscript tests/fixtures/make_fixtures.R
# The files mimic the structure of real provider outputs (STAR N_* rows, StringTie duplicated gene IDs, Excel
# workbooks with several sheets, a sheet name containing ">", a wrong used range <dimension ref="A1"/>, a signed
# FoldChange column, a NOISeq table) but contain no real measurements. Gene symbols are real human symbols so that
# identifier mapping and ORA can be exercised; pathways and gene sets are synthetic.
suppressMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, quiet = TRUE))
root <- "tests/fixtures"
set.seed(1)

## genes -------------------------------------------------------------------------------------------------
eg <- org.Hs.eg.db::org.Hs.eg.db
entrez <- sort(as.integer(AnnotationDbi::keys(eg, "ENTREZID")))[1:4000]
tab <- suppressWarnings(suppressMessages(AnnotationDbi::select(eg, as.character(entrez), c("SYMBOL", "ENSEMBL", "CHR", "GENETYPE"), "ENTREZID")))
tab <- tab[!is.na(tab$ENSEMBL) & !is.na(tab$CHR) & tab$CHR %in% as.character(1:22) & tab$GENETYPE == "protein-coding", ]
tab <- tab[!duplicated(tab$ENTREZID) & !duplicated(tab$ENSEMBL) & !duplicated(tab$SYMBOL), ]
tab <- tab[seq_len(240), ]
special <- data.frame(
  SYMBOL = c("XIST", "TSIX", "RPS4Y1", "DDX3Y", "MT-CO1", "MT-ND1", "RNA5-8SN1"),
  CHR = c("X", "X", "Y", "Y", "MT", "MT", "21"),
  biotype = c("lncRNA", "lncRNA", "protein_coding", "protein_coding", "protein_coding", "protein_coding", "rRNA")
)
special$ENSEMBL <- sprintf("ENSG9%010d", seq_len(nrow(special)))
genes <- data.frame(
  ensembl = c(tab$ENSEMBL, special$ENSEMBL, sprintf("ENSG8%010d", 1:12)),
  symbol = c(tab$SYMBOL, special$SYMBOL, sprintf("LINC9%04d", 1:12)),
  chr = c(tab$CHR, special$CHR, rep(c("1", "2", "3"), 4)),
  biotype = c(rep("protein_coding", nrow(tab)), special$biotype, rep("lncRNA", 12)),
  stringsAsFactors = FALSE
)
genes$symbol[nrow(genes)] <- NA # one gene without a name in the GTF
genes$start <- as.integer(seq(10000, by = 25000, length.out = nrow(genes)))
genes$end <- genes$start + as.integer(sample(800:6000, nrow(genes), replace = TRUE))
genes$strand <- sample(c("+", "-"), nrow(genes), replace = TRUE)
genes$length <- genes$end - genes$start + 1L
n <- nrow(genes)

## resources ---------------------------------------------------------------------------------------------
res <- file.path(root, "resources")
dir.create(res, recursive = TRUE, showWarnings = FALSE)
attr9 <- function(i, extra = "") {
  nm <- if (is.na(genes$symbol[i])) "" else sprintf(' gene_name "%s";', genes$symbol[i])
  sprintf('gene_id "%s";%s gene_biotype "%s";%s', genes$ensembl[i], nm, genes$biotype[i], extra)
}
gtf <- unlist(lapply(seq_len(n), function(i) {
  tx <- sprintf(' transcript_id "%s"; transcript_biotype "%s";', sub("ENSG", "ENST", genes$ensembl[i]), genes$biotype[i])
  mid <- genes$start[i] + (genes$end[i] - genes$start[i]) %/% 2
  c(
    paste(genes$chr[i], "synthetic", "gene", genes$start[i], genes$end[i], ".", genes$strand[i], ".", attr9(i), sep = "\t"),
    paste(genes$chr[i], "synthetic", "exon", genes$start[i], mid - 200, ".", genes$strand[i], ".", attr9(i, tx), sep = "\t"),
    paste(genes$chr[i], "synthetic", "exon", mid, genes$end[i], ".", genes$strand[i], ".", attr9(i, tx), sep = "\t")
  )
}))
tmp <- tempfile(fileext = ".gtf")
writeLines(c("#!synthetic annotation for tests", gtf), tmp)
gzip_deterministic(tmp, file.path(res, "annotation.gtf.gz"))

ent <- tab$ENTREZID
kegg <- data.frame(
  pathway_id = rep(sprintf("hsa990%02d", 1:6), each = 20),
  pathway_name = rep(sprintf("Synthetic pathway %d", 1:6), each = 20),
  entrez_id = ent[1:120], stringsAsFactors = FALSE
)
data.table::fwrite(kegg, tmp <- tempfile(fileext = ".tsv"), sep = "\t")
gzip_deterministic(tmp, file.path(res, "kegg_pathways.tsv.gz"))

sets <- data.frame(
  gs_collection = "H", gs_subcollection = "", gs_name = rep(sprintf("SYNTHETIC_SET_%02d", 1:8), each = 25),
  ensembl_gene = genes$ensembl[1:200], gene_symbol = genes$symbol[1:200], db_version = "synthetic"
)
data.table::fwrite(sets, tmp <- tempfile(fileext = ".tsv"), sep = "\t")
gzip_deterministic(tmp, file.path(res, "gene_sets.tsv.gz"))

## counts ------------------------------------------------------------------------------------------------
simulate <- function(n_samples, depth, effect_rows, effect_cols, lfc) {
  mu <- matrix(rgamma(n, shape = 0.6, rate = 0.01), n, n_samples)
  mu[effect_rows, effect_cols] <- mu[effect_rows, effect_cols] * 2^lfc
  mu <- sweep(mu, 2, depth / colSums(mu), "*")
  matrix(rnbinom(length(mu), mu = mu, size = 8), n, n_samples)
}
up <- 1:15
down <- 16:30

## scenario a: STAR + StringTie ---------------------------------------------------------------------------
a <- file.path(root, "scenario_a")
for (d in c("counts", "abundance")) dir.create(file.path(a, d), recursive = TRUE, showWarnings = FALSE)
libs <- c("lib_A1", "lib_A2", "lib_B1", "lib_B2") # A = treated, B = control
cnt_a <- simulate(4, c(4e5, 3e5, 5e5, 3.5e5), c(up, down), 1:2, rep(c(2, -2), each = 15))
for (j in seq_along(libs)) {
  rev <- cnt_a[, j]
  fwd <- rbinom(n, rev, 0.02)
  uns <- rev + fwd + rbinom(n, rev, 0.05)
  nrows <- c(N_unmapped = 900 + j, N_multimapping = 500 + j, N_noFeature = 300 + j, N_ambiguous = 40 + j)
  lines <- c(
    sprintf("%s\t%d\t%d\t%d", names(nrows), nrows, nrows + 7L, nrows + 11L),
    sprintf("%s\t%d\t%d\t%d", genes$ensembl, uns, fwd, rev)
  )
  writeLines(lines, file.path(a, "counts", paste0(libs[j], "_ReadsPerGene.out.tab")))
  tpm <- (rev / genes$length) / sum(rev / genes$length) * 1e6
  fpkm <- (rev / genes$length * 1e3) / (sum(rev) / 1e6)
  st <- data.frame(
    `Gene ID` = genes$ensembl, `Gene Name` = ifelse(is.na(genes$symbol), "-", genes$symbol), Reference = genes$chr,
    Strand = genes$strand, Start = genes$start, End = genes$end, Coverage = round(rev / 10, 6),
    FPKM = round(fpkm, 6), TPM = round(tpm, 6), check.names = FALSE
  )
  # StringTie may report a gene twice (two loci); the second row is split off the first gene's values.
  dup <- st[1, ]
  dup$Start <- dup$Start + 5000L
  dup$FPKM <- 0.5
  dup$TPM <- 0.5
  st <- rbind(st[1, ], dup, st[-1, ])
  data.table::fwrite(st, file.path(a, "abundance", paste0(libs[j], "_gene_abundance.tab")), sep = "\t")
}
# s5: a library that failed before alignment, present only in the MultiQC report (make_multiqc_fixture.py).
writeLines(c(
  "sample_id\tinput_column\tcondition\tbatch\tlabel\tinputs_available", "s1\tlib_A1\ttreated\tb1\tA1 treated\t",
  "s2\tlib_A2\ttreated\tb2\tA2 treated\t", "s3\tlib_B1\tcontrol\tb1\tB1 control\t",
  "s4\tlib_B2\tcontrol\tb2\tB2 control\t", "s5\tlib_C1\t\t\tC1 QC only\tqc"
), file.path(a, "samples.tsv"))
files_a <- setdiff(sort(list.files(a, recursive = TRUE)), c("md5sum.txt", "samples.tsv"))
writeLines(sprintf("%s  ./%s", vapply(file.path(a, files_a), md5_file, ""), files_a), file.path(a, "md5sum.txt"))

## format variants for adapter tests -----------------------------------------------------------------------
fmt <- file.path(root, "formats")
dir.create(fmt, showWarnings = FALSE)
writeLines(sprintf("%s\t%d", genes$ensembl[1:20], cnt_a[1:20, 1]), file.path(fmt, "star_2col.txt"))
writeLines(sprintf("%s\t%d\t%d\t%d", c(genes$ensembl[1:5], genes$ensembl[2]), 1:6, 0L, 1:6), file.path(fmt, "star_duplicated.txt"))

## scenario b: raw and normalised counts in Excel ------------------------------------------------------------
b <- file.path(root, "scenario_b")
dir.create(b, showWarnings = FALSE)
cols_b <- c("D1_control", "D2_control", "D3_control", "D1_treated", "D2_treated", "D3_treated")
cnt_b <- simulate(6, c(2e5, 2.5e5, 1.8e5, 2.2e5, 2.6e5, 2e5), c(up, down), 4:6, rep(c(1.5, -1.5), each = 15))
sym_b <- ifelse(is.na(genes$symbol), genes$ensembl, genes$symbol)
raw <- data.frame(Geneid = sym_b, cnt_b, check.names = FALSE)
names(raw)[-1] <- cols_b
keep <- rowSums(cnt_b >= 10) >= 3
sf <- DESeq2::estimateSizeFactorsForMatrix(cnt_b[keep, ])
norm <- data.frame(Gene_ID = sym_b[keep], sweep(cnt_b[keep, ], 2, sf, "/"), check.names = FALSE)
names(norm)[-1] <- cols_b
write_xlsx_fixture <- function(sheets, path) {
  wb <- openxlsx::createWorkbook(creator = "fixture")
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]])
  }
  tmp <- tempfile(fileext = ".xlsx")
  openxlsx::saveWorkbook(wb, tmp, overwrite = TRUE)
  ex <- tempfile("fx_")
  utils::unzip(tmp, exdir = ex)
  # Providers' workbooks often declare a wrong used range; readers that trust it see one row only.
  for (f in list.files(file.path(ex, "xl", "worksheets"), pattern = "\\.xml$", full.names = TRUE)) {
    x <- readLines(f, warn = FALSE)
    writeLines(gsub('<dimension ref="[^"]*"/>', '<dimension ref="A1"/>', x), f)
  }
  rezip <- tempfile(fileext = ".xlsx")
  withr::with_dir(ex, zip::zip(rezip, files = sort(list.files(ex, recursive = TRUE, all.files = TRUE)), mode = "mirror"))
  normalise_zip(rezip, path)
}
write_xlsx_fixture(list(Sheet1 = raw), file.path(b, "raw_counts.xlsx"))
donor <- sub("^D([0-9])_.*$", "\\1", cols_b)
writeLines(c(
  "sample_id\tinput_column\tcondition\tdonor\tdonor_source",
  sprintf(
    "d%s_%s\t%s\t%s\t%s\tinferred from sample name [I]", donor, ifelse(grepl("control", cols_b), "ctrl", "trt"),
    cols_b, sub("^D[0-9]_", "", cols_b), donor
  )
), file.path(b, "samples.tsv"))
write_xlsx_fixture(list(NormalizedCount = norm), file.path(b, "normalized_counts.xlsx"))

## scenario c: DE tables only -----------------------------------------------------------------------------------
cc <- file.path(root, "scenario_c")
dir.create(cc, showWarnings = FALSE)
cd <- data.frame(condition = factor(rep(c("control", "treated"), each = 3)), row.names = cols_b)
dds <- DESeq2::DESeqDataSetFromMatrix(cnt_b[keep, ], cd, ~condition)
dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
r <- as.data.frame(DESeq2::results(dds, contrast = c("condition", "treated", "control")))
fc <- 2^r$log2FoldChange
de <- data.frame(
  Gene_ID = sym_b[keep], baseMean = r$baseMean, log2FoldChange = r$log2FoldChange, FC = fc,
  FoldChange = ifelse(r$log2FoldChange < 0, -1 / fc, fc), lfcSE = r$lfcSE, stat = r$stat,
  pvalue = r$pvalue, padj = r$padj
)
de <- de[!is.na(de$padj), ]
sig <- de[de$padj <= 0.05, ]
write_xlsx_fixture(
  list(
    DESeq2_results = de, "padj<=0.05" = sig,
    "padj<=0.05+FoldChange>=1.5" = sig[sig$FoldChange >= 1.5, ],
    "padj<=0.05+FoldChange<=-1.5" = sig[sig$FoldChange <= -1.5, ],
    Statistic = data.frame(padj0.05 = nrow(sig), UP = sum(sig$FoldChange >= 1.5), Down = sum(sig$FoldChange <= -1.5))
  ),
  file.path(cc, "comparison.xlsx")
)
m_lfc <- r$log2FoldChange[!is.na(r$padj)]
noiseq <- data.frame(
  GeneID = de$Gene_ID, control_mean = round(de$baseMean * 0.9, 4), treated_mean = round(de$baseMean * 1.1, 4),
  M = round(m_lfc, 6), D = round(abs(m_lfc) * 10, 4), prob = round(stats::plogis(abs(m_lfc) * 3 - 2), 6),
  ranking = round(m_lfc * 5, 4)
)
data.table::fwrite(noiseq, file.path(cc, "noiseq_results.tsv"), sep = "\t")

cat("fixtures written:\n")
for (f in sort(list.files(root, recursive = TRUE))) cat(sprintf("%8d  %s\n", file.size(file.path(root, f)), f))
