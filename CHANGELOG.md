# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/) (0.x: the config schema may still change between minor versions).

## [Unreleased]

### Added

- `qc` module: expression QC per sample on all genes (library size, detected genes, shares of mitochondrial,
  rRNA, ribosomal-protein, protein-coding, lncRNA and hemoglobin reads and of the most expressed genes; TPM
  metrics when abundances are present), most expressed genes per sample, a per-sample summary that joins
  imported QC metrics, and the `qc_metrics` figure.
- `pca` module: a grid of PCAs over sex-gene variants, symbol exclusions, original or thinned counts, gene
  subsets, transforms (VST, log CPM TMM) and numbers of most variable genes, with percent variance and a
  descriptive "separates" label per PC; PC scores, PCA on log2(TPM + 1), sample correlations and distances, PC1
  drivers (correlation with QC metrics, loadings vs gene span, lncRNA status and expression) and figures.
- `de_test` module: DESeq2 Wald tests for every design with a `contrast`, guarded by `assert_df()`, with
  `inference_scope = sample`, a design summary, and `compare_with` to check a fit against an imported provider
  table (`provider_comparison`).
- `check_normalization()`: at ingest, when raw and provider-normalised counts are both configured, reports the
  implied size factors, whether DESeq2 median-of-ratios reproduces them, and which "at least k reads in at least
  n samples" rule yields the provider's gene set (`01_ingest/normalization_check.tsv`, ingest report).
- Generic gene list `hemoglobin` (`inst/extdata/gene_lists/expression_qc.tsv`).
- Modules may return `figure_sizes`; the ingest QC table is passed to modules as `inputs$qc`.
- README, this changelog and a CI workflow (data guard, ruff, pytest, styler, lintr, testthat, R CMD check,
  synthetic scenarios run twice with an md5 comparison).

### Changed

- Config schema: `modules.pca.exclude_symbol_regex` is replaced by `symbol_exclusion_variants` (a list: each
  regex adds a variant to the grid); new `modules.pca.orient_pc1`, `abundance_sex_genes`,
  `modules.qc.figure_metrics`, `modules.de_test.results_alpha`; `modules.de_test.compare_with` items are now
  `{design, contrast_table}` objects.
- `de_test` writes its results as `contrast_stats` (registry outputs updated), like the other producers of
  contrasts.

### Fixed

- `concordance` failed when contrasts from several producers (e.g. `contrasts` and `de_test`) were present: a
  pair is now evaluated for every method that estimates both of its contrasts, and fails only when no method
  covers both.

## [0.1.0] - 2026-10-21

### Added

- First release: `bulkrnaseq` R package with adapters (STAR ReadsPerGene, StringTie, matrix tables in TSV/CSV/
  Excel, DESeq2/NOISeq/generic DE tables, GTF), canonical objects and validators, capability resolution from
  `inst/modules.yaml`, the `contrasts`, `concordance`, `ora` and `volcano` modules, deterministic writers
  (TSV, gzip, Excel, PNG), provenance, a Snakemake workflow with explicit resource fetching, the `bin/bulkrnaseq`
  launcher (`init`, `plan`, `fetch`, `run`, `version`) and toolkit version pinning.
