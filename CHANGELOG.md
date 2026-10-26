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
- `bulkrnaseq plan` and `run` stop before any analysis when an input file is missing, listing the paths.
- When a step fails, `bulkrnaseq run` prints the last lines of its log (the R error) after Snakemake's report.
- README, this changelog and a CI workflow (data guard, ruff, pytest, styler, lintr, testthat, R CMD check,
  synthetic scenarios run twice with an md5 comparison).

### Changed

- Config schema: `modules.pca.exclude_symbol_regex` is replaced by `symbol_exclusion_variants` (a list: each
  regex adds a variant to the grid); new `modules.pca.orient_pc1`, `abundance_sex_genes`,
  `modules.qc.figure_metrics`, `modules.de_test.results_alpha`; `modules.de_test.compare_with` items are now
  `{design, contrast_table}` objects.
- `de_test` writes its results as `contrast_stats` (registry outputs updated), like the other producers of
  contrasts.
- Config schema is strict where the code is implemented: unknown parameters of the implemented modules
  (`qc`, `pca`, `contrasts`, `de_test`, `concordance`, `ora`, `volcano`) and of `designs` / `contrast_pairs`
  items are errors instead of being ignored. Parameters that no code read are removed until their analyses are
  ported (`modules.contrasts.norm_variants`, `prior_variants`, `thinning`, `key_genes`;
  `modules.concordance.strata`, `deciles`, `concordant_gene_lists`, `power_simulation`); `designs[].engine`
  accepts only `deseq2`; `modules.contrasts.spread_subsets` is declared.
- `inputs.qc`, `inputs.coverage` and `inputs.checksums` are refused with "not implemented" instead of being
  ignored by the ingest (scenario a no longer declares them).
- A module that is not implemented is reported as such first, before any missing resource or condition.
- Module defaults live only in the schema (`with_defaults()`), also when module functions are called directly.
- The module registry lists only the tables and figures the modules write; a test checks it.
- `write_module_outputs()` refuses a table with p-value columns and no `inference_scope` column.
- Methodological: imported DE tables with a significance metric need `df_residual` (residual degrees of freedom
  of the provider's model). Their p-values or posterior probabilities are `inference_scope = sample` only with
  1 or more, and `conditional_on_contrast` with 0 (e.g. NOISeq without replicates); before, every imported
  table was labelled `sample`. Configs must add `df_residual` to `inputs.contrast_tables`.
- Methodological: contrast tables carry an `exploratory` column (`TRUE` when the model had 1 residual degree of
  freedom: `de_test`, the `deseq2_apeglm` estimator, imported tables), copied into `volcano_points` and the ORA
  tables, so results built on a 1-df fit say so where they are read; before, the flag was only in
  `design_summary`.
- Methodological: `concordance_summary` has a `shared_samples` column (samples behind both contrasts of a
  pair; `NA` when one is imported) and the module a caveat: contrasts that share samples share their noise,
  which inflates concordance.

### Removed

- `read_gtf_exons()` (unused until the coverage and covariates modules are ported).
- The `repo` argument of `write_run_manifest()` and `run_step_cli()`: the toolkit identity comes from
  `00_provenance/toolkit.json`.

### Fixed

- `pca`: without `modules.pca.separation_factors` the `pc*_separates` labels and the colours of the PCA figure
  were empty, because the schema default (an empty list) was taken as "no factor"; it now means every factor of
  `samples.factors`, as documented. Configs that set `separation_factors` are unaffected.
- Workflow: the `ingest` rule did not depend on `00_provenance/toolkit.json`, so a change to the ingest or
  adapter code alone did not re-run it and its outputs could come from code other than the recorded commit.
- `run_manifest.json`: the adapter parameters of each input were written with doubled quotes.
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
