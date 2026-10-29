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
- MultiQC reports as QC input (`inputs.qc.multiqc`): `read_multiqc()` decodes the plot data embedded in MultiQC
  1.x HTML reports with `inst/python/multiqc_decode.py` (Python standard library; no R package decodes LZString),
  writes it to `01_ingest/multiqc/` and derives per-sample metrics from the STAR, RSeQC and raw-read FastQC
  sections (`raw_read_pairs`, `pct_unique`, `assigned_reverse`, `genebody_3to5`, `mean_q_r1`,
  `adapter_pct_last_pos_r1`, ...) into the QC table, hence `qc_summary` and `pc1_vs_qc`. Sample names are matched
  through `input_column` or `inputs.qc.sample_aliases`; an ambiguous or unmatched name is an error.
- Samples with `inputs_available: qc` (no expression data, e.g. libraries that never reached the counts) may leave
  factor columns empty and appear only in the QC tables.
- `pc1_vs_qc.meta.json` records `n_samples` and the smallest two-sided Spearman permutation p attainable
  (`2 / n!`); the registry caveats say that with 4 samples the table lists candidates, not drivers.
- `pairing` module: sex genes (XIST, chrY genes) and expression fingerprints of genotype-sensitive genes per
  pair of samples (`fingerprint_variants`: gene lists, CPM or TMM-CPM, minimum reads), male RNA fraction relative
  to configured all-male reference samples with optional denominators (`sex_denominators`), and
  `pairing_summary`: whether the pairs of a hypothesised donor or pool (`group_column`) are the nearest across
  `across_factor`, with the chance of that outcome for an unrelated grouping (at most 1/k! for k groups). The
  outputs are expression clues and say so; variants without genes are skipped and noted. New gene lists
  `fingerprint_genes`, `fingerprint_genes_extended`, `pan_t`, `xist` (sources: earlier analysis code, no citation).
- `power` module: BCV of the data with edgeR on a design of the config (common, tagwise, per-gene without
  shrinkage, bootstrap over genes with the project seed), the 95% interval factor of a single gene's BCV for the
  residual df, and sample sizes per group from Hart et al. 2013 over a grid of CV, mean count, fold change and
  alpha (`n_hart()`); `grid.cv_labels` names the assumption behind each CV.
- `contamination` module: profile scores per sample (mean log2(CPM + 1) and CPM sum over a gene list, CPM of all
  genes; `rel_to_max` for `reference_profile`), per-gene CPM of the profiles, `ratio_checks` (log2 CPM ratios between
  two samples over a profile's genes), `excess_checks` (CPM of genes above baseline samples per unit of a profile's
  CPM, and the ratio between two samples) and `contaminant_mask`: genes of `mask.include` gene-set groups of an
  MSigDB export and in no `mask.exclude` group, plus the genes of `mask.add_profiles`. The registry caveats say that
  the mask removes known markers only and that profiles chosen from the same samples make it partly data-derived.
- `reference` module: reference expression of configured cell types from Human Protein Atlas exports
  (`hpa_immune_cell`, nTPM; `hpa_single_cell`, nCPM) for every gene, `target_expressed`, `contam_sensitive` and
  `reference_universe` (target expressed and not contamination-sensitive), gene counts per criterion
  (`reference_universe`) and, for `exclusion_candidates`, the first criterion that excludes each gene
  (`universe_exclusions`, which uses the contamination mask when available).
- Data levels `gene_mask` (from `contamination`) and `gene_universe` (from `reference`): per-gene flags
  (`objects/gene_flags.rds`) that `pca`, `contrasts`, `concordance` and `power` use in gene subset expressions, e.g.
  `no_masked` or `in_reference_universe+no_masked`. The universe does not include the mask: subsets combine the two
  explicitly, so every table states which genes it used.
- `modules.power.universes`: the BCV is estimated on each configured gene subset, in order, on one random stream.
- `read_gene_sets()` (MSigDB export) and `read_reference_profiles()` (HPA downloads, zipped or not; a missing cell
  type is an error).
- Registry: a module parameter named `resource` names a resource the module needs (checked by the capability
  resolution and verified before the step); the condition `configured` skips a module whose listed parameters are
  all empty (`contamination`, `reference`), with the reason.
- Generic gene list `proliferation` (`inst/extdata/gene_lists/contamination.tsv`).
- Scenario a has a synthetic MultiQC report (`tests/fixtures/make_multiqc_fixture.py`), a QC-only library, and
  synthetic gene sets and HPA exports (`tests/fixtures/make_reference_fixture.py`) for `contamination`, `reference`
  and the gene subsets of `pca`, `contrasts`, `concordance` and `power`.
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
- Config schema is strict where the code is implemented: unknown parameters of the implemented modules (`qc`, `pca`,
  `contrasts`, `de_test`, `concordance`, `ora`, `volcano`, `pairing`, `power`) and of `designs` / `contrast_pairs`
  items are errors instead of being ignored. Parameters that no code read are removed until their analyses are
  ported (`modules.contrasts.norm_variants`, `prior_variants`, `thinning`, `key_genes`;
  `modules.concordance.strata`, `deciles`, `concordant_gene_lists`, `power_simulation`); `designs[].engine` accepts
  only `deseq2`; `modules.contrasts.spread_subsets` is declared.
- Config schema is strict for `modules.contamination` and `modules.reference`; their draft parameters are replaced
  (`contamination.profiles` is a list of `{id, genes}`, `ratio_checks` / `excess_checks` / `mask` are structured,
  `mask.resource` names the gene-set export; `reference.target_cell`, `contaminant_cells`, `axes` and
  `universe_grid` give way to `profiles`, `target`, `contaminants` and `exclusion_candidates`).
- `reference` no longer lists `reference_axis` among its outputs nor `contrast_stats` among its inputs: the axis
  analysis reads fold changes and the universe, and will run downstream of `contrasts` (which now reads the universe).
- `inputs.coverage` and `inputs.checksums` are refused with "not implemented" instead of being ignored by the
  ingest. `inputs.qc` loses `sample_name_pattern` and `tables` (never read) and gains `multiqc[].read` and
  `sample_aliases`.
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

- `gene_subset()`: with a gene without expression value (e.g. a gene tested only by `de_test`) the `ge_median`,
  `ge_q25` and `ge_q75` strata were empty, because the median was `NA`; such genes are now left out of the strata.
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
