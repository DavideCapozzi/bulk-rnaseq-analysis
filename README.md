# bulkrnaseq

Reproducible analysis of the files an RNA-seq provider delivers: gene counts, abundances, normalised counts,
differential-expression tables, coverage tracks and QC reports. An R package plus a Snakemake workflow turn them
into a fixed, documented set of tables and figures.

One YAML file and one sample sheet describe a project. Each analysis module declares the data it needs, so the same
workflow serves a project with raw STAR counts and a project with only a DESeq2 spreadsheet. Modules whose inputs
are missing are skipped with a stated reason, never silently.

Statistical guardrails live in the code, not in documentation:

- no sample-level p-value without residual degrees of freedom (`assert_df()`: error at 0, exploratory flag at 1);
- every p-value carries an `inference_scope`: `sample` (samples are the replicates), `conditional_on_contrast`
  (genes are the replicates, e.g. ORA or NOISeq without replicates) or `ranking_only` (gene permutations, fgsea);
- every random step is seeded, and two runs give byte-identical outputs;
- sample metadata come only from the sample sheet, never from file or column names.

Status: early development (0.x). [`docs/MODULES.md`](docs/MODULES.md), generated from the module registry
([`inst/modules.yaml`](inst/modules.yaml)) and the config schema, lists every module with its status, the data it
needs, its outputs, caveats and parameters; [`docs/CONFIG.md`](docs/CONFIG.md) documents every other section of
the config. Modules not implemented yet are skipped with the reason "not
implemented"; inputs not supported yet (coverage tracks, md5sum files) are refused rather than ignored. MultiQC HTML
reports (1.x with embedded plot data) are read as QC input; decoding them needs `python3`. Design notes:
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md). Changes: [`CHANGELOG.md`](CHANGELOG.md).

## Requirements

Linux, and a conda/micromamba environment with R 4.4, Bioconductor (DESeq2, edgeR, SummarizedExperiment) and
Snakemake. The exact environment used for development and CI is pinned in
[`envs/bulk-rnaseq-analysis-dev.explicit.txt`](envs/bulk-rnaseq-analysis-dev.explicit.txt):

```bash
micromamba create -n bulkrnaseq -f envs/bulk-rnaseq-analysis-dev.explicit.txt
export BULKRNASEQ_ENV=$(micromamba env list --json | python3 -c 'import json,sys; print([p for p in json.load(sys.stdin)["envs"] if p.endswith("/bulkrnaseq")][0])')
```

The launcher takes Snakemake and Rscript from `BULKRNASEQ_ENV`. The workflow loads the package from the checkout,
so installing it is needed only for interactive use (`R CMD INSTALL .`).

## Quickstart on the synthetic fixtures

Three example configs in [`config/examples/`](config/examples/) run on small synthetic data shipped in
`tests/fixtures/` (a test checks that the last column matches what `plan` resolves):

| Config | Inputs | Modules that run |
|---|---|---|
| `scenario_a_star.yaml` | STAR counts, StringTie abundances, a MultiQC report | `qc`, `pca`, `markers`, `contamination`, `reference`, `contrasts`, `de_test`, `concordance`, `geneset_test`, `volcano`, `signatures`, `pairing`, `power`, `reference_axis`, `de_heatmap` |
| `scenario_b_excel.yaml` | raw and normalised counts in Excel | `qc`, `pca`, `markers`, `de_test`, `geneset_test`, `volcano`, `pairing`, `power`, `de_heatmap` |
| `scenario_c_de_only.yaml` | a DESeq2 and a NOISeq table only | `geneset_test`, `ora`, `volcano` |

```bash
bin/bulkrnaseq plan  config/examples/scenario_b_excel.yaml   # what runs, what is skipped and why
bin/bulkrnaseq fetch config/examples/scenario_b_excel.yaml   # copy resources into the cache (explicit step)
bin/bulkrnaseq run   config/examples/scenario_b_excel.yaml
```

Outputs go to `results/<project name>/`: one folder per step (`01_ingest`, `02_qc`, …) with TSV tables, a
`.meta.json` sidecar per table (description, units, filters, seed), PNG figures, one Excel workbook per module in
`report/`, `report/summary.md` with the status of every module and the reason for each skip, and
`report/methods.md`, a methods text of the modules that ran. With `export.researcher_workbook` the main results
also go to one workbook for the researcher, with a README sheet that explains every column.
`00_provenance/` records the resolved config, the toolkit commit, the md5 of every input and, in `benchmarks/`,
the wall time and peak memory of every step. `--cores N` (before the config) runs independent modules in
parallel, by default up to 4; the outputs do not depend on it.

### Scenario a: provider counts from STAR

Four libraries (two conditions in two batches) with STAR `ReadsPerGene` counts and StringTie abundances, plus a
fifth library that has only MultiQC metrics (`inputs_available: qc`): it appears in `02_qc/tables/qc_summary.tsv`
and nowhere else. The config shows most of the toolkit on one project:

- `contrasts` estimates one-sample-per-group contrasts within each batch, without p-values; `concordance`
  compares them (`11_concordance/tables/concordance_summary.tsv`, with the samples the two contrasts share);
- `de_test` fits `~ batch + condition`: 1 residual degree of freedom, so its results are flagged `exploratory`
  in `contrast_stats` and in `volcano_points`;
- `contamination` builds a gene mask from synthetic gene sets and profiles, and `reference` a gene universe from
  synthetic HPA exports; `pca`, `contrasts`, `concordance` and `power` then repeat their analysis on gene subsets
  such as `no_masked` and `in_reference_universe+no_masked` (each table has a column naming the subset);
- `geneset_test` runs cameraPR and fgsea on synthetic gene sets for each contrast and lists the sets significant
  in both contrasts of a pair; fgsea p-values carry `inference_scope = ranking_only` (genes are permuted, not
  samples);
- `reference_axis` correlates the fold changes with an axis built from two synthetic HPA cell types, averages them
  over two gene sets, tests the sign agreement of the contrast pair and scores each sample on the axis;
- `markers` reports two panels per sample (one symbol is missing on purpose and listed in the metadata) and a
  ratio of one gene between two samples; `signatures` scores synthetic sets and a custom set with three methods
  in two universes, with their difference per contrast and their sign agreement in the pair;
- `de_heatmap` draws the most significant genes of `de_test`;
- the researcher workbook (`report/results_workbook.xlsx`) collects the results with a README sheet;
- `ora` is skipped: the config has no `kegg_pathways` resource.

### Scenario b: count matrices in Excel

Six samples (three donors, two conditions) as raw and provider-normalised count matrices with gene symbols as
IDs. At ingest, `01_ingest/normalization_check.tsv` checks whether DESeq2 median-of-ratios on a filtered gene set
reproduces the provider's normalised counts. `de_test` fits an unpaired (`~ condition`) and a paired
(`~ donor + condition`) design; `10_de_test/tables/design_summary.tsv` gives the residual degrees of freedom of
each. `geneset_test` ranks the genes by the DESeq2 Wald statistic (`ranking_statistic: stat`) of `de_test`
(cameraPR, the default). `markers` reports one panel with the provider's normalised counts next to the CPM;
`de_heatmap` shows the provider's normalised counts of the top genes, half up and half down; the researcher
workbook has the DE, gene set, marker and QC sheets. Without `contrasts:`, `contrast_pairs:` or contamination and reference profiles in the config, those modules
are skipped and `report/summary.md` says why.

### Scenario c: differential-expression tables only

No counts: a DESeq2 workbook (`df_residual: 4`) and a NOISeq table without replicates (`df_residual: 0`). Every
module that needs counts is skipped with the reason. `volcano` and `ora` (GO, KEGG, Reactome) run on both tables,
`geneset_test` on the DESeq2 one (the first method; `modules.geneset_test.estimator` chooses another);
the NOISeq probabilities carry `inference_scope = conditional_on_contrast`, because with 0 residual degrees of
freedom they say nothing about new samples. The researcher workbook has one set of DE sheets per table.

## A new project

```bash
bin/bulkrnaseq init ~/projects/my_study --scenario excel --git
```

This creates `input/` (provider files, kept read-only), `analysis/pipeline_v1/config/` (config and sample sheet)
and, with `--git`, a local repository that tracks the config only. Then:

1. copy the provider files into `input/` (every run records their md5 in `00_provenance/run_manifest.json`);
2. fill in the sample sheet: `sample_id`, `input_column` (the file name part or matrix column of each sample) and
   the experimental factors; nothing is parsed from sample names;
3. point `inputs` at the files ([`docs/CONFIG.md`](docs/CONFIG.md)); a provider DE table needs `df_residual`, the
   residual degrees of freedom of the provider's model, which decides whether its p-values are sample-level;
4. `plan` shows which modules run and why the others are skipped; add `designs`, `contrasts` or module parameters
   ([`docs/MODULES.md`](docs/MODULES.md)) until the plan matches the questions of the study;
5. `fetch` retrieves the resources (GTF, gene sets, HPA exports) into the cache, then `run`.

A typical case: a provider delivers raw counts, normalised counts and a DESeq2 workbook in Excel, for two sorted
fractions of three donors. The `matrix` and `deseq2` adapters read them; `01_ingest/normalization_check.tsv`
reconstructs the provider's gene filter and size factors; `de_test` refits the provider's design
(`compare_with` reports gene by gene whether it is reproduced) and, once the pairing of the samples is confirmed,
a paired `~ donor + fraction` design; `volcano`, `ora` and `geneset_test` run on the provider table and on the
refitted ones.

Each config pins `toolkit.version`. The run stops if it differs from this checkout's `DESCRIPTION` or if the
checkout has uncommitted changes, so every results folder maps to one commit. `run --dev` turns both checks into
warnings and records that in the run manifest.

## Development

```bash
Rscript -e 'testthat::test_local()'      # unit tests (synthetic data only)
python -m pytest -q                      # launcher and workflow helpers
pre-commit install                       # data guard, ruff, styler, lintr on every commit
```

The data guard (`tools/check_no_data.py`) rejects data files, absolute paths, files over 500 KB and
project-specific patterns kept outside the repository. Real projects and their equivalence checks against earlier
analyses live outside this repository.

## License

MIT, see [`LICENSE`](LICENSE).
