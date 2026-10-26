# bulkrnaseq

Reproducible analysis of the files an RNA-seq provider delivers: gene counts, abundances, normalised counts,
differential-expression tables, coverage tracks and QC reports. An R package plus a Snakemake workflow turn them
into a fixed, documented set of tables and figures.

One YAML file and one sample sheet describe a project. Each analysis module declares the data it needs, so the same
workflow serves a project with raw STAR counts and a project with only a DESeq2 spreadsheet. Modules whose inputs
are missing are skipped with a stated reason, never silently.

Statistical guardrails live in the code, not in documentation:

- no sample-level p-value without residual degrees of freedom (`assert_df()`: error at 0, exploratory flag at 1);
- every p-value carries an `inference_scope` (`sample` or `conditional_on_contrast`);
- every random step is seeded, and two runs give byte-identical outputs;
- sample metadata come only from the sample sheet, never from file or column names.

Status: early development (0.x). Implemented modules: `qc`, `pca`, `contrasts`, `de_test`, `concordance`, `ora`,
`volcano`. The others listed in [`inst/modules.yaml`](inst/modules.yaml) are skipped with the reason "not
implemented" until they are ported; inputs not supported yet (MultiQC reports, coverage tracks, md5sum files) are
refused rather than ignored. Design notes: [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md). Changes: [`CHANGELOG.md`](CHANGELOG.md).

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
`tests/fixtures/`:

| Scenario | Inputs | Modules that run |
|---|---|---|
| `scenario_a_star.yaml` | STAR counts and StringTie abundances | qc, pca, contrasts, de_test, concordance, volcano |
| `scenario_b_excel.yaml` | raw and normalised counts in Excel | qc, pca, de_test (unpaired and paired), volcano |
| `scenario_c_de_only.yaml` | a DESeq2 and a NOISeq table only | ora, volcano |

```bash
bin/bulkrnaseq plan  config/examples/scenario_b_excel.yaml   # what runs, what is skipped and why
bin/bulkrnaseq fetch config/examples/scenario_b_excel.yaml   # copy resources into the cache (explicit step)
bin/bulkrnaseq run   config/examples/scenario_b_excel.yaml
```

Outputs go to `results/example_excel/`: one folder per step (`01_ingest`, `02_qc`, …) with TSV tables, a
`.meta.json` sidecar per table (description, units, filters, seed), PNG figures, and `report/summary.md` with the
status of every module. `00_provenance/` records the resolved config, the toolkit commit and the md5 of every
input.

## A new project

```bash
bin/bulkrnaseq init ~/projects/my_study --scenario excel --git
```

This creates `input/` (provider files, kept read-only), `analysis/pipeline_v1/config/` (config and sample sheet)
and, with `--git`, a local repository that tracks the config only. Fill in the sample sheet (`sample_id`,
`input_column` and the experimental factors), point `inputs` at the provider files, then `plan`, `fetch`, `run`.

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
