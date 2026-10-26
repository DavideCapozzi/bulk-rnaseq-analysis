# Architecture

Status: 0.x, partly implemented. Modules implemented: `qc`, `pca`, `contrasts`, `de_test`, `concordance`, `ora`,
`volcano`; the others are registered and skipped as "not implemented". Sections marked **[I]** are design intentions
not yet implemented or verified; **[V]** statements cite a file in this repository.

`bulkrnaseq` is an R package plus a Snakemake workflow that turns the outputs of an RNA-seq provider
(gene counts, abundances, normalised counts, differential-expression tables, coverage tracks, QC reports)
into a fixed, documented set of tables and figures per project. One YAML file and one sample sheet describe a
project; the workflow activates only the analysis modules whose input requirements are met.

The design goals, in order:

1. **Reproduce first, change later.** A refactored analysis must reproduce the reference tables of the analysis
   it replaces before any methodological change. Changes are separate commits with a CHANGELOG entry.
2. **Statistical guardrails in code.** Designs without residual degrees of freedom cannot produce sample-level
   p-values; gene-permutation tests are labelled as ranking tools; every random step is seeded and recorded.
3. **Input-level agnostic.** A project may provide anything from raw STAR counts to a single DE spreadsheet.
4. **No hidden state.** No downloads inside analysis rules, no absolute paths, no sample metadata parsed from
   column names, no constants in code: thresholds, gene lists and seeds live in config or resource files.

---

## 1. Layers

```
                 ┌────────────────────────── Snakemake (workflow/) ──────────────────────────┐
config.yaml ───► │ validate config ─► resolve capabilities ─► rules per module ─► report/export │
samples.tsv      └──────┬───────────────────────────┬─────────────────────────┬──────────────┘
                        │ thin wrappers (workflow/scripts/*.R, *.py) call package functions
                        ▼                           ▼                         ▼
   ┌──────────────┐   ┌──────────────────┐   ┌────────────────────┐   ┌──────────────────────┐
   │  Adapters    │──►│ Canonical objects│──►│  Analysis modules  │──►│ Writers and reports  │
   │ (I/O only)   │   │ + validate_*()   │   │ (pure functions)   │   │ (I/O only)           │
   └──────────────┘   └──────────────────┘   └────────────────────┘   └──────────────────────┘
    read_star_counts    SummarizedExperiment   run_<module>(objects,     write_table() + .meta.json
    read_stringtie      contrast table         params) → list(tables,    write_figure()
    read_matrix_table   QC long table          figures, meta)            write_module_xlsx()
    read_de_table       coverage query                                   write_run_manifest()
    read_multiqc        resource objects                                 write_ingest_report()
    read_bedgraph_regions
    check_md5sum, read_gtf, read_gene_sets, read_reference_profiles
```

Rules:

- **Adapters** turn one file format into one canonical object. They are the only package functions besides
  writers that touch the file system. They never infer sample metadata from file or column names: the sample
  sheet maps every input file or column to a `sample_id` through `input_column`.
- **Canonical objects** have a validator (`validate_se()`, `validate_contrast_table()`, `validate_qc_table()`,
  `validate_coverage_query()`) that every adapter calls before returning and every module calls on entry.
- **Modules** are pure: inputs are canonical objects, a parameter list and resource objects; output is
  `list(tables = <named list of data.frames>, figures = <named list of ggplot/grob>, meta = <list>)`.
  No `setwd`, no global options, no file access, no `library()` calls, no column access by position.
- **Writers** serialise tables (TSV, gzip for matrices) with a sidecar `<name>.meta.json`, figures (PNG with
  fixed device and size), per-module Excel workbooks with an `Info` sheet, and the run manifest.
- **Wrappers** in `workflow/scripts/` contain only argument handling and one call into the package
  (`bulkrnaseq::run_module_cli()`, `bulkrnaseq::run_ingest_cli()`), so all logic is testable with `testthat`.

## 2. Canonical objects

### 2.1 Sample level: `SummarizedExperiment`

| Slot | Content |
|---|---|
| assays | Named subset of `counts` (integer, never derived from normalised values), `tpm`, `fpkm`, `normalized` (provider-normalised values, kept as delivered). Missing assays are absent, not filled. |
| rowData | `feature_id` (input ID, unique), `gene_symbol`, `ensembl_gene_id`, `entrez_id`, `chromosome`, `biotype`, `start`, `end`, `strand`, `length` when available; mapping columns come from the local annotation. |
| colData | Sample sheet columns, `sample_id` as row names; factor levels as declared in the config. |
| metadata | `provenance`: list of `{file, md5, adapter, params}`; `organism`; `id_type`; `ingest_report` (mapping losses, duplicates, filters). |

Rule: normalised values are **never** converted back to counts. A module that needs counts declares the
`counts` level and is skipped when only `normalized` is available.

### 2.2 Contrast level: long contrast table (`data.frame`)

Fixed columns, in this order:

| Column | Type | Meaning |
|---|---|---|
| `contrast_id` | chr | ID from the config |
| `feature_id` | chr | input gene ID |
| `gene_symbol` | chr | mapped symbol (NA if unmapped) |
| `lfc` | dbl | log2 fold change, numerator over denominator |
| `lfc_se` | dbl | standard error or NA |
| `stat` | dbl | test statistic or NA |
| `pvalue` | dbl | NA unless a valid test produced it |
| `padj` | dbl | adjusted p-value; never a rescaled posterior probability |
| `base_mean` | dbl | mean of normalised counts or NA |
| `significance_metric` | chr | `padj`, `posterior_prob` or `none` |
| `method` | chr | e.g. `logcpm_difference`, `deseq2_wald`, `deseq2_apeglm`, `imported_deseq2`, `imported_noiseq` |
| `df_residual` | int | residual degrees of freedom of the producing model; NA if unknown (imported) |
| `source` | chr | `computed` or `imported:<file>` |

Extension columns are allowed after the fixed ones when the producer documents them (e.g. `posterior_prob`
for NOISeq, `ave_log2cpm` for the log-CPM estimator, `fold_change_signed` kept from provider tables).
`significance_values(ct)` returns `padj` or `posterior_prob` according to `significance_metric`, so ORA and
volcano plots never need to relabel one as the other.

### 2.3 Coverage

Coverage tracks (STAR bedGraph, one file per sample × multiplicity × strand) are never loaded whole.
`read_bedgraph_regions(path, regions)` streams the file once, keeps only intervals overlapping the requested
regions (half-open, 0-based, as in bedGraph) and returns the total signal of the file in the same pass.
Implementation: `data.table::fread(cmd = <awk filter>)`, so R receives only matching lines; `awk` is a
documented system dependency. Strand assignment (`str1`/`str2` → genomic strand) is either declared in the
config or inferred from reference genes on known strands and written to a table.

### 2.4 QC: long table

`sample_id`, `metric`, `value` (dbl), `source` (`multiqc:<report>:<plot>`, `star_counts:N_*`, `expression`,
`imported:<file>`). Wide views are produced by writers only.

### 2.5 Resources

Annotation (gene and exon tables from a GTF), gene sets (MSigDB export), reference expression profiles
(HPA immune cell / single cell), transcript FASTA, KEGG pathway tables, curated gene lists. All are read from a
per-project cache directory listed in `00_provenance/resources_manifest.tsv` (path, md5, version, retrieval date,
source URL). Only the `fetch_resources` rule writes the cache; analysis rules fail if a resource is missing or its
md5 differs from the manifest.

## 3. Data levels and capabilities

Every module declares in `inst/modules.yaml` the data levels it requires, the resources it needs, extra
conditions (e.g. at least two contrasts) and the levels it provides. Levels:

| Level | Provided by |
|---|---|
| `counts` | STAR ReadsPerGene (2 or 4 columns), generic matrix with `kind: counts` (integers checked) |
| `abundance` | StringTie gene_abundance (TPM, FPKM), generic matrix with `kind: tpm`/`fpkm` |
| `normalized` | generic matrix with `kind: normalized`; never upgraded to `counts` |
| `contrast_stats` | imported DE tables; modules `contrasts` and `de_test` |
| `ranked` | derived from any `contrast_stats` with a per-gene statistic (`lfc` or `stat`) |
| `coverage` | bedGraph |
| `qc` | MultiQC, STAR `N_*` rows, imported QC tables |

The Snakefile resolves capabilities before building the DAG: start from the levels provided by the inputs,
add levels provided by enabled modules whose requirements are met, repeat until stable. Modules whose
requirements are not met are written to the run manifest as
`{"module": "concordance", "status": "skipped", "reason": "needs >= 2 contrasts with per-gene statistics; input provides 1 (<file>)"}`;
they never fail the run silently or abort it.

### 3.1 Module × level matrix

`R` = required, `o` = optional (used if present), `—` = not used. Resources in the last column.

| NN | Module | counts | abundance | normalized | contrast_stats | ranked | coverage | qc | Resources / conditions |
|---|---|---|---|---|---|---|---|---|---|
| 02 | `qc` | R¹ | o | — | — | — | — | R¹ | annotation; ¹ either level suffices |
| 03 | `pca` | R | o | — | — | — | — | o | annotation; ≥ 3 samples |
| 04 | `markers` | R | o | — | — | — | — | — | annotation, marker panels |
| 05 | `contamination` | R | — | — | — | — | — | — | annotation, contaminant profiles, gene sets (mask) |
| 06 | `contamination_sets` | — | — | — | — | R | — | — | contaminant profiles |
| 07 | `reference` | R² | — | — | o | — | — | — | annotation, reference profiles; ² or any gene list from `contrast_stats` |
| 08 | `deconvolution` | o³ | R³ | — | — | — | — | — | reference profiles; ³ TPM default, CPM variant |
| 09 | `contrasts` | R | — | — | — | — | — | — | contrasts in config; provides `contrast_stats`, `ranked` |
| 10 | `de_test` | R | — | — | — | — | — | — | design with df ≥ 1 (`assert_df`); provides `contrast_stats` |
| 11 | `concordance` | o | — | — | R | — | — | — | ≥ 2 contrasts and ≥ 1 contrast pair; power simulation needs `counts` |
| 12 | `geneset_test` | — | — | — | — | R | — | — | gene sets |
| 13 | `ora` | — | — | — | R | — | — | — | OrgDb, GO/Reactome packages, KEGG table in cache |
| 14 | `volcano` | — | — | — | R | — | — | — | — |
| 15 | `signatures` | R | — | — | — | — | — | — | gene sets |
| 16 | `pairing` | R | o | — | — | — | — | — | annotation (chromosomes), sex and fingerprint gene lists |
| 17 | `power` | R | — | — | — | — | — | — | design with df ≥ 1 for BCV; Hart grid always |
| 18 | `coverage` | — | — | — | — | — | R | — | annotation exons |
| 19 | `covariates` | — | — | — | R | — | — | — | annotation exons, transcript FASTA |

Example (scenario c, one DESeq2 table): `contrast_stats` and `ranked` are available; `contamination_sets`,
`geneset_test`, `ora`, `volcano` run; `concordance` runs only if two or more contrast tables are given;
`pca`, `markers`, `contrasts`, `signatures`, `power` (BCV) and per-sample scores are skipped with a reason.

## 4. Statistical guardrails (in code)

- **`assert_df(design, what)`** is called before any sample-level test. The residual degrees of freedom are
  computed from the design matrix (`n samples − rank(X)`), never from the number of groups:
  `df = 0` → error (`"<what>: 0 residual degrees of freedom; only descriptive estimators are allowed"`);
  `df = 1` → warning and `exploratory = TRUE`; `df ≥ 2` → silent.
  Descriptive estimators (log-CPM differences, PCA, scores) do not call it.
- **Inference scope as a column.** Every output table that contains a p-value has a column `inference_scope`
  (also written to the `Info` sheet of the Excel export and to `.meta.json`):
  `sample` — a test whose units are samples (DESeq2/edgeR on a design with df ≥ 1);
  `conditional_on_contrast` — a test whose units are genes for one given contrast (cameraPR, fgsea, ORA, binomial
  sign tests). cameraPR accounts for inter-gene correlation but still tests genes for a given contrast, not
  samples, so on a contrast between single libraries it stays `conditional_on_contrast`.
- **Gene-permutation tests.** Adjusted p-values of pre-ranked fgsea are labelled `ranking_only` (CHANGELOG,
  methodological change): they order sets within a contrast and are not used to select sets downstream.
- **Randomness.** Every function with a random step takes `seed`, runs inside `withr::with_seed()`, and
  writes `seed` and `resampling_unit` (e.g. `genes`, `reads (binomial thinning)`) to the metadata.
- **Imported DE tables are not recomputed.** Thresholds (`padj`, `|fold change|`, posterior probability) come
  from the config; the significance metric keeps its name.
- **Units are declared**: every numeric table column has a unit or scale in its `.meta.json` (`cpm`,
  `log2_cpm_tmm_prior2`, `tpm`, `counts`, `log2_ratio`, `fraction`, `correlation`).

## 5. Gene identifiers and organisms

- The input ID is kept as `feature_id`. Symbols, Ensembl and Entrez IDs are mapped from the cached annotation
  (GTF) and the organism's OrgDb; losses (unmapped, one-to-many, many-to-one) go to `01_ingest/ingest_report.md`.
- Duplicated IDs in count inputs follow `on_duplicate: error | sum | first` (default `error`).
- `inst/organisms.yaml` holds OrgDb, KEGG code, MSigDB species name, Y and mitochondrial chromosome names
  and the mitochondrial symbol prefix. Only `human` is implemented and tested. No other file in the package may
  contain organism-specific strings; a lint test greps for them.

## 6. Output layout

Under `project.output_dir`:

```
00_provenance/
  run_manifest.json        config hash, toolkit identity, package versions, input md5, modules run/skipped (+reason)
  toolkit.json             toolkit version, pinned version, git commit/describe/dirty, code md5, development mode
  resources_manifest.tsv   resource, path, md5, version, retrieved, source
  config_resolved.yaml     config after schema defaults
  run_timestamps.json      wall-clock times (the only file excluded from the md5 identity check)
01_ingest/
  counts.tsv.gz, tpm.tsv.gz, fpkm.tsv.gz, normalized.tsv.gz   wide: feature_id + one column per sample_id
  features.tsv.gz          rowData
  samples.tsv              colData
  contrasts.tsv.gz         imported contrast tables (canonical)
  qc_metrics.tsv           QC long table
  multiqc/                 decoded MultiQC plot data (JSON) and general statistics
  checksums.tsv            md5sum verification (status per file)
  ingest_report.md
NN_<module>/
  tables/<name>.tsv[.gz] + <name>.meta.json
  figures/<name>.png
  logs/<module>.log
report/
  summary.md
  NN_<module>.xlsx         when export.xlsx is true; first sheet "Info" with parameters
```

Module numbers are fixed in `inst/modules.yaml`, so a skipped module leaves a gap rather than renumbering.
Tables use English snake_case column names. Figures are deterministic: fixed sizes and DPI, explicit legend order,
seeded label repulsion. Excel files are written with fixed document properties and zip timestamps so that two
runs produce identical files **[I]**.

## 7. Configuration

- One YAML per project plus a sample sheet TSV, validated with the JSON Schema in `config/schema.yaml` through
  `snakemake.utils.validate` (defaults filled from the schema). The package also validates in R
  (`validate_config()`) so functions can be called without Snakemake.
- Thresholds, priors, cutoffs, seeds, numbers of replicates and gene lists live in the config or in resource files.
- Analysis variants (sensitivity analyses) are parameter lists, not code: e.g.
  `contrasts.norm_methods: [TMM, RLE, upperquartile, none, reference]`, `pca.gene_subsets`, `concordance.strata`.
- `docs/CONFIG.md` is generated from the schema.

### 7.1 Projects and the launcher

The toolkit is one checkout; a project is a folder with data, config and results, and no code:

```
<project>/
  input/                              provider files (read-only)
  analysis/<analysis>/                e.g. pipeline_v1; the Snakemake working folder (.snakemake/)
    config/config.yaml, samples.tsv   versioned by an optional local git repository (config only)
    resources/                        resource cache written by fetch_resources
    results/                          outputs (section 6)
```

`bin/bulkrnaseq` (Python standard library) drives it:

| Command | Effect |
|---|---|
| `init <dir> --scenario star\|excel\|de_only [--git]` | creates the layout from `workflow/templates/`, with `toolkit.version` set to this checkout; `--git` adds a local repository whose `.gitignore` tracks only `analysis/*/config/` |
| `plan <config>` | validates the config and prints, per module, run / skipped / disabled with the reason |
| `fetch <config>` | runs `fetch_resources` |
| `run <config> [-n] [-- snakemake args]` | runs the workflow |
| `version` | toolkit version and git state |

Snakemake is taken from `BULKRNASEQ_ENV` (environment prefix), then from the launcher's or Python's folder, then
from `PATH`; that folder is put first on `PATH`, so `Rscript` comes from the same environment
(`BULKRNASEQ_RSCRIPT` overrides it).

### 7.2 Toolkit version and re-runs

- Every project config pins `toolkit.version`. Before building the DAG the Snakefile stops if it differs from
  `DESCRIPTION` or if the checkout has uncommitted changes (`git status --porcelain`, untracked files included),
  so every results folder maps to one commit. Development mode (`run --dev`, `BULKRNASEQ_DEV=1`) turns both
  into warnings and is recorded as `dev_mode` in `toolkit.json` and in the run manifest.
- Releases: bump `Version` in `DESCRIPTION`, add a CHANGELOG entry, tag `v<version>`. A project moves to a new
  release by changing `toolkit.version` and re-running; results of an older release are reproduced by checking
  out its tag.
- `config_resolved.yaml` and `toolkit.json` are written at parse time only when their content changes
  (`write_if_changed`), and every rule depends on them. A config change or any change to the toolkit code
  (`code_md5` over `DESCRIPTION`, `NAMESPACE`, `R/`, `inst/`, `workflow/`) or to its commit therefore re-runs the
  whole workflow, independently of Snakemake's metadata in `.snakemake/` (which is lost when the working folder
  changes). A dry run refreshes these two files as well.

## 8. Extending

Adding an adapter or a module = one exported function + one registry entry (`inst/modules.yaml` or the
adapter table in the schema) + one schema fragment + tests on fixtures + a section in `docs/MODULES.md`.
The generic module rule and wrapper need no change. `CONTRIBUTING.md` walks through an example.

## 9. Testing strategy

| Layer | Tests |
|---|---|
| Adapters | One fixture per format generated by `tests/fixtures/make_fixtures.R` (seeded, < 200 KB each): STAR 2- and 4-column with `N_*` rows, StringTie with a duplicated gene ID, Excel with several sheets, a sheet name containing `>` and `<dimension ref="A1"/>`, signed FoldChange, a NOISeq table, a bedGraph pair, a small MultiQC HTML with LZString-compressed data, an md5sum file. |
| Canonical objects | Validators reject each malformed variant with an informative error naming the object and column. |
| Modules | Unit tests per function on synthetic data with known answers (e.g. NNLS on an exactly solvable mixture, Hart's formula against a hand-computed value, concordance on constructed vectors). |
| Guardrails | `assert_df()` errors at df = 0 and warns at df = 1; a module requesting p-values on a 1-vs-1 contrast fails with the documented message. |
| Capabilities | Resolution on the three scenarios: incompatible modules are skipped with the expected reason. |
| Integration | `snakemake` on fixtures: (a) STAR + StringTie + bedGraph + MultiQC; (b) raw and normalised counts from Excel; (c) DE tables only (DESeq2 and NOISeq). Two runs give identical md5 for every output except `run_timestamps.json`. |
| Equivalence | Outside this repository: a project instance compares its outputs with the reference tables of the analysis being replaced, column by column (relative tolerance 1e-8, absolute 1e-12 near 0). |
| Repository hygiene | `lintr`, `styler`, `R CMD check` (0 errors, 0 warnings), `snakemake --lint`, `ruff`, `pytest`, and a data guard: no absolute paths, no file > 500 KB, no data extensions outside `tests/fixtures/`. |

Coverage target: ≥ 80% for adapters and statistical modules (`covr`).

## 10. Out of scope, and why

- **FASTQ processing, alignment and quantification** (including FASTQ checksums): the toolkit starts from provider
  outputs. `docs/FASTQ_BOUNDARY.md` specifies the files an upstream step must deliver. Light FASTQ metrics
  computed elsewhere can be imported as a QC table.
- **Project narrative steps** (checking numbers quoted in a report, bibliography retrieval, process audits,
  parsing a provider's PDF, decision rules written for one study): they encode one project's narrative,
  not a reusable computation.
- **Organisms other than human**: the table is extensible, but only human is implemented and tested.
- **Interactive reports** (HTML dashboards): `report/summary.md` and Excel exports only.
