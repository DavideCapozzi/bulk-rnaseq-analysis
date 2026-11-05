# Configuration

<!-- Generated from inst/modules.yaml and inst/schema/schema.yaml by `Rscript tools/make_docs.R`; do not edit by hand (a test fails when this file is out of date). -->

A project is one YAML file validated against `inst/schema/schema.yaml` (`config/schema.yaml` links to it). Defaults are filled from the schema and recorded in `00_provenance/config_resolved.yaml`; relative paths resolve against the folder of the config file. `bulkrnaseq init` writes a skeleton and `config/examples/` holds one complete config per input scenario.

Sections: [`project`](#project), [`toolkit`](#toolkit), [`samples`](#samples), [`inputs`](#inputs), [`resources`](#resources), [`annotation`](#annotation), [`gene_filter`](#gene_filter), [`contrasts`](#contrasts), [`contrast_pairs`](#contrast_pairs), [`designs`](#designs), [`thresholds`](#thresholds), [`modules`](#modules), [`export`](#export). Parameters of the analysis modules are in [MODULES.md](MODULES.md).

Schema types used below:

- `file_or_pattern`: Path relative to inputs.base_dir. May contain the placeholder {input}, replaced by each value of the sample sheet column input_column (one file per sample).
- `gene_list_ref`: Key of a gene list in resources.gene_lists, or an inline list of gene symbols.
- `seed`: Non-negative integer seed of a random number generator.
- `threshold_list`: Non-negative thresholds (e.g. \|log2 fold change\| cut-offs).

## project

Project identity, output folder and default seed. Required.

Required fields: `name`, `output_dir`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `project.name` | string |  | Short project identifier (snake_case). |
| `project.organism` | string (human) | `"human"` | Key of inst/organisms.yaml. |
| `project.output_dir` | string |  | Output folder; relative paths resolve against the config file. |
| `project.seed` | seed | `1` | Default seed for every random step. |

## toolkit

Toolkit release this config is written for. The workflow stops before running when it differs from the Version field of the checked-out DESCRIPTION or when the checkout has uncommitted changes; development mode (BULKRNASEQ_DEV=1, bulkrnaseq run --dev) turns both into warnings recorded in the run manifest. Required.

Required fields: `version`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `toolkit.version` | string |  | Toolkit version, e.g. 0.1.0. |

## samples

Sample sheet and experimental factors.

Required fields: `sheet`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `samples.sheet` | string |  | TSV with at least sample_id and input_column. sample_id is the identifier used in every output; input_column links the sample to a file ({input} placeholder) or to a column of a matrix input. Other columns are sample metadata; metadata is never parsed from names. Optional column inputs_available (comma-separated data levels, default: all configured) marks samples present only in some inputs, e.g. 'qc' for libraries with QC reports but no counts. |
| `samples.factors` | object | `{}` | Sample sheet columns used as experimental factors, with level order (reference first). |
| `samples.label_column` | string |  | Optional sample sheet column used as plot label. |

## inputs

Provider files. Every data level is optional; modules that need a missing level are skipped with the reason. Required.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `inputs.base_dir` | string |  | Folder all input paths are relative to. |
| `inputs.counts` | object |  | Raw read counts per gene and sample (integers). |
| `inputs.counts.adapter` | string (star_reads_per_gene, matrix) |  | star_reads_per_gene: one STAR ReadsPerGene.out.tab per sample (files); matrix: one table, genes x samples (file, id_column). |
| `inputs.counts.files` | file_or_pattern |  | STAR files, one per sample (star_reads_per_gene). |
| `inputs.counts.file` | string |  | TSV or Excel table, genes x samples (matrix); value columns are matched to the sample sheet's input_column. |
| `inputs.counts.sheet` | string or null | `null` | Excel sheet (matrix; null: the first sheet). |
| `inputs.counts.id_column` | string |  | Column with the gene identifiers (matrix). |
| `inputs.counts.strand_column` | string (unstranded, forward, reverse, single) | `"single"` | STAR ReadsPerGene.out.tab has 4 columns (unstranded, forward, reverse); 'single' is for the 2-column 'deseq2-ready' files where the provider already chose the strand. |
| `inputs.counts.on_duplicate` | string (error, sum, first) | `"error"` | Rows with the same gene identifier: error stops, sum adds them, first keeps the first. |
| `inputs.abundance` | object |  | Gene abundances (TPM, FPKM) per sample. |
| `inputs.abundance.adapter` | string (stringtie_gene_abundance, matrix) |  | stringtie_gene_abundance: one StringTie gene abundance table per sample (files); matrix: one table, genes x samples (file, id_column, kind). |
| `inputs.abundance.files` | file_or_pattern |  | StringTie gene abundance files, one per sample (stringtie_gene_abundance). |
| `inputs.abundance.file` | string |  | TSV or Excel table, genes x samples (matrix); value columns are matched to the sample sheet's input_column. |
| `inputs.abundance.sheet` | string or null | `null` | Excel sheet (matrix; null: the first sheet). |
| `inputs.abundance.id_column` | string |  | Column with the gene identifiers (matrix). |
| `inputs.abundance.kind` | string (tpm, fpkm) |  | For the matrix adapter only. |
| `inputs.abundance.on_duplicate` | string (error, sum, first) | `"error"` | Rows with the same gene identifier: error stops, sum adds them, first keeps the first. |
| `inputs.normalized` | object |  | Provider-normalised expression, kept as an assay and checked against the raw counts in 01_ingest (normalization_check). |
| `inputs.normalized.adapter` | string (matrix) |  | matrix: one table, genes x samples. |
| `inputs.normalized.file` | string |  | TSV or Excel table, genes x samples; value columns are matched to the sample sheet's input_column. |
| `inputs.normalized.sheet` | string or null | `null` | Excel sheet (null: the first sheet). |
| `inputs.normalized.id_column` | string |  | Column with the gene identifiers. |
| `inputs.normalized.description` | string |  | How the provider normalised (free text, copied to the ingest report). |
| `inputs.contrast_tables` | array of object |  | Provider differential-expression tables, imported as canonical contrast tables. |
| `inputs.contrast_tables[].id` | string |  | Contrast id used by the modules and in every output. |
| `inputs.contrast_tables[].adapter` | string (deseq2, noiseq, generic) |  | deseq2 and noiseq: the tools' column names; generic: columns given in columns. |
| `inputs.contrast_tables[].file` | string |  | TSV or Excel file. |
| `inputs.contrast_tables[].sheet` | string or null | `null` | Excel sheet (null: the first sheet). |
| `inputs.contrast_tables[].columns` | object |  | Mapping canonical column -> input column; required for 'generic'. |
| `inputs.contrast_tables[].numerator` | string |  | Numerator group of the table, for the reader (not used by any step). |
| `inputs.contrast_tables[].denominator` | string |  | Denominator group of the table, for the reader (not used by any step). |
| `inputs.contrast_tables[].df_residual` | integer or null | `null` | Residual degrees of freedom of the provider's model (n samples - rank of the design; 0 without replicates). Required when the table has padj or posterior_prob: it decides whether its p-values are sample-level. |
| `inputs.coverage` | object |  | Input of the bedGraph coverage adapter, not implemented yet; validate_config() rejects it until it is ported. |
| `inputs.qc` | object |  | Imported QC reports. |
| `inputs.qc.multiqc` | array of object |  | MultiQC HTML reports, decoded into 01_ingest/multiqc/ and joined to the QC summary. |
| `inputs.qc.multiqc[].id` | string |  | Report id (names the decoded files in 01_ingest/multiqc/). |
| `inputs.qc.multiqc[].file` | string |  | MultiQC HTML report with embedded compressed plot data (MultiQC 1.x; tested with 1.12 and 1.13). |
| `inputs.qc.multiqc[].read` | string (r1, r2) |  | The report holds FastQC on raw reads of this mate: its FastQC metrics get the suffix _r1/_r2 and r1 gives raw_read_pairs. FastQC sections of reports without `read` are not read. |
| `inputs.qc.sample_aliases` | object |  | MultiQC sample name -> input_column, for names that do not contain the input_column (by default a MultiQC name belongs to the one sample whose input_column it contains). |
| `inputs.checksums` | object |  | Input of the md5sum checksum file adapter, not implemented yet; validate_config() rejects it until it is ported. |

## resources

Reference files (GTF, gene sets, HPA exports, KEGG) and gene lists.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `resources.cache_dir` | string | `"resources"` | Folder of the retrieved resources; a relative path resolves against the config file. |
| `resources.items` | object |  | Resources fetched by the fetch_resources rule into cache_dir. 'source' is a URL or file:// path; 'md5' (optional) is checked after retrieval. Analysis rules only read the cache. |
| `resources.items.<key>.source` | string |  | URL, file:// path, or a retrieval keyword: kegg_rest (KEGG pathways of the organism) or msigdbr (MSigDB gene sets of the organism from the msigdbr package, restricted to collections). |
| `resources.items.<key>.collections` | array of string |  | For source msigdbr: collections to export, each "collection" or "collection:subcollection" (e.g. H, C2:CP:REACTOME, C7:IMMUNESIGDB). |
| `resources.items.<key>.md5` | string |  | Expected md5, checked after retrieval and before use. |
| `resources.items.<key>.version` | string |  | Release recorded in the resources manifest (e.g. Ensembl 115). |
| `resources.items.<key>.kind` | string (gtf, msigdb_export, hpa_immune_cell, hpa_single_cell, fasta, kegg_pathways, tsv) |  | File type; decides the reader. |
| `resources.gene_lists` | array of string | `[]` | Named gene lists (TSV with columns list, gene_symbol, source). Package defaults are in inst/extdata/gene_lists/. |

## annotation

Gene annotation source and identifier type of the inputs.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `annotation.gtf` | string | `"gtf"` | Key of the GTF in resources.items. |
| `annotation.id_type` | string (ensembl_gene_id, gene_symbol, entrez_id) | `"ensembl_gene_id"` | Type of the gene identifiers in the input files. |
| `annotation.symbol_source` | string (annotation, abundance) | `"annotation"` | Where gene_symbol comes from when both are present. |

## gene_filter

Default gene filter for count-based modules (modules may override).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gene_filter.biotypes` | array of string | `["protein_coding","lncRNA"]` | Gene biotypes kept (null: all). |
| `gene_filter.exclude_mitochondrial` | boolean | `true` | Drop the genes of the mitochondrial chromosome. |
| `gene_filter.min_cpm` | number | `2` | CPM (library size only) that a gene needs in at least min_samples samples. |
| `gene_filter.min_samples` | integer | `2` | Number of samples in which a gene needs min_cpm. |
| `gene_filter.sex_genes` | string (keep, drop_xist_y, drop_xy) | `"drop_xist_y"` | keep: all genes; drop_xist_y: drop chrY genes and the gene list x_inactivation; drop_xy: drop chrX and chrY genes. |

## contrasts

Contrasts computed from counts.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `contrasts[].id` | string |  | Contrast id used in every output. |
| `contrasts[].numerator` | array of string |  | sample_ids, or "factor=level". |
| `contrasts[].denominator` | array of string |  | sample_ids, or "factor=level". |
| `contrasts[].label` | string |  | Free-text label shown in reports. |

## contrast_pairs

Pairs of contrasts expected to share an effect (used by concordance and gene-set concordance).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `contrast_pairs[].id` | string |  | Pair id used in outputs. |
| `contrast_pairs[].first` | string |  | Contrast id (keys x/y are avoided because YAML 1.1 reads 'y' as a boolean). |
| `contrast_pairs[].second` | string |  | Contrast id. |
| `contrast_pairs[].label` | string |  | Free-text label shown in reports. |

## designs

Model designs for inferential DE (module de_test) and dispersion estimates.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `designs[].id` | string |  | Design id referenced by modules (de_test, contrasts.dispersion_design, power.design). |
| `designs[].formula` | string |  | e.g. "~ patient + condition" |
| `designs[].contrast` | array of string |  | [factor, numerator level, denominator level] tested by de_test; absent: the design only provides dispersions. |
| `designs[].engine` | string (deseq2) | `"deseq2"` | Fitting engine (DESeq2 only). |

## thresholds

Significance and fold-change thresholds of de_test, volcano and ORA.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `thresholds.padj` | number | `0.05` | Adjusted p-value threshold (tables whose significance_metric is padj). |
| `thresholds.abs_fold_change` | number | `1.5` | Linear \|fold change\|; log2 threshold is derived. |
| `thresholds.posterior_prob` | number | `0.8` | Posterior probability threshold (tables whose significance_metric is posterior_prob, e.g. NOISeq). |

## export

Extra export formats.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `export.xlsx` | boolean | `false` | One Excel workbook per module in report/, one sheet per table. |
| `export.methods` | boolean | `true` | Write report/methods.md, a methods text of the modules that ran, with versions and references. |
| `export.methods_notes` | string or null | `null` | Free text added to the methods text under Notes (e.g. what the provider did before delivery). |
| `export.researcher_workbook` | object | `{}` | One workbook in report/ with the main results for the researcher (see the README sheet). |
| `export.researcher_workbook.enabled` | boolean | `false` | Write the researcher workbook. |
| `export.researcher_workbook.file` | string | `"results_workbook.xlsx"` | File name in report/. |
| `export.researcher_workbook.sheets` | array of string (de, ora, geneset_test, markers, signatures, qc, pca) | `["de","ora","geneset_test","markers","signatures","qc","pca"]` | Result types included, in this order, when their module ran: de (per contrast with a significance metric: all genes, UP, DOWN), ora, geneset_test, markers, signatures, qc, pca. |

## modules

Parameters of the analysis modules: one section per module, each with `enabled` (default `true`), documented in [MODULES.md](MODULES.md): [`qc`](MODULES.md#02-qc), [`pca`](MODULES.md#03-pca), [`markers`](MODULES.md#04-markers), [`contamination`](MODULES.md#05-contamination), [`contamination_sets`](MODULES.md#06-contamination_sets), [`reference`](MODULES.md#07-reference), [`deconvolution`](MODULES.md#08-deconvolution), [`contrasts`](MODULES.md#09-contrasts), [`de_test`](MODULES.md#10-de_test), [`concordance`](MODULES.md#11-concordance), [`geneset_test`](MODULES.md#12-geneset_test), [`ora`](MODULES.md#13-ora), [`volcano`](MODULES.md#14-volcano), [`signatures`](MODULES.md#15-signatures), [`pairing`](MODULES.md#16-pairing), [`power`](MODULES.md#17-power), [`coverage`](MODULES.md#18-coverage), [`covariates`](MODULES.md#19-covariates), [`reference_axis`](MODULES.md#20-reference_axis), [`de_heatmap`](MODULES.md#21-de_heatmap).
