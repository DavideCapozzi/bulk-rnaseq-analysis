# Modules

<!-- Generated from inst/modules.yaml and inst/schema/schema.yaml by `Rscript tools/make_docs.R`; do not edit by hand (a test fails when this file is out of date). -->

Every module declares the data levels it needs. A module whose needs are not met, or that is not implemented yet, is skipped with the reason in `00_provenance/capabilities.json` and `report/summary.md`. Defaults come from the config schema; `bulkrnaseq plan <config>` shows what a project runs. The other sections of the config are in [CONFIG.md](CONFIG.md).

20 modules: 16 implemented, 4 not implemented.

| NN | Module | Status | Data levels | Provides |
|---|---|---|---|---|
| 02 | [`qc`](#02-qc) | implemented | one of `counts`, `qc` | none |
| 03 | [`pca`](#03-pca) | implemented | `counts` | none |
| 04 | [`markers`](#04-markers) | implemented | `counts` | none |
| 05 | [`contamination`](#05-contamination) | implemented | `counts` | `gene_mask` |
| 06 | [`contamination_sets`](#06-contamination_sets) | not implemented | `ranked` | none |
| 07 | [`reference`](#07-reference) | implemented | `counts` | `gene_universe` |
| 08 | [`deconvolution`](#08-deconvolution) | not implemented | one of `abundance`, `counts` | none |
| 09 | [`contrasts`](#09-contrasts) | implemented | `counts` | `contrast_stats`, `ranked` |
| 10 | [`de_test`](#10-de_test) | implemented | `counts` | `contrast_stats`, `ranked` |
| 11 | [`concordance`](#11-concordance) | implemented | `contrast_stats` | none |
| 12 | [`geneset_test`](#12-geneset_test) | implemented | `ranked` | none |
| 13 | [`ora`](#13-ora) | implemented | `contrast_stats` | none |
| 14 | [`volcano`](#14-volcano) | implemented | `contrast_stats` | none |
| 15 | [`signatures`](#15-signatures) | implemented | `counts` | none |
| 16 | [`pairing`](#16-pairing) | implemented | `counts` | none |
| 17 | [`power`](#17-power) | implemented | `counts` | none |
| 18 | [`coverage`](#18-coverage) | not implemented | `coverage` | none |
| 19 | [`covariates`](#19-covariates) | not implemented | `contrast_stats` | none |
| 20 | [`reference_axis`](#20-reference_axis) | implemented | `contrast_stats` | none |
| 21 | [`de_heatmap`](#21-de_heatmap) | implemented | `contrast_stats`, `counts` | none |

Schema types used below:

- `file_or_pattern`: Path relative to inputs.base_dir. May contain the placeholder {input}, replaced by each value of the sample sheet column input_column (one file per sample).
- `gene_list_ref`: Key of a gene list in resources.gene_lists, or an inline list of gene symbols.
- `seed`: Non-negative integer seed of a random number generator.
- `threshold_list`: Non-negative thresholds (e.g. \|log2 fold change\| cut-offs).

## 02 qc

Status: implemented (`run_qc`). Output folder: `02_qc/`.

- Data levels: one of `counts`, `qc`
- Optional levels: `abundance`
- Resources: `annotation`
- Conditions: none
- Provides: none
- Tables: `qc_summary`, `expression_qc`, `top_genes`
- Figures: `qc_metrics`

Caveats:

- Expression QC uses assigned reads per sample; percentages are of the library size, not of raw reads.
- MultiQC metrics are taken or derived from the provider's MultiQC reports (STAR, RSeQC, FastQC on raw reads); samples with inputs_available: qc appear only in qc_summary.
- assigned_reverse and assigned_same are STAR GeneCounts 'Reverse Stranded' and 'Same Stranded' columns whatever the library protocol; strand_ratio_star = same / reverse.

Parameters (`modules.qc`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `top_n_share` | integer | `50` | N of pct_counts_top<N>_genes, the share of reads in the N most expressed genes. |
| `top_n_list` | integer | `15` | Most expressed genes listed per sample in top_genes. |
| `ribosomal_protein_regex` | string | `"^RP[SL][0-9]+[A-Z]?[0-9]*$"` | Regular expression on the symbols of protein-coding ribosomal-protein genes. |
| `hemoglobin_genes` | gene_list_ref | `"hemoglobin"` | Hemoglobin genes. |
| `rrna_biotypes` | array of string | `["rRNA","Mt_rRNA","rRNA_pseudogene"]` | Biotypes counted as rRNA. |
| `figure_metrics` | array of string | `["library_size","genes_count_ge10","pct_counts_mito","pct_counts_rrna","pct_counts_protein_coding","pct_counts_top50_genes"]` | Metrics shown in the qc_metrics figure (absent ones are skipped). |

## 03 pca

Status: implemented (`run_pca`). Output folder: `03_pca/`.

- Data levels: `counts`
- Optional levels: `abundance`, `qc`, `gene_mask`, `gene_universe`
- Resources: `annotation`
- Conditions: `min_samples: 3`
- Provides: none
- Tables: `pca_summary`, `pca_scores`, `pca_scores_log2tpm`, `sample_correlations`, `sample_distances`, `pc1_vs_qc`, `pc1_loading_vs_gene_features`, `pc1_top_loadings`
- Figures: `pca`, `sample_correlation`, `pc1_loading_vs_gene_span`

Caveats:

- With n samples the PCA has at most n - 1 components; correlations with QC metrics are descriptive.
- pc1_vs_qc lists candidates, not drivers: with n samples even a perfect rank correlation has an exact two-sided permutation p of 2/n! (1/12 with 4 samples), and with 4 samples |Spearman rho| >= 0.8 occurs by chance 1 time in 3.
- 'separates' labels the factor with the largest mean difference on a PC; it is not a test.

Parameters (`modules.pca`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `transforms` | array of string (vst, logcpm_tmm) | `["vst","logcpm_tmm"]` | Transforms of the grid: vst (DESeq2) and logcpm_tmm (log2 CPM after TMM). |
| `ntop` | array of integer or string | `[500,1000,2000,"all"]` | Numbers of most variable genes of the grid (all: every filtered gene). |
| `sex_gene_variants` | array of string | `["keep","drop_xist_y","drop_xy"]` | gene_filter.sex_genes values of the grid. |
| `gene_subsets` | array of string | `"all"` | Gene subset expressions of the grid. |
| `symbol_exclusion_variants` | array of string | `[]` | Regexes on gene symbols (e.g. "^HLA-"); each adds a variant of the grid without the matching genes, besides the variant without exclusion. |
| `thinning` | boolean | `false` | Adds a variant on counts thinned to the smallest library (edgeR::thinCounts, project.seed). |
| `separation_factors` | array of string | `[]` | Factors labelled in pca_summary (default: all samples.factors, in order). |
| `orient_pc1` | object or null | `null` | Orients PC1 so that samples at this level have positive mean score (drivers, pc1_order). |
| `orient_pc1.factor` | string |  | Sample sheet factor. |
| `orient_pc1.level` | string |  | Level whose samples get a positive mean PC1 score. |
| `abundance_sex_genes` | string (keep, drop_xist_y, drop_xy) | `"keep"` | sex_genes rule of the PCA on log2(TPM + 1). |
| `scores_ntop` | integer or string | `1000` | Number of most variable genes of pca_scores (a value of ntop). |
| `pc_drivers` | boolean | `true` | Report PC1 drivers on VST (correlation with QC metrics, loadings against gene features). |
| `qc_metrics_for_drivers` | array of string | `[]` | QC metrics correlated with the PCs in pc1_vs_qc (absent ones are skipped; empty: none). |
| `top_loadings` | integer | `40` | Genes with the largest absolute PC1 loading in pc1_top_loadings. |

## 04 markers

Status: implemented (`run_markers`). Output folder: `04_markers/`.

- Data levels: `counts`
- Optional levels: `abundance`, `normalized`
- Resources: none
- Conditions: `configured: panels`
- Provides: none
- Tables: `marker_expression`, `marker_panel_scores`, `marker_ratios`
- Figures: `markers_heatmap`, `marker_panels`

Caveats:

- Panels are starting markers, not validated signatures. CPM are library-size only (no TMM) unless stated.

Parameters (`modules.markers`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `panels` | array of object | `[]` | Marker panels, shown in this order; the module runs only when at least one is configured. |
| `panels[].id` | string |  | Panel id used in outputs. |
| `panels[].label` | string |  | Free-text label shown in reports. |
| `panels[].genes` | gene_list_ref |  | Genes of the panel (a gene list name or inline symbols). |
| `biotypes` | array or null | `null` | A panel symbol matches the features with that symbol whose biotype is listed (e.g. [protein_coding]); null: any biotype. |
| `ratios` | array of object | `[]` | log2((CPM numerator + pseudocount) / (CPM denominator + pseudocount)) of a gene in a sample. |
| `ratios[].id` | string |  | Ratio id used in outputs. |
| `ratios[].label` | string |  | Free-text label shown in reports. |
| `ratios[].numerator` | object |  | Gene symbol (one feature after the biotype filter) and sample_id of the numerator. |
| `ratios[].numerator.gene` | string |  | Gene symbol. |
| `ratios[].numerator.sample` | string |  | sample_id. |
| `ratios[].denominator` | object |  | Gene symbol and sample_id of the denominator. |
| `ratios[].denominator.gene` | string |  | Gene symbol. |
| `ratios[].denominator.sample` | string |  | sample_id. |
| `pseudocount` | number | `1` | Added to both CPM values of each ratio. |

## 05 contamination

Status: implemented (`run_contamination`). Output folder: `05_contamination/`.

- Data levels: `counts`
- Optional levels: none
- Resources: `annotation`
- Conditions: `configured: profiles, mask.include`
- Provides: `gene_mask`
- Tables: `contamination_scores`, `contamination_profile_genes`, `profile_ratio_consistency`, `profile_excess`, `contaminant_mask`
- Figures: `contamination_profiles`

Caveats:

- Profile scores, ratios and excesses are descriptive; ratio consistency across a profile's genes is compatible with, not proof of, a single source.
- The mask removes genes that the configured gene sets (e.g. a lung cell atlas) assign to contaminating cell types, plus the configured profiles; it cannot remove contaminant signal from other genes, nor show that the remaining signal is intrinsic to the target cells.
- Profiles chosen from the same samples (e.g. PC1 drivers) make the mask partly data-derived: say so where it is used.

Parameters (`modules.contamination`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `profiles` | array of object | `[]` | Gene profiles of contaminating (or reference) cells, scored per sample on CPM of all genes. |
| `profiles[].id` | string |  | Profile id. |
| `profiles[].genes` | gene_list_ref |  | Genes of the profile. |
| `reference_profile` | string or null | `null` | Profile id whose CPM sum is also reported relative to the sample with the largest value (rel_to_max). |
| `ratio_pseudocount` | number | `1` | CPM added before the log2 ratios of ratio_checks. |
| `ratio_checks` | array of object | `[]` | Per-gene log2 CPM ratios between two samples over a profile's genes. |
| `ratio_checks[].id` | string |  | Check id. |
| `ratio_checks[].profile` | string |  | Profile id whose genes are compared. |
| `ratio_checks[].numerator_sample` | string |  | sample_id of the numerator. |
| `ratio_checks[].denominator_sample` | string |  | sample_id of the denominator. |
| `excess_checks` | array of object | `[]` | CPM of genes above the mean of baseline_samples, divided by the CPM sum of a profile, in two samples, and the ratio of the two excesses (a similar ratio for genes and profile is compatible with the genes coming from the profile's cells). |
| `excess_checks[].id` | string |  | Check id. |
| `excess_checks[].genes` | gene_list_ref |  | Genes whose excess is measured. |
| `excess_checks[].profile` | string |  | Profile id whose CPM sum divides the excess. |
| `excess_checks[].baseline_samples` | array of string |  | sample_ids whose mean CPM is the baseline. |
| `excess_checks[].samples` | array of string |  | The two sample_ids compared. |
| `mask` | object | `{}` | masked = (in any include group and in no exclude group) or among the genes of add_profiles. Groups are gene sets of an MSigDB export (resource); their genes match by Ensembl ID (by symbol when annotation.id_type is gene_symbol). |
| `mask.resource` | string or null | `null` | Key of an msigdb_export in resources.items. |
| `mask.include` | array of object | `[]` | Groups of gene sets whose genes are masked. |
| `mask.include[].id` | string |  | Group id (gene flag column of contaminant_mask). |
| `mask.include[].sets` | array of string |  | Gene set names (gs_name) of the MSigDB export. |
| `mask.exclude` | array of object | `[]` | Groups of gene sets whose genes are not masked, even when in an include group. |
| `mask.exclude[].id` | string |  | Group id (gene flag column of contaminant_mask). |
| `mask.exclude[].sets` | array of string |  | Gene set names (gs_name) of the MSigDB export. |
| `mask.add_profiles` | array of string | `[]` | Profile ids whose genes are always masked. |

## 06 contamination_sets

Status: not implemented (`run_contamination_sets`). Output folder: `06_contamination_sets/`.

- Data levels: `ranked`
- Optional levels: none
- Resources: `contaminant_profiles`
- Conditions: none
- Provides: none
- Tables: `contamination_set_tests`
- Figures: none

Caveats:

- Gene-set test of contaminant profiles on one contrast; p-values follow the contrast's inference_scope.

Parameters (`modules.contamination_sets`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |

## 07 reference

Status: implemented (`run_reference`). Output folder: `07_reference/`.

- Data levels: `counts`
- Optional levels: `gene_mask`
- Resources: `annotation`
- Conditions: `configured: profiles`
- Provides: `gene_universe`
- Tables: `reference_profiles`, `reference_universe`, `universe_exclusions`
- Figures: none

Caveats:

- Reference profiles come from other platforms, tissues and units (nTPM, nCPM): thresholds are coarse filters.
- reference_universe uses the references only; combine it with the contamination mask explicitly (in_reference_universe+no_masked).

Parameters (`modules.reference`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `profiles` | array of object | `[]` | Reference expression of one cell type from an HPA export: nTPM for kind hpa_immune_cell, nCPM for hpa_single_cell (kind of the resource in resources.items). |
| `profiles[].id` | string |  | Profile id. |
| `profiles[].resource` | string |  | Key in resources.items. |
| `profiles[].cell` | string |  | Value of the 'Immune cell' or 'Cell type' column. |
| `target` | string or null | `null` | Profile id of the target cells. |
| `contaminants` | array of string | `[]` | Profile ids of contaminating cells. |
| `min_target_expr` | number | `1` | target_expressed: target value >= this. |
| `contam_min` | number | `20` | contam_sensitive: max contaminant value >= this ... |
| `contam_fold` | number | `5` | ... and >= contam_fold x (target + 1). |
| `exclusion_candidates` | gene_list_ref | `[]` | Genes (protein-coding) whose inclusion in the universe, or the reason for exclusion, is reported. |

## 08 deconvolution

Status: not implemented (`run_deconvolution`). Output folder: `08_deconvolution/`.

- Data levels: one of `abundance`, `counts`
- Optional levels: none
- Resources: `reference_profiles`
- Conditions: none
- Provides: none
- Tables: `deconvolution`, `deconvolution_variants`
- Figures: none

Caveats:

- Relative proportions between reference components; weights do not sum to 1 across platforms.
- Report the range across variants, not a point estimate.

Parameters (`modules.deconvolution`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |

## 09 contrasts

Status: implemented (`run_contrasts`). Output folder: `09_contrasts/`.

- Data levels: `counts`
- Optional levels: `gene_mask`, `gene_universe`
- Resources: `annotation`
- Conditions: `contrasts_defined: true`
- Provides: `contrast_stats`, `ranked`
- Tables: `contrast_stats`, `contrast_spread`, `normalization_factors`
- Figures: none

Caveats:

- Log-CPM differences are descriptive estimators; no p-values are produced (assert_df not involved).
- Shrunken estimators fitted with 1 residual df are flagged exploratory and can be non-informative.

Parameters (`modules.contrasts`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `estimators` | array of string (logcpm_difference, deseq2_apeglm) | `"logcpm_difference"` | logcpm_difference: difference of mean log2 CPM (any contrast, no p-value); deseq2_apeglm: apeglm-shrunken DESeq2 fold changes of one-vs-one contrasts (needs dispersion_design). |
| `prior_count` | number | `2` | Prior count of the log2 CPM of logcpm_difference. |
| `norm_method` | string (TMM, RLE, upperquartile, none, reference) | `"TMM"` | Library normalisation of logcpm_difference (reference: factors computed on reference_genes). |
| `reference_genes` | null or | `null` | Genes of norm_method reference. |
| `dispersion_design` | string or null | `null` | Id of the designs[] entry whose dispersions deseq2_apeglm uses (at least 1 residual degree of freedom). |
| `spread_subsets` | array of string | `"all"` | Gene subsets of contrast_spread. |
| `spread_threshold` | number | `1` | \|log2 fold change\| counted in contrast_spread (n_abs_lfc_gt_threshold). |

## 10 de_test

Status: implemented (`run_de_test`). Output folder: `10_de_test/`.

- Data levels: `counts`
- Optional levels: none
- Resources: `annotation`
- Conditions: `design_df_min: 1`
- Provides: `contrast_stats`, `ranked`
- Tables: `contrast_stats`, `design_summary`, `provider_comparison`
- Figures: none

Caveats:

- assert_df(): error at 0 residual df, warning and exploratory flag at 1.
- Fold changes are maximum-likelihood estimates (not shrunken); p-values have inference_scope sample.

Parameters (`modules.de_test`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `prefilter` | object | `{"min_count":10,"min_samples":3}` | Genes with at least min_count reads in at least min_samples samples are tested. |
| `prefilter.min_count` | integer |  | Reads a gene needs in at least min_samples samples. |
| `prefilter.min_samples` | integer |  | Number of samples in which a gene needs min_count reads. |
| `results_alpha` | number | `0.1` | alpha of DESeq2::results(), the FDR target of independent filtering (DESeq2 default 0.1). |
| `compare_with` | array of object | `[]` | Computed design vs imported contrast table, reported in provider_comparison. |
| `compare_with[].design` | string |  | Design id. |
| `compare_with[].contrast_table` | string |  | Id of an inputs.contrast_tables entry. |

## 11 concordance

Status: implemented (`run_concordance`). Output folder: `11_concordance/`.

- Data levels: `contrast_stats`
- Optional levels: `counts`, `gene_mask`, `gene_universe`
- Resources: none
- Conditions: `min_contrasts: 2`, `contrast_pairs_defined: true`
- Provides: none
- Tables: `concordance_summary`
- Figures: `concordance_scatter`

Caveats:

- Sign concordance between two contrasts measures shared effects; its power depends on the expression stratum.
- Contrasts that share samples (shared_samples > 0) share their noise, which inflates concordance.

Parameters (`modules.concordance`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `thresholds` | threshold_list | `[0.58,1]` | \|log2 fold change\| thresholds of the concordance counts. |
| `subsets` | array of string | `"all"` | Gene subset expressions summarised. |

## 12 geneset_test

Status: implemented (`run_geneset_test`). Output folder: `12_geneset_test/`.

- Data levels: `ranked`
- Optional levels: `gene_mask`, `gene_universe`
- Resources: none
- Conditions: none
- Provides: none
- Tables: `geneset_results`, `geneset_method_comparison`, `geneset_concordance`, `geneset_concordant_sets`
- Figures: `nes_concordance`, `nes_heatmap`

Caveats:

- fgsea pre-ranked permutes genes, not samples: its p-values only rank sets within one contrast (inference_scope ranking_only).
- cameraPR accounts for inter-gene correlation but still tests genes for a given contrast, not samples: with one library per condition it is conditional_on_contrast.

Parameters (`modules.geneset_test`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `resource` | string | `"gene_sets"` | Key of an msigdb_export in resources.items. |
| `estimator` | string or null | `null` | Method of the contrast tables whose statistic ranks the genes (e.g. logcpm_difference); null: the first method available. |
| `ranking_statistic` | string (lfc, stat) | `"lfc"` | Column that ranks the genes. lfc: log2 fold change; stat: test statistic of DESeq2 tables (Wald stat, imported or from de_test), which also weighs the precision of each fold change. |
| `contrasts` | array of string | `[]` | Contrast ids to test, in this order (empty = every contrast of the estimator, config order). |
| `methods` | array of string (camera_pr, fgsea) | `"camera_pr"` | camera_pr: limma::cameraPR; fgsea: pre-ranked fgsea (multilevel). |
| `collections` | array of object | `[{"id":"H","collection":"H"}]` | Gene set collections of the export: gs_collection, optional gs_subcollection and a regular expression on gs_name. P-values are adjusted over all collections together. |
| `collections[].id` | string |  | Collection id used in outputs. |
| `collections[].collection` | string |  | gs_collection of the export (e.g. H, C2, C7). |
| `collections[].subcollection` | string or null | `null` | gs_subcollection (e.g. CP:REACTOME; null: any). |
| `collections[].set_pattern` | string or null | `null` | Regular expression on gs_name (null: every set). |
| `universes` | array of object | `[{"id":"all","subset":"all"}]` | Ranked genes, as subset expressions over the genes of the contrast tables; tested in this order. |
| `universes[].id` | string |  | Universe id used in outputs. |
| `universes[].subset` | string |  | Subset expression over the gene flags (e.g. all, no_masked). |
| `min_size` | integer | `15` | Smallest set size within the universe. |
| `max_size` | integer | `500` | Largest set size within the universe. |
| `inter_gene_cor` | array of number | `0.01` | cameraPR variants. |
| `use_ranks` | array of boolean | `false` | cameraPR variants. |
| `fgsea_tie_jitter_sd` | number | `1e-06` | SD of the normal noise added to the ranking statistic before fgsea, to break ties reproducibly (0: none). Genes are ordered by decreasing expression, then feature_id, before the draw. |
| `seed` | integer or null | `null` | Seed set once before fgsea (universes outer, contrasts inner, on one random stream); null: project.seed. |
| `leading_edge_max` | integer | `15` | Leading-edge symbols reported per set. |
| `alpha` | number | `0.05` | Adjusted p-value threshold of the counts, the concordance and the figures. |
| `concordance_collections` | array of string | `[]` | Collection ids summarised in geneset_concordance (empty = all). |
| `figure_universe` | string or null | `null` | Universe of nes_concordance (null: the first). |
| `figure_collections` | array of string | `[]` | Collection ids drawn in nes_concordance (empty = all). |
| `heatmap_collection` | string or null | `null` | Collection of nes_heatmap (null: the first). |

## 13 ora

Status: implemented (`run_ora`). Output folder: `13_ora/`.

- Data levels: `contrast_stats`
- Optional levels: none
- Resources: `kegg_pathways`
- Conditions: `significance_available: true`
- Provides: none
- Tables: `ora_input_genes`, `ora_go_bp`, `ora_go_cc`, `ora_go_mf`, `ora_kegg`, `ora_reactome`
- Figures: `ora_go_bp`, `ora_go_cc`, `ora_go_mf`, `ora_kegg`, `ora_reactome`

Caveats:

- Hypergeometric p-values treat genes as independent units; the universe is the set of tested genes.

Parameters (`modules.ora`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `databases` | array of string (GO_BP, GO_CC, GO_MF, KEGG, Reactome) | `["GO_BP","GO_CC","GO_MF","KEGG","Reactome"]` | Databases tested (KEGG needs the kegg_pathways resource). |
| `p_adjust_method` | string | `"BH"` | pAdjustMethod of clusterProfiler. |
| `pvalue_cutoff` | number | `0.05` | pvalueCutoff of clusterProfiler (terms are kept when p-value and adjusted p-value are below it). |
| `show_category` | integer | `15` | Terms per contrast in the ORA figures (lowest adjusted p-value first). |

## 14 volcano

Status: implemented (`run_volcano`). Output folder: `14_volcano/`.

- Data levels: `contrast_stats`
- Optional levels: none
- Resources: none
- Conditions: `significance_available: true`
- Provides: none
- Tables: `volcano_points`
- Figures: `volcano`

Caveats:

- The y axis shows the significance metric under its own name (padj or posterior probability).

Parameters (`modules.volcano`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `label_top` | integer | `15` | Significant genes labelled per volcano (most significant first). |

## 15 signatures

Status: implemented (`run_signatures`). Output folder: `15_signatures/`.

- Data levels: `counts`
- Optional levels: `gene_mask`, `gene_universe`
- Resources: none
- Conditions: `configured: msigdb_sets, custom_sets`
- Provides: none
- Tables: `signature_scores`, `signature_differences`, `signature_method_comparison`
- Figures: `signature_scores`

Caveats:

- z-score is relative to the samples in the run; singscore and ssGSEA are within-sample ranks.
- ssGSEA scores are normalised over the sets scored together: adding or removing a set changes the others.
- Differences between contrast sides and their sign agreement within a pair are descriptive (no p-values).

Parameters (`modules.signatures`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `resource` | string or null | `null` | Key of the msigdb_export in resources.items that holds msigdb_sets. |
| `msigdb_sets` | array of string | `[]` | MSigDB set names (gs_name) scored; the set name is its id in outputs. |
| `custom_sets` | array of object | `[]` | Gene sets given as gene symbols (a gene list name or inline symbols). |
| `custom_sets[].id` | string |  | Set id used in outputs. |
| `custom_sets[].genes` | gene_list_ref |  | Genes of the set. |
| `composites` | array of object | `[]` | Two-sided scores: score(up) - score(down) for zscore and ssgsea; simpleScore(up, down) for singscore. |
| `composites[].id` | string |  | Composite id used in outputs. |
| `composites[].up` | string |  | Id of the set expected higher. |
| `composites[].down` | string |  | Id of the set expected lower. |
| `methods` | array of string (zscore, ssgsea, singscore) | `"zscore"` | zscore: mean per-gene z-score across the samples of the run; ssgsea: GSVA ssGSEA; singscore: singscore simpleScore on within-sample ranks. |
| `universes` | array of object | `[{"id":"all","subset":"all"}]` | Genes scored, as subset expressions over the genes kept by gene_filter; scored in this order. |
| `universes[].id` | string |  | Universe id used in outputs. |
| `universes[].subset` | string |  | Subset expression over the gene flags (e.g. all, no_masked). |
| `min_size` | integer | `10` | Sets with fewer genes in a universe are not scored there (also ssGSEA minSize). |
| `prior_count` | number | `2` | Prior count of the log2 CPM. |
| `norm_method` | string (TMM, RLE, upperquartile, none) | `"TMM"` | Library normalisation of the log2 CPM. |
| `figure_universe` | string or null | `null` | Universe of the signature_scores figure (null: the first). |
| `figure_method` | string or null | `null` | Method of the signature_scores figure (null: the first). |

## 16 pairing

Status: implemented (`run_pairing`). Output folder: `16_pairing/`.

- Data levels: `counts`
- Optional levels: `abundance`
- Resources: `annotation`, `sex_genes`, `fingerprint_genes`
- Conditions: none
- Provides: none
- Tables: `sex_genes`, `sex_summary`, `fingerprint_genes`, `fingerprint_distances`, `pairing_summary`, `male_fraction`, `male_fraction_genes`
- Figures: `pairing_sex_fingerprint`

Caveats:

- Expression-based fingerprints are clues, not genotype evidence; expression of polymorphic genes also depends on cell composition.
- pairing_summary orients the question of who the donors are to whoever holds the sample records; it does not resolve it. With k donors (or pools) the same-group pairs are the nearest by chance with probability up to 1/k! (1/2 with 2, 1/6 with 3).
- A sex mismatch within a hypothesised pair (XIST against chrY genes) is the strongest expression clue against the pairing; matching sex is weak evidence for it.
- The male fraction assumes similar chrY expression across samples and donors and is relative to the male reference samples of the config.

Parameters (`modules.pairing`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `sex_genes` | gene_list_ref | `"sex_genes"` | Sex-chromosome genes of the sex check (those on chrX or chrY). |
| `xist_genes` | gene_list_ref | `"xist"` | XIST gene(s). |
| `group_column` | string or null | `null` | Sample sheet column with the hypothesised donor or pool (equal values = same donor). Its pairs are labelled in fingerprint_distances and checked in pairing_summary. |
| `across_factor` | string or null | `null` | Two-level factor that each donor is expected to span (e.g. the sorted fraction). |
| `fingerprint_variants` | array of object | `[{"id":"fingerprint","genes":"fingerprint_genes"}]` | Gene sets of the expression fingerprint, each giving pairwise sample distances. |
| `fingerprint_variants[].id` | string |  | Variant id. |
| `fingerprint_variants[].genes` | gene_list_ref |  | Genes of the fingerprint. |
| `fingerprint_variants[].normalization` | string (cpm, cpm_tmm) | `"cpm"` | CPM on all genes, or TMM-normalised CPM on the genes kept by gene_filter (sex genes kept). |
| `fingerprint_variants[].min_total_count` | integer | `0` | Genes with fewer reads over all samples are dropped. |
| `male_reference_samples` | array of string | `[]` | sample_ids taken as all-male material (XIST about 0); male fractions are relative to their mean. Empty: the male fraction is not estimated. |
| `sex_denominators` | array of object | `[{"id":"chr_y"}]` | Male-fraction estimates, relative to male_reference_samples. |
| `sex_denominators[].id` | string |  | Estimate id (column male_fraction_<id>). |
| `sex_denominators[].numerator` | gene_list_ref |  | Genes summed in the numerator (default the chrY genes of sex_genes). |
| `sex_denominators[].denominator` | gene_list_ref |  | Genes summed in the denominator (default none). |

## 17 power

Status: implemented (`run_power`). Output folder: `17_power/`.

- Data levels: `counts`
- Optional levels: `gene_mask`, `gene_universe`
- Resources: none
- Conditions: `design_df_min: 1`
- Provides: none
- Tables: `bcv`, `sample_size_grid`
- Figures: none

Caveats:

- BCV from 1 residual df has a per-gene 95% interval spanning more than an order of magnitude (x0.45 to x31.9).
- The bootstrap of the common BCV resamples genes: it measures stability to the choice of genes, not biological variability.
- Sample sizes come from the Hart et al. 2013 approximation and the assumed CVs of the config; they are indicative.

Parameters (`modules.power`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `design` | string or null | `null` | Id of a designs[] entry (default: the first). |
| `universes` | array of object | `[{"id":"all","subset":"all"}]` | Gene sets of the BCV, as subset expressions over the filtered genes (e.g. in_reference_universe+no_masked); estimated in this order on one random stream. |
| `universes[].id` | string |  | Universe id used in outputs. |
| `universes[].subset` | string |  | Subset expression over the filtered genes. |
| `bootstrap_replicates` | integer | `200` | Bootstrap samples of genes for the interval of the common BCV. |
| `reference_scenario` | object | `{}` | Scenario of the n_hart_* columns of the bcv table. |
| `reference_scenario.mu` | number | `100` | Mean count. |
| `reference_scenario.fold_change` | number | `2` | Fold change to detect. |
| `reference_scenario.alpha` | number | `0.01` | Two-sided significance level. |
| `grid` | object | `{}` | Grid of sample_size_grid. |
| `grid.cv` | array of number | `[0.2,0.3,0.4,0.6]` | Biological coefficients of variation. |
| `grid.cv_labels` | array of string | `[]` | One label per cv value (the assumption behind it); empty for none. |
| `grid.mu` | array of number | `[20,100]` | Mean counts. |
| `grid.fold_change` | array of number | `[1.5,2]` | Fold changes to detect. |
| `grid.alpha` | array of number | `[0.01,0.001]` | Two-sided significance levels. |
| `grid.power` | number | `0.8` | Target power. |

## 18 coverage

Status: not implemented (`run_coverage`). Output folder: `18_coverage/`.

- Data levels: `coverage`
- Optional levels: none
- Resources: `annotation`
- Conditions: none
- Provides: none
- Tables: `strand_assignment`, `coverage_features`, `coverage_targets`, `coverage_target_summary`, `coverage_target_exons`, `coverage_totals`
- Figures: `coverage_target`

Caveats:

- bedGraph values from STAR are RPM of mapped reads; ratios between samples depend on that normalisation.

Parameters (`modules.coverage`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |

## 19 covariates

Status: not implemented (`run_covariates`). Output folder: `19_covariates/`.

- Data levels: `contrast_stats`
- Optional levels: none
- Resources: `annotation`, `transcript_fasta`
- Conditions: `min_contrasts: 2`
- Provides: none
- Tables: `gene_covariates`, `covariate_residual_concordance`, `covariate_correlations`
- Figures: none

Caveats:

- Regressing lfc on length, GC and expression also removes biology correlated with them.

Parameters (`modules.covariates`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |

## 20 reference_axis

Status: implemented (`run_reference_axis`). Output folder: `20_reference_axis/`.

- Data levels: `contrast_stats`
- Optional levels: `counts`, `gene_mask`, `gene_universe`
- Resources: none
- Conditions: `configured: axes, gene_sets.sets, sign_test.threshold`
- Provides: none
- Tables: `reference_axis`, `axis_genes`
- Figures: none

Caveats:

- Axes come from external references (other platforms and tissues): correlations with them are descriptive.
- The binomial p-value of the sign test treats genes as independent, which they are not: conditional_on_contrast and descriptive.

Parameters (`modules.reference_axis`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `estimator` | string or null | `null` | Method of the contrast tables whose log2 fold changes are used (null: the first method available). |
| `contrasts` | array of string | `[]` | Contrast ids, in this order (empty = every contrast of the estimator, config order). |
| `subset` | string | `"all"` | Genes (subset expression over the genes of the contrast tables) of the axis correlations, the sign test and the axis-high and axis-low genes of the sample scores. |
| `set_subset` | string | `"all"` | Genes of the mean fold change per gene set. |
| `axes` | array of object | `[]` | Reference axes: log2((mean numerator + pseudocount) / (mean denominator + pseudocount)) of the cell types of an HPA export, per gene. |
| `axes[].id` | string |  | Axis id used in outputs. |
| `axes[].resource` | string |  | Key of an hpa_immune_cell or hpa_single_cell export. |
| `axes[].numerator` | array of string |  | Cell types of the numerator (values of the export's cell column). |
| `axes[].denominator` | array of string |  | Cell types of the denominator. |
| `axes[].pseudocount` | number | `1` | Added to the mean expressions before the log2 ratio. |
| `axes[].threshold` | number | `1` | Genes with axis >= threshold (high) and <= -threshold (low) define the sample scores. |
| `gene_sets` | object | `{}` | Gene sets of an MSigDB export whose mean fold change is reported. |
| `gene_sets.resource` | string or null | `null` | Key of an msigdb_export in resources.items. |
| `gene_sets.sets` | array of string | `[]` | Sets whose mean fold change is reported. |
| `sign_test` | object | `{}` | Sign agreement of contrast pairs. |
| `sign_test.threshold` | number or null | `null` | \|log2 fold change\| that both contrasts of each contrast pair must exceed; the fraction with the same sign gets a two-sided binomial p-value (descriptive: genes are not independent). null: no test. |
| `sample_scores` | boolean | `true` | Per-sample axis scores (needs counts). |
| `prior_count` | number | `2` | Prior count of the log2 CPM of the sample scores. |
| `norm_method` | string (TMM, RLE, upperquartile, none) | `"TMM"` | Library normalisation of the sample scores. |

## 21 de_heatmap

Status: implemented (`run_de_heatmap`). Output folder: `21_de_heatmap/`.

- Data levels: `contrast_stats`, `counts`
- Optional levels: `normalized`
- Resources: none
- Conditions: `significance_available: true`
- Provides: none
- Tables: `de_heatmap_genes`, `de_heatmap_values`
- Figures: `de_heatmap`

Caveats:

- Rows are z-scores across the samples shown: colours compare samples within a gene, not genes with each other.
- Genes are selected on the same samples, so their separation of the groups is expected and only descriptive.

Parameters (`modules.de_heatmap`):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` | Run the module when its data levels are present. |
| `contrasts` | array of string | `[]` | Contrast ids drawn, in this order (empty = every contrast with a significance metric). |
| `n_top` | integer | `50` | Most significant genes classified UP or DOWN by the project thresholds, per contrast. |
| `per_direction` | boolean | `false` | Take up to half of n_top from each direction (UP and DOWN), so that the strongest direction does not fill the heatmap. |
| `values` | string (auto, normalized, median_of_ratios) | `"auto"` | Normalised counts shown. normalized: the normalized input; median_of_ratios: counts divided by DESeq2 size factors computed on the genes of the contrast table; auto: normalized when present, else median_of_ratios. |
| `cluster_rows` | boolean | `true` | Order the genes by hierarchical clustering (complete linkage) of their z-scores. |
| `z_limit` | number | `2.5` | Limit of the colour scale; larger absolute z-scores are drawn at the limit. |
