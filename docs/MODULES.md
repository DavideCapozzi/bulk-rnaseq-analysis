# Modules

<!-- Generated from inst/modules.yaml and inst/schema/schema.yaml by `Rscript tools/make_modules_doc.R`; do not edit by hand (a test fails when this file is out of date). -->

Every module declares the data levels it needs. A module whose needs are not met, or that is not implemented yet, is skipped with the reason in `00_provenance/capabilities.json` and `report/summary.md`. Defaults come from the config schema; `bulkrnaseq plan <config>` shows what a project runs.

18 modules: 12 implemented, 6 not implemented.

| NN | Module | Status | Data levels | Provides |
|---|---|---|---|---|
| 02 | [`qc`](#02-qc) | implemented | one of `counts`, `qc` | none |
| 03 | [`pca`](#03-pca) | implemented | `counts` | none |
| 04 | [`markers`](#04-markers) | not implemented | `counts` | none |
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
| 15 | [`signatures`](#15-signatures) | not implemented | `counts` | none |
| 16 | [`pairing`](#16-pairing) | implemented | `counts` | none |
| 17 | [`power`](#17-power) | implemented | `counts` | none |
| 18 | [`coverage`](#18-coverage) | not implemented | `coverage` | none |
| 19 | [`covariates`](#19-covariates) | not implemented | `contrast_stats` | none |

Schema types used below:

- `file_or_pattern`: Path relative to inputs.base_dir. May contain the placeholder {input}, replaced by each value of the sample sheet column input_column (one file per sample).
- `gene_list_ref`: Key of a gene list in resources.gene_lists, or an inline list of gene symbols.
- `seed`: integer
- `threshold_list`: array of number

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
| `enabled` | boolean | `true` |  |
| `top_n_share` | integer | `50` |  |
| `top_n_list` | integer | `15` |  |
| `ribosomal_protein_regex` | string | `"^RP[SL][0-9]+[A-Z]?[0-9]*$"` |  |
| `hemoglobin_genes` | gene_list_ref | `"hemoglobin"` |  |
| `rrna_biotypes` | array of string | `["rRNA","Mt_rRNA","rRNA_pseudogene"]` |  |
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
| `enabled` | boolean | `true` |  |
| `transforms` | array of string (vst, logcpm_tmm) | `["vst","logcpm_tmm"]` |  |
| `ntop` | array of integer or string | `[500,1000,2000,"all"]` |  |
| `sex_gene_variants` | array of string | `["keep","drop_xist_y","drop_xy"]` |  |
| `gene_subsets` | array of string | `"all"` |  |
| `symbol_exclusion_variants` | array of string | `[]` | Regexes on gene symbols (e.g. "^HLA-"); each adds a variant of the grid without the matching genes, besides the variant without exclusion. |
| `thinning` | boolean | `false` | Adds a variant on counts thinned to the smallest library (edgeR::thinCounts, project.seed). |
| `separation_factors` | array of string | `[]` | Factors labelled in pca_summary (default: all samples.factors, in order). |
| `orient_pc1` | object or null | `null` | Orients PC1 so that samples at this level have positive mean score (drivers, pc1_order). |
| `abundance_sex_genes` | string (keep, drop_xist_y, drop_xy) | `"keep"` |  |
| `scores_ntop` | integer or string | `1000` |  |
| `pc_drivers` | boolean | `true` |  |
| `qc_metrics_for_drivers` | array of string | `[]` |  |
| `top_loadings` | integer | `40` |  |

## 04 markers

Status: not implemented (`run_markers`). Output folder: `04_markers/`.

- Data levels: `counts`
- Optional levels: `abundance`
- Resources: `annotation`, `marker_panels`
- Conditions: none
- Provides: none
- Tables: `marker_expression`, `marker_panel_scores`, `marker_ratios`
- Figures: `markers_heatmap`, `marker_panels`

Caveats:

- Panels are starting markers, not validated signatures. CPM are library-size only (no TMM) unless stated.

Parameters (`modules.markers`; draft, not read by any code yet):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` |  |
| `panels` | string | `"marker_panels"` |  |
| `pseudocount` | number | `1` |  |
| `ratios` | array of object | `[]` |  |
| `ratios[].id` | string |  |  |
| `ratios[].numerator_gene` | string |  |  |
| `ratios[].denominator_gene` | string |  |  |
| `ratios[].numerator_sample` | string |  |  |
| `ratios[].denominator_sample` | string |  |  |

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
| `enabled` | boolean | `true` |  |
| `profiles` | array of object | `[]` | Gene profiles of contaminating (or reference) cells, scored per sample on CPM of all genes. |
| `profiles[].id` | string |  |  |
| `profiles[].genes` | gene_list_ref |  |  |
| `reference_profile` | string or null | `null` | Profile id whose CPM sum is also reported relative to the sample with the largest value (rel_to_max). |
| `ratio_pseudocount` | number | `1` | CPM added before the log2 ratios of ratio_checks. |
| `ratio_checks` | array of object | `[]` | Per-gene log2 CPM ratios between two samples over a profile's genes. |
| `ratio_checks[].id` | string |  |  |
| `ratio_checks[].profile` | string |  |  |
| `ratio_checks[].numerator_sample` | string |  |  |
| `ratio_checks[].denominator_sample` | string |  |  |
| `excess_checks` | array of object | `[]` | CPM of genes above the mean of baseline_samples, divided by the CPM sum of a profile, in two samples, and the ratio of the two excesses (a similar ratio for genes and profile is compatible with the genes coming from the profile's cells). |
| `excess_checks[].id` | string |  |  |
| `excess_checks[].genes` | gene_list_ref |  |  |
| `excess_checks[].profile` | string |  |  |
| `excess_checks[].baseline_samples` | array of string |  |  |
| `excess_checks[].samples` | array of string |  |  |
| `mask` | object | `{}` | masked = (in any include group and in no exclude group) or among the genes of add_profiles. Groups are gene sets of an MSigDB export (resource); their genes match by Ensembl ID (by symbol when annotation.id_type is gene_symbol). |
| `mask.resource` | string or null | `null` | Key of an msigdb_export in resources.items. |
| `mask.include` | array of object | `[]` |  |
| `mask.include[].id` | string |  |  |
| `mask.include[].sets` | array of string |  |  |
| `mask.exclude` | array of object | `[]` |  |
| `mask.exclude[].id` | string |  |  |
| `mask.exclude[].sets` | array of string |  |  |
| `mask.add_profiles` | array of string | `[]` |  |

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

Parameters (`modules.contamination_sets`; draft, not read by any code yet):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` |  |
| `method` | string (camera_pr) | `"camera_pr"` |  |
| `inter_gene_cor` | number | `0.01` |  |

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
| `enabled` | boolean | `true` |  |
| `profiles` | array of object | `[]` | Reference expression of one cell type from an HPA export: nTPM for kind hpa_immune_cell, nCPM for hpa_single_cell (kind of the resource in resources.items). |
| `profiles[].id` | string |  |  |
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

Parameters (`modules.deconvolution`; draft, not read by any code yet):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` |  |
| `components` | object | `{}` |  |
| `marker_min_expr` | number | `50` |  |
| `marker_specificity` | number | `10` |  |
| `solver` | string (optim_lbfgsb, lawson_hanson) | `"optim_lbfgsb"` |  |
| `variants` | array | `[]` |  |

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
| `enabled` | boolean | `true` |  |
| `estimators` | array of string (logcpm_difference, deseq2_apeglm) | `"logcpm_difference"` |  |
| `prior_count` | number | `2` |  |
| `norm_method` | string (TMM, RLE, upperquartile, none, reference) | `"TMM"` |  |
| `reference_genes` | null or | `null` |  |
| `dispersion_design` | string or null | `null` |  |
| `spread_subsets` | array of string | `"all"` | Gene subsets of contrast_spread. |
| `spread_threshold` | number | `1` |  |

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
| `enabled` | boolean | `true` |  |
| `prefilter` | object | `{"min_count":10,"min_samples":3}` | Genes with at least min_count reads in at least min_samples samples are tested. |
| `prefilter.min_count` | integer |  |  |
| `prefilter.min_samples` | integer |  |  |
| `results_alpha` | number | `0.1` | alpha of DESeq2::results(), the FDR target of independent filtering (DESeq2 default 0.1). |
| `compare_with` | array of object | `[]` | Computed design vs imported contrast table, reported in provider_comparison. |
| `compare_with[].design` | string |  |  |
| `compare_with[].contrast_table` | string |  |  |

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
| `enabled` | boolean | `true` |  |
| `thresholds` | threshold_list | `[0.58,1]` |  |
| `subsets` | array of string | `"all"` |  |

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
| `enabled` | boolean | `true` |  |
| `resource` | string | `"gene_sets"` | Key of an msigdb_export in resources.items. |
| `estimator` | string or null | `null` | Method of the contrast tables whose log2 fold changes rank the genes (e.g. logcpm_difference); null: the first method available. |
| `contrasts` | array of string | `[]` | Contrast ids to test, in this order (empty = every contrast of the estimator, config order). |
| `methods` | array of string (camera_pr, fgsea) | `"camera_pr"` |  |
| `collections` | array of object | `[{"id":"H","collection":"H"}]` | Gene set collections of the export: gs_collection, optional gs_subcollection and a regular expression on gs_name. P-values are adjusted over all collections together. |
| `collections[].id` | string |  |  |
| `collections[].collection` | string |  |  |
| `collections[].subcollection` | string or null | `null` |  |
| `collections[].set_pattern` | string or null | `null` |  |
| `universes` | array of object | `[{"id":"all","subset":"all"}]` | Ranked genes, as subset expressions over the genes of the contrast tables; tested in this order. |
| `universes[].id` | string |  |  |
| `universes[].subset` | string |  |  |
| `min_size` | integer | `15` |  |
| `max_size` | integer | `500` |  |
| `inter_gene_cor` | array of number | `0.01` | cameraPR variants. |
| `use_ranks` | array of boolean | `false` | cameraPR variants. |
| `fgsea_tie_jitter_sd` | number | `1e-06` | SD of the normal noise added to the ranking statistic before fgsea, to break ties reproducibly (0: none). Genes are ordered by decreasing expression, then feature_id, before the draw. |
| `seed` | integer or null | `null` | Seed set once before fgsea (universes outer, contrasts inner, on one random stream); null: project.seed. |
| `leading_edge_max` | integer | `15` | Leading-edge symbols reported per set. |
| `alpha` | number | `0.05` |  |
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
| `enabled` | boolean | `true` |  |
| `databases` | array of string (GO_BP, GO_CC, GO_MF, KEGG, Reactome) | `["GO_BP","GO_CC","GO_MF","KEGG","Reactome"]` |  |
| `p_adjust_method` | string | `"BH"` |  |
| `pvalue_cutoff` | number | `0.05` |  |
| `show_category` | integer | `15` |  |

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
| `enabled` | boolean | `true` |  |
| `label_top` | integer | `15` |  |

## 15 signatures

Status: not implemented (`run_signatures`). Output folder: `15_signatures/`.

- Data levels: `counts`
- Optional levels: none
- Resources: `gene_sets`
- Conditions: none
- Provides: none
- Tables: `signature_scores`
- Figures: `signature_scores`

Caveats:

- z-score is relative to the samples in the run; singscore and ssGSEA are within-sample ranks.

Parameters (`modules.signatures`; draft, not read by any code yet):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` |  |
| `methods` | array of string (zscore, ssgsea, singscore) | `["zscore","singscore"]` |  |
| `sets` | array of string | `[]` |  |
| `composites` | array | `[]` |  |
| `min_size` | integer | `10` |  |
| `universes` | array of string | `"all"` |  |

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
| `enabled` | boolean | `true` |  |
| `sex_genes` | gene_list_ref | `"sex_genes"` |  |
| `xist_genes` | gene_list_ref | `"xist"` |  |
| `group_column` | string or null | `null` | Sample sheet column with the hypothesised donor or pool (equal values = same donor). Its pairs are labelled in fingerprint_distances and checked in pairing_summary. |
| `across_factor` | string or null | `null` | Two-level factor that each donor is expected to span (e.g. the sorted fraction). |
| `fingerprint_variants` | array of object | `[{"id":"fingerprint","genes":"fingerprint_genes"}]` |  |
| `fingerprint_variants[].id` | string |  |  |
| `fingerprint_variants[].genes` | gene_list_ref |  |  |
| `fingerprint_variants[].normalization` | string (cpm, cpm_tmm) | `"cpm"` | CPM on all genes, or TMM-normalised CPM on the genes kept by gene_filter (sex genes kept). |
| `fingerprint_variants[].min_total_count` | integer | `0` | Genes with fewer reads over all samples are dropped. |
| `male_reference_samples` | array of string | `[]` | sample_ids taken as all-male material (XIST about 0); male fractions are relative to their mean. Empty: the male fraction is not estimated. |
| `sex_denominators` | array of object | `[{"id":"chr_y"}]` |  |
| `sex_denominators[].id` | string |  |  |
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
| `enabled` | boolean | `true` |  |
| `design` | string or null | `null` | Id of a designs[] entry (default: the first). |
| `universes` | array of object | `[{"id":"all","subset":"all"}]` | Gene sets of the BCV, as subset expressions over the filtered genes (e.g. in_reference_universe+no_masked); estimated in this order on one random stream. |
| `universes[].id` | string |  |  |
| `universes[].subset` | string |  |  |
| `bootstrap_replicates` | integer | `200` |  |
| `reference_scenario` | object | `{}` | Scenario of the n_hart_* columns of the bcv table. |
| `reference_scenario.mu` | number | `100` |  |
| `reference_scenario.fold_change` | number | `2` |  |
| `reference_scenario.alpha` | number | `0.01` |  |
| `grid` | object | `{}` |  |
| `grid.cv` | array of number | `[0.2,0.3,0.4,0.6]` |  |
| `grid.cv_labels` | array of string | `[]` | One label per cv value (the assumption behind it); empty for none. |
| `grid.mu` | array of number | `[20,100]` |  |
| `grid.fold_change` | array of number | `[1.5,2]` |  |
| `grid.alpha` | array of number | `[0.01,0.001]` |  |
| `grid.power` | number | `0.8` |  |

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

Parameters (`modules.coverage`; draft, not read by any code yet):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` |  |
| `target_genes` | gene_list_ref | `[]` |  |
| `normalizer_genes` | gene_list_ref | `[]` |  |
| `flank_bp` | integer | `10000` |  |
| `window_bp` | integer | `500` |  |

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

Parameters (`modules.covariates`; draft, not read by any code yet):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enabled` | boolean | `true` |  |
| `spline_df` | integer | `4` |  |
| `strata` | array of string | `["ge_median","ge_q75","all"]` |  |
