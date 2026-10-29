# Four samples of 1e6 reads, so CPM = counts and the expected values are computed by hand.
contam_se <- function() {
  genes <- c("P1", "M1", "T1", "X1", "M2", "L1", "FILLER")
  m <- cbind(
    s1 = c(10, 100, 1000, 50, 5, 20, 0),
    s2 = c(20, 300, 500, 90, 5, 20, 0),
    s3 = c(30, 100, 1000, 50, 5, 20, 0),
    s4 = c(40, 300, 250, 70, 5, 20, 0)
  )
  m[length(genes), ] <- 1e6 - colSums(m)
  rownames(m) <- sprintf("ENSG%011d", seq_along(genes))
  storage.mode(m) <- "integer"
  rd <- S4Vectors::DataFrame(
    feature_id = rownames(m), gene_symbol = genes, biotype = "protein_coding", chromosome = "1",
    row.names = rownames(m)
  )
  cd <- S4Vectors::DataFrame(sample_id = colnames(m), g = c("a", "b", "a", "b"), row.names = colnames(m))
  SummarizedExperiment::SummarizedExperiment(list(counts = m), rowData = rd, colData = cd)
}
id_of <- function(se, sym) rownames(se)[match(sym, SummarizedExperiment::rowData(se)$gene_symbol)]

contam_context <- function(se) {
  sets <- tempfile(fileext = ".tsv")
  ids <- function(s) id_of(se, s)
  writeLines(c(
    "gs_collection\tgs_subcollection\tgs_name\tensembl_gene\tgene_symbol\tdb_version",
    sprintf("C8\t\tSET_MY\t%s\t%s\tv1", ids(c("M1", "M2")), c("M1", "M2")),
    sprintf("C8\t\tSET_T\t%s\t%s\tv1", ids(c("M2", "T1")), c("M2", "T1"))
  ), sets)
  imm <- tempfile(fileext = ".tsv")
  # T-reg / classical monocyte per gene; X1 is absent from the reference.
  v <- data.frame(
    sym = c("P1", "M1", "T1", "M2", "L1", "FILLER"), treg = c(0.5, 50, 100, 5, 0.2, 2),
    mono = c(0, 400, 10, 0, 0, 0)
  )
  writeLines(c(
    "Gene\tGene name\tImmune cell\tTPM\tpTPM\tnTPM",
    sprintf("%s\t%s\tT-reg\t0\t0\t%s", ids(v$sym), v$sym, v$treg),
    sprintf("%s\t%s\tclassical monocyte\t0\t0\t%s", ids(v$sym), v$sym, v$mono)
  ), imm)
  sc <- tempfile(fileext = ".tsv")
  writeLines(c("Gene\tGene name\tCell type\tnCPM", sprintf("%s\tM2\tmacrophages\t40", ids("M2"))), sc)
  items <- list(
    sets = list(source = "file://x", kind = "msigdb_export", version = "test sets"),
    imm = list(source = "file://x", kind = "hpa_immune_cell"), sc = list(source = "file://x", kind = "hpa_single_cell")
  )
  cfg <- list(resources = list(items = items), gene_filter = list(min_cpm = 2, min_samples = 2, sex_genes = "keep"))
  list(
    cfg = cfg, organism = "human", gene_lists = read_gene_lists(),
    resources = list(sets = sets, imm = imm, sc = sc)
  )
}
contam_params <- function() {
  list(
    profiles = list(list(id = "prof_a", genes = c("P1", "M1")), list(id = "prof_t", genes = "T1")),
    reference_profile = "prof_t",
    ratio_checks = list(list(id = "s3_over_s1", profile = "prof_a", numerator_sample = "s3", denominator_sample = "s1")),
    excess_checks = list(list(
      id = "x", genes = "X1", profile = "prof_a", baseline_samples = c("s1", "s3"), samples = c("s2", "s4")
    )),
    mask = list(
      resource = "sets", include = list(list(id = "myeloid", sets = "SET_MY")),
      exclude = list(list(id = "t_nk", sets = "SET_T")), add_profiles = "prof_a"
    )
  )
}

test_that("contamination scores profiles on CPM, checks ratios and excesses, and builds the mask", {
  se <- contam_se()
  res <- run_contamination(list(se = se), contam_params(), contam_context(se))
  expect_registry_outputs(res, "contamination")
  sc <- res$tables$contamination_scores
  a <- sc[sc$profile_id == "prof_a", ]
  expect_equal(a$sum_cpm, c(110, 320, 130, 340))
  expect_equal(a$mean_log2_cpm_p1[[1]], mean(log2(c(11, 101))))
  expect_true(all(is.na(a$rel_to_max)))
  expect_equal(sc$rel_to_max[sc$profile_id == "prof_t"], c(1, 0.5, 1, 0.25))
  rc <- res$tables$profile_ratio_consistency
  expect_equal(rc$log2_ratio, c(log2(31 / 11), 0))
  ex <- res$tables$profile_excess
  g <- ex[ex$row_type == "gene", ]
  expect_equal(c(g$excess_per_profile_s2, g$excess_per_profile_s4), c(40 / 320, 20 / 340))
  expect_equal(g$ratio, 2)
  expect_equal(ex$ratio[ex$row_type == "profile"], (320 - 120) / (340 - 120))
  expect_equal(ex$cpm_s2[ex$row_type == "profile"], 320)
  mk <- res$tables$contaminant_mask
  expect_setequal(mk$gene_symbol, c("P1", "M1", "T1", "M2"))
  # M2 is in the myeloid set but also in the T/NK set: kept. P1 is masked as a profile gene.
  expect_equal(mk$masked[match(c("P1", "M1", "T1", "M2"), mk$gene_symbol)], c(TRUE, TRUE, FALSE, FALSE))
  fl <- res$gene_flags
  expect_equal(nrow(fl), nrow(se))
  expect_equal(fl$masked_by_sets[match(id_of(se, c("M1", "P1")), fl$feature_id)], c(TRUE, FALSE))
  expect_equal(fl$masked_by_profiles[match(id_of(se, c("M1", "P1")), fl$feature_id)], c(TRUE, TRUE))
  expect_equal(res$meta$contaminant_mask$n_kept_by_exclude, 1)
  expect_equal(res$meta$contaminant_mask$db_version, "v1")
  expect_s3_class(res$figures$contamination_profiles, "ggplot")
})

test_that("contamination refuses unknown profiles, samples and gene sets", {
  se <- contam_se()
  ctx <- contam_context(se)
  p <- contam_params()
  p$reference_profile <- "nope"
  expect_error(run_contamination(list(se = se), p, ctx), "reference_profile 'nope' is not among the profile ids")
  p <- contam_params()
  p$ratio_checks[[1]]$numerator_sample <- "s9"
  expect_error(run_contamination(list(se = se), p, ctx), "not among the samples with counts: s9")
  p <- contam_params()
  p$mask$include[[1]]$sets <- "SET_ZZ"
  expect_error(run_contamination(list(se = se), p, ctx), "set\\(s\\) not found: SET_ZZ")
  expect_error(run_contamination(list(se = se), list(), ctx), "no profiles and no mask gene sets configured")
  p <- contam_params()
  p$profiles[[2]]$id <- "prof_a"
  expect_error(run_contamination(list(se = se), p, ctx), "duplicated profiles id")
})

test_that("reference builds the universe from the references only and explains candidate exclusions", {
  se <- contam_se()
  ctx <- contam_context(se)
  mask <- run_contamination(list(se = se), contam_params(), ctx)$gene_flags
  params <- list(
    profiles = list(
      list(id = "t_reg", resource = "imm", cell = "T-reg"), list(id = "mono", resource = "imm", cell = "classical monocyte"),
      list(id = "macro", resource = "sc", cell = "macrophages")
    ),
    target = "t_reg", contaminants = c("mono", "macro"),
    exclusion_candidates = c("T1", "M1", "P1", "X1", "M2", "L1", "NOPE")
  )
  res <- run_reference(list(se = se, gene_flags = combine_gene_flags(list(contamination = mask))), params, ctx)
  expect_registry_outputs(res, "reference")
  rp <- res$tables$reference_profiles
  at <- function(col, sym) rp[[col]][match(sym, rp$gene_symbol)]
  expect_equal(at("contam_max", c("M1", "M2", "X1")), c(400, 40, NA))
  expect_equal(at("contam_sensitive", c("M1", "M2", "T1")), c(TRUE, TRUE, FALSE))
  # The universe ignores the mask: M1 is masked but its exclusion here is the contaminant reference.
  expect_equal(rp$gene_symbol[rp$reference_universe], c("T1", "FILLER"))
  ex <- res$tables$universe_exclusions
  why <- stats::setNames(ex$exclusion_reason, ex$gene_symbol)
  expect_equal(
    why[c("T1", "M1", "P1", "X1", "M2", "L1", "NOPE")],
    c(
      T1 = "included", M1 = "mask_gene_sets", P1 = "mask_profiles", X1 = "absent_in_reference",
      M2 = "contam_sensitive", L1 = "target_below_min", NOPE = "not_filtered"
    )
  )
  expect_true(res$meta$universe_exclusions$mask_applied)
  u <- res$tables$reference_universe
  expect_equal(u$n_reference_universe_unmasked[u$genes == "all_features"], 2)
  expect_equal(names(res$gene_flags), c("feature_id", "target_expressed", "contam_sensitive", "reference_universe"))
  # Without the mask the reasons use the references only.
  res2 <- run_reference(list(se = se), params, ctx)
  ex2 <- res2$tables$universe_exclusions
  expect_equal(ex2$exclusion_reason[ex2$gene_symbol == "M1"], "contam_sensitive")
  expect_false(res2$meta$universe_exclusions$mask_applied)
  expect_true(is.na(res2$tables$reference_universe$n_reference_universe_unmasked[[1]]))
  params$target <- "zz"
  expect_error(run_reference(list(se = se), params, ctx), "target 'zz' is not among the profile ids")
})

test_that("reference profile reader keeps the requested cells and refuses missing ones", {
  f <- fixture("resources", "hpa_immune_cell.tsv.zip")
  x <- read_reference_profiles(f, "hpa_immune_cell", c("T-reg", "classical monocyte"))
  expect_setequal(unique(x$cell), c("T-reg", "classical monocyte"))
  expect_type(x$value, "double")
  expect_error(read_reference_profiles(f, "hpa_immune_cell", "T-regs"), "cell type\\(s\\) not found: T-regs")
})

test_that("gene flags of several modules are combined, absent genes FALSE, duplicated flags refused", {
  a <- data.frame(feature_id = c("g1", "g2"), masked = c(TRUE, FALSE))
  b <- data.frame(feature_id = c("g2", "g3"), reference_universe = c(TRUE, NA))
  x <- combine_gene_flags(list(a = a, b = b))
  expect_equal(x$feature_id, c("g1", "g2", "g3"))
  expect_equal(x$masked, c(TRUE, FALSE, FALSE))
  expect_equal(x$reference_universe, c(FALSE, TRUE, FALSE))
  expect_error(combine_gene_flags(list(a = a, c = a)), "provided by more than one module: masked")
})

test_that("gene subsets leave genes without expression out of the strata", {
  g <- data.frame(expression = c(1, 2, 3, NA), masked = c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(gene_subset(g, "ge_median"), c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(gene_subset(g, "no_masked+ge_median"), c(FALSE, TRUE, FALSE, FALSE))
})

test_that("power estimates the BCV on each universe of the config, in order", {
  se <- nb_se(groups = c("a", "b", "a", "b"))
  SummarizedExperiment::colData(se)$h <- c("x", "x", "y", "y")
  flags <- data.frame(feature_id = rownames(se), keep_me = seq_len(nrow(se)) %% 2 == 0)
  cfg <- list(designs = list(list(id = "add", formula = "~ g + h")), gene_filter = list(min_cpm = 1, min_samples = 2))
  ctx <- list(cfg = cfg, organism = "human", gene_lists = list(), seed = 3)
  params <- list(
    bootstrap_replicates = 5,
    universes = list(list(id = "all", subset = "all"), list(id = "half", subset = "in_keep_me"))
  )
  res <- run_power(list(se = se, gene_flags = flags), params, ctx)
  b <- res$tables$bcv
  expect_equal(b$universe, c("all", "half"))
  expect_true(b$n_genes[[2]] < b$n_genes[[1]])
  # The first universe draws the same bootstrap as when it is alone.
  alone <- run_power(list(se = se), list(bootstrap_replicates = 5), ctx)$tables$bcv
  expect_identical(b[1, names(alone)], alone)
  expect_equal(res$meta$bcv$universes[[2]]$subset, "in_keep_me")
  params$universes[[2]]$id <- "all"
  expect_error(run_power(list(se = se, gene_flags = flags), params, ctx), "duplicated universe id")
})

test_that("contamination and reference run only when configured, and need the resources their parameters name", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("sample_id\tinput_column", "a\tA", "b\tB", "c\tC"), f)
  base <- list(
    project = list(name = "x", output_dir = "o"), samples = list(sheet = f),
    inputs = list(counts = list(adapter = "matrix", file = "x.tsv")), resources = list(items = list(gtf = list(source = "file://g")))
  )
  caps <- resolve_capabilities(make_cfg(base))
  expect_match(
    caps$reason[caps$module == "contamination"],
    "needs modules.contamination.profiles or modules.contamination.mask.include in the config"
  )
  expect_match(caps$reason[caps$module == "reference"], "needs modules.reference.profiles in the config")
  cfg <- base
  cfg$modules <- list(
    contamination = list(profiles = list(list(id = "p", genes = "pan_t"))),
    reference = list(profiles = list(list(id = "t", resource = "hpa", cell = "T-reg")), target = "t")
  )
  caps <- resolve_capabilities(make_cfg(cfg))
  expect_equal(caps$status[caps$module == "contamination"], "run")
  expect_match(caps$reason[caps$module == "reference"], "needs resource\\(s\\) 'hpa' named in modules.reference")
  cfg$resources$items$hpa <- list(source = "file://h", kind = "hpa_immune_cell")
  caps <- resolve_capabilities(make_cfg(cfg), read_sample_sheet(f, list()))
  expect_equal(caps$depends_on[caps$module == "reference"], "ingest,contamination")
  expect_equal(caps$depends_on[caps$module == "pca"], "ingest,contamination,reference")
})
