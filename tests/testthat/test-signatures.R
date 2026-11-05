signatures_fixture <- function() {
  se <- nb_se(n_genes = 400, groups = c("a", "a", "b", "b"), seed = 3)
  ids <- rownames(se)
  SummarizedExperiment::colData(se)$group <- c("a", "a", "b", "b")
  f <- tempfile(fileext = ".tsv")
  data.table::fwrite(data.frame(
    gs_collection = "H", gs_subcollection = "", gs_name = rep(c("SET_UP", "SET_OTHER", "SET_TINY"), c(30, 30, 3)),
    ensembl_gene = c(ids[1:30], ids[201:230], ids[301:303]), gene_symbol = c(ids[1:30], ids[201:230], ids[301:303]),
    db_version = "v1"
  ), f, sep = "\t")
  cfg <- list(
    samples = list(factors = list(group = c("a", "b"))),
    gene_filter = list(min_cpm = 0, min_samples = 1, biotypes = NULL, exclude_mitochondrial = FALSE, sex_genes = "keep"),
    contrasts = list(
      list(id = "b_vs_a", numerator = list("group=b"), denominator = list("group=a")),
      list(id = "s3_vs_s1", numerator = list("s3"), denominator = list("s1"))
    ),
    contrast_pairs = list(list(id = "p", first = "b_vs_a", second = "s3_vs_s1"))
  )
  ctx <- list(cfg = cfg, organism = "human", resources = list(gene_sets = f), gene_lists = list())
  list(se = se, ctx = ctx, ids = ids)
}

signatures_params <- list(
  resource = "gene_sets", msigdb_sets = list("SET_UP", "SET_OTHER", "SET_TINY"),
  composites = list(list(id = "up_minus_other", up = "SET_UP", down = "SET_OTHER")),
  methods = list("zscore", "ssgsea", "singscore"), min_size = 5
)

test_that("signatures scores sets, composites, contrast differences and pair agreement", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("singscore")
  x <- signatures_fixture()
  res <- run_signatures(
    list(se = x$se, features = as.data.frame(SummarizedExperiment::rowData(x$se))),
    signatures_params, x$ctx
  )
  expect_registry_outputs(res, "signatures")
  s <- res$tables$signature_scores
  expect_equal(unique(s$method), c("zscore", "ssgsea", "singscore"))
  expect_equal(unique(s$set_id), c("SET_UP", "SET_OTHER", "up_minus_other"))
  expect_equal(res$meta$signature_scores$sets_below_min_size, "all:SET_TINY")
  expect_equal(nrow(s), 3 * 3 * 4)
  # The first 30 genes are up in group b: the set score separates the groups.
  z <- s[s$method == "zscore" & s$set_id == "SET_UP", ]
  expect_true(all(z$score[3:4] > z$score[1:2]))
  zo <- s[s$method == "zscore" & s$set_id == "SET_OTHER", ]
  zc <- s[s$method == "zscore" & s$set_id == "up_minus_other", ]
  expect_equal(zc$score, z$score - zo$score)
  expect_true(all(is.na(zc$n_genes)))
  expect_equal(z$n_genes, rep(30L, 4))
  d <- res$tables$signature_differences
  dz <- d[d$method == "zscore" & d$set_id == "SET_UP", ]
  expect_equal(dz$contrast_id, c("b_vs_a", "s3_vs_s1"))
  expect_equal(dz$difference, c(mean(z$score[3:4]) - mean(z$score[1:2]), z$score[3] - z$score[1]))
  expect_equal(dz$n_numerator, c(2L, 1L))
  pc <- res$tables$signature_method_comparison
  expect_equal(nrow(pc), 3 * 3)
  expect_equal(pc$concordant, sign(pc$first_difference) == sign(pc$second_difference))
  expect_s3_class(res$figures$signature_scores, "ggplot")
})

test_that("signatures matches GSVA and singscore called directly", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("singscore")
  x <- signatures_fixture()
  res <- run_signatures(list(se = x$se), signatures_params, x$ctx)
  lc <- log_cpm(SummarizedExperiment::assay(x$se, "counts"))
  sets <- list(SET_UP = x$ids[1:30], SET_OTHER = x$ids[201:230])
  gs <- GSVA::gsva(GSVA::ssgseaParam(lc, sets, minSize = 5), verbose = FALSE)
  s <- res$tables$signature_scores
  got <- s$score[s$method == "ssgsea" & s$set_id == "SET_UP"]
  expect_equal(got, unname(gs["SET_UP", ]))
  rk <- singscore::rankGenes(lc)
  both <- singscore::simpleScore(rk, upSet = sets$SET_UP, downSet = sets$SET_OTHER)$TotalScore
  expect_equal(s$score[s$method == "singscore" & s$set_id == "up_minus_other"], both)
})

test_that("signatures takes custom sets and universes from gene flags", {
  x <- signatures_fixture()
  flags <- data.frame(feature_id = x$ids, masked = seq_along(x$ids) %in% 1:10)
  p <- list(
    custom_sets = list(list(id = "custom", genes = x$ids[1:30])),
    universes = list(list(id = "all", subset = "all"), list(id = "unmasked", subset = "no_masked"))
  )
  res <- run_signatures(list(se = x$se, gene_flags = flags), p, x$ctx)
  s <- res$tables$signature_scores
  expect_equal(unique(s$n_genes[s$universe == "unmasked"]), 20L)
  expect_equal(unique(s$n_genes[s$universe == "all"]), 30L)
})

test_that("signatures fails loudly on bad configuration and runs only when configured", {
  x <- signatures_fixture()
  run <- function(p) run_signatures(list(se = x$se), p, x$ctx)
  expect_error(run(list()), "no sets configured")
  expect_error(run(list(msigdb_sets = list("SET_UP"))), "needs resource")
  expect_error(run(list(resource = "zz", msigdb_sets = list("SET_UP"))), "resource 'zz'")
  expect_error(run(list(resource = "gene_sets", msigdb_sets = list("NOPE"))), "not found: NOPE")
  expect_error(
    run(list(custom_sets = list(list(id = "a", genes = "g0001")), composites = list(list(id = "c", up = "a", down = "b")))),
    "unknown set\\(s\\): b"
  )
  expect_error(run(list(custom_sets = list(list(id = "a", genes = "g0001")), figure_method = "ssgsea")), "figure_method")
  cfg <- make_cfg(list(
    project = list(name = "x", output_dir = "o"), inputs = list(counts = list(adapter = "matrix", file = "c.tsv")),
    resources = list(items = list())
  ))
  caps <- resolve_capabilities(cfg)
  expect_match(caps$reason[caps$module == "signatures"], "needs modules.signatures.msigdb_sets or")
})
