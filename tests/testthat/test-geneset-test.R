geneset_fixture <- function(n = 300, seed = 1) {
  set.seed(seed)
  ids <- sprintf("ENSG%011d", seq_len(n))
  up <- seq_len(30)
  lfc_a <- rnorm(n)
  lfc_a[up] <- lfc_a[up] + 2
  lfc_b <- rnorm(n)
  lfc_b[up] <- lfc_b[up] + 2
  ct <- rbind(
    new_contrast_table("a", ids, lfc_a,
      gene_symbol = paste0("G", seq_len(n)), method = "m", source = "computed",
      extra = data.frame(ave_log2cpm = seq(10, 1, length.out = n))
    ),
    new_contrast_table("b", ids, lfc_b,
      gene_symbol = paste0("G", seq_len(n)), method = "m", source = "computed",
      extra = data.frame(ave_log2cpm = seq(10, 1, length.out = n))
    )
  )
  gs <- data.frame(
    gs_collection = c(rep("H", 60), rep("C2", 40), rep("C7", 40)),
    gs_subcollection = c(rep("", 60), rep("CP:REACTOME", 40), rep("IMMUNESIGDB", 40)),
    gs_name = c(
      rep(c("HALLMARK_UP", "HALLMARK_RANDOM"), each = 30), rep(c("REACTOME_X", "REACTOME_Y"), each = 20),
      rep(c("GSE1_KEEP", "GSE2_DROP"), each = 20)
    ),
    ensembl_gene = ids[c(1:30, 151:180, 31:70, 71:110)], gene_symbol = "", db_version = "v1"
  )
  f <- tempfile(fileext = ".tsv")
  data.table::fwrite(gs, f, sep = "\t")
  ctx <- list(
    cfg = list(
      contrasts = list(list(id = "b"), list(id = "a")),
      contrast_pairs = list(list(id = "ab", first = "a", second = "b"))
    ),
    organism = "human", seed = 7, resources = list(gene_sets = f)
  )
  list(ct = ct, ctx = ctx)
}

geneset_params <- list(
  methods = c("fgsea", "camera_pr"), min_size = 10, inter_gene_cor = c(0.01, 0.05), use_ranks = c(FALSE, TRUE),
  collections = list(
    list(id = "H", collection = "H"),
    list(id = "REACTOME", collection = "C2", subcollection = "CP:REACTOME"),
    list(id = "C7_GSE1", collection = "C7", subcollection = "IMMUNESIGDB", set_pattern = "^GSE1_")
  )
)

test_that("geneset_test runs fgsea and cameraPR variants with their inference scopes", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("limma")
  x <- geneset_fixture()
  res <- run_geneset_test(list(contrast_stats = x$ct), geneset_params, x$ctx)
  expect_registry_outputs(res, "geneset_test")
  r <- res$tables$geneset_results
  expect_setequal(unique(r$set_id), c("HALLMARK_UP", "HALLMARK_RANDOM", "REACTOME_X", "REACTOME_Y", "GSE1_KEEP"))
  expect_equal(unique(r$contrast_id), c("b", "a"))
  expect_equal(unique(r$inference_scope[r$method == "fgsea"]), "ranking_only")
  expect_equal(unique(r$inference_scope[r$method == "camera_pr"]), "conditional_on_contrast")
  expect_equal(nrow(r), 5 * 2 * (1 + 4))
  expect_equal(r$collection[r$set_id == "REACTOME_X"][[1]], "REACTOME")
  up <- r[r$set_id == "HALLMARK_UP", ]
  expect_true(all(up$direction == "Up"))
  expect_true(all(up$padj < 0.05))
  # padj is BH over all collections of one block
  blk <- r[r$method == "camera_pr" & r$inter_gene_cor == 0.05 & r$use_ranks & r$contrast_id == "a", ]
  expect_equal(blk$padj, stats::p.adjust(blk$pvalue, "BH"))
  fg <- r[r$method == "fgsea" & r$set_id == "HALLMARK_UP", ]
  expect_true(all(lengths(strsplit(fg$leading_edge, ",")) <= 15))
  expect_true(all(startsWith(fg$leading_edge, "G")))
  cmpn <- res$tables$geneset_method_comparison
  expect_equal(nrow(cmpn), 2 * 5)
  expect_true(all(cmpn$n_sets == 5))
  conc <- res$tables$geneset_concordance
  expect_equal(nrow(conc), 5 * 3)
  h <- conc[conc$method == "fgsea" & conc$collection == "H", ]
  expect_equal(h$n_concordant, 1)
  expect_equal(h$concordant_sets, "HALLMARK_UP")
  sets <- res$tables$geneset_concordant_sets
  expect_true("HALLMARK_UP" %in% sets$set_id)
  expect_named(res$figures, c("nes_concordance", "nes_heatmap"))
  expect_s3_class(res$figures$nes_concordance, "ggplot")
  expect_equal(res$meta$geneset_results$seed, 7)
  expect_equal(res$meta$geneset_results$collections[[3]]$n_sets, 1)
  write_module_outputs(res, tempfile(), list())
})

test_that("geneset_test is reproducible with the seed and honours its parameters", {
  skip_if_not_installed("fgsea")
  x <- geneset_fixture()
  p <- utils::modifyList(geneset_params, list(methods = "fgsea", contrasts = "a"))
  r1 <- run_geneset_test(list(contrast_stats = x$ct), p, x$ctx)$tables$geneset_results
  r2 <- run_geneset_test(list(contrast_stats = x$ct), p, x$ctx)$tables$geneset_results
  expect_identical(r1, r2)
  expect_equal(unique(r1$contrast_id), "a")
  r3 <- run_geneset_test(list(contrast_stats = x$ct), utils::modifyList(p, list(seed = 99)), x$ctx)
  expect_equal(r3$meta$geneset_results$seed, 99)
  expect_null(r3$tables$geneset_concordance)
})

test_that("geneset_test selects universes from gene flags", {
  skip_if_not_installed("limma")
  x <- geneset_fixture()
  flags <- data.frame(feature_id = unique(x$ct$feature_id), masked = FALSE)
  flags$masked[1:20] <- TRUE
  p <- utils::modifyList(geneset_params, list(
    methods = "camera_pr", inter_gene_cor = 0.01, use_ranks = FALSE,
    universes = list(list(id = "all", subset = "all"), list(id = "unmasked", subset = "no_masked"))
  ))
  res <- run_geneset_test(list(contrast_stats = x$ct, gene_flags = flags), p, x$ctx)
  r <- res$tables$geneset_results
  expect_equal(unique(r$n_ranked[r$universe == "unmasked"]), 280)
  expect_equal(r$size[r$universe == "unmasked" & r$set_id == "HALLMARK_UP"][[1]], 10)
  expect_equal(res$meta$geneset_results$universes[[2]]$n_genes, 280)
  expect_null(res$figures$nes_concordance)
})

test_that("geneset_test fails loudly on bad configuration", {
  x <- geneset_fixture()
  run <- function(p, ctx = x$ctx) run_geneset_test(list(contrast_stats = x$ct), p, ctx)
  expect_error(run(list(estimator = "nope")), "estimator 'nope'")
  expect_error(run(list(contrasts = "zz")), "without estimates: zz")
  expect_error(run(list(collections = list(list(id = "X", collection = "C5")))), "'X' selects no gene set")
  expect_error(
    run(list(collections = list(list(id = "A", collection = "H"), list(id = "B", collection = "H")))),
    "more than one collection"
  )
  expect_error(run(list(resource = "other")), "resource 'other'")
  ctx <- x$ctx
  ctx$cfg$annotation$id_type <- "entrez_id"
  expect_error(run(list(), ctx), "Ensembl gene IDs or gene symbols")
  expect_error(run(list(universes = list(list(id = "u", subset = "all"), list(id = "u", subset = "all")))), "duplicated")
})

test_that("module tables accept only known inference scopes", {
  res <- list(tables = list(t = data.frame(pvalue = 0.1, inference_scope = "whatever")), figures = list(), meta = list())
  expect_error(write_module_outputs(res, tempfile()), "invalid inference_scope whatever")
  res$tables$t$inference_scope <- "ranking_only"
  expect_silent(write_module_outputs(res, tempfile()))
})

test_that("camera_pr_ranks() is identical to limma::cameraPR(use.ranks = TRUE)", {
  set.seed(3)
  with_ties <- stats::setNames(round(stats::rnorm(500), 1), paste0("g", 1:500))
  no_ties <- stats::setNames(stats::rnorm(500), paste0("g", 1:500))
  index <- list(one = 7L, small = c(1L, 5L, 9L), large = sample(500, 120), shifted = order(no_ties)[1:40])
  for (stat in list(with_ties, no_ties)) {
    for (igc in c(0, 0.01, 0.05)) {
      expect_identical(
        camera_pr_ranks(stat, index, igc),
        limma::cameraPR(stat, index, inter.gene.cor = igc, use.ranks = TRUE, sort = FALSE)
      )
    }
  }
})
