#' Gene set test module
#'
#' Tests gene set collections of an MSigDB export on the log2 fold changes of each contrast (one estimator), in
#' one or more gene universes:
#' * `fgsea` (pre-ranked GSEA, multilevel, `eps = 0`, one worker): genes are ordered by decreasing expression, then
#'   `feature_id`, normal noise of SD `fgsea_tie_jitter_sd` is added to break ties, and the seed is set once
#'   before the loop (universes outer, contrasts inner). Gene permutations make its p-values a ranking of sets
#'   within one contrast (`inference_scope = "ranking_only"`).
#' * `camera_pr` (`limma::cameraPR`) for each `inter_gene_cor` and `use_ranks`: it allows for inter-gene
#'   correlation but still tests genes for a given contrast (`"conditional_on_contrast"`).
#'
#' P-values are BH-adjusted over all sets of all collections, per method variant, universe and contrast. For each
#' `contrast_pairs` entry whose two contrasts are tested, sets significant in both contrasts with the same
#' direction are listed (`geneset_concordant_sets`) and counted per collection (`geneset_concordance`).
#'
#' @param inputs `list` with `contrast_stats` (canonical table), optional `features` and `gene_flags`.
#' @param params `list`. `modules.geneset_test` from the config.
#' @param context `list` with `cfg`, `organism`, `seed` and `resources` (path of `params$resource`).
#' @return `list(tables, figures, meta, figure_sizes)`; tables `geneset_results`, `geneset_method_comparison`,
#'   and with contrast pairs `geneset_concordance`, `geneset_concordant_sets`; figures `nes_concordance` and
#'   `nes_heatmap` when fgsea runs.
#' @examples
#' \donttest{
#' set.seed(1)
#' ids <- sprintf("ENSG%011d", 1:200)
#' ct <- new_contrast_table("a", ids, rnorm(200), method = "m", source = "computed")
#' f <- tempfile(fileext = ".tsv")
#' data.table::fwrite(data.frame(
#'   gs_collection = "H", gs_subcollection = "", gs_name = rep(c("SET_1", "SET_2"), each = 20),
#'   ensembl_gene = ids[1:40], gene_symbol = "", db_version = "v1"
#' ), f, sep = "\t")
#' res <- run_geneset_test(
#'   list(contrast_stats = ct), list(min_size = 5),
#'   list(cfg = list(), organism = "human", seed = 1, resources = list(gene_sets = f))
#' )
#' res$tables$geneset_results
#' }
#' @export
run_geneset_test <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "geneset_test"))
  ct <- inputs$contrast_stats
  if (is.null(ct) || nrow(ct) == 0) abort("geneset_test: no contrast tables")
  est <- params$estimator %||% unique(ct$method)[[1]]
  if (!est %in% ct$method) {
    abort("geneset_test: estimator '", est, "' not among the contrast tables' methods: ", toString(unique(ct$method)))
  }
  ct <- ct[ct$method == est, , drop = FALSE]
  contrasts <- geneset_contrast_order(unlist(params$contrasts), ct, context$cfg)
  id_type <- context$cfg$annotation$id_type %||% "ensembl_gene_id"
  if (!id_type %in% c("ensembl_gene_id", "gene_symbol")) {
    abort("geneset_test: needs Ensembl gene IDs or gene symbols as feature IDs (annotation.id_type is ", id_type, ")")
  }
  path <- context$resources[[params$resource]]
  if (is.null(path)) abort("geneset_test: resource '", params$resource, "' is not available in the cache")
  gs <- read_gene_sets(path)
  sets <- geneset_members(gs, params$collections, if (id_type == "gene_symbol") "gene_symbol" else "ensembl_gene")
  # split() as in the original scripts: one list element per set name, in sorted order.
  gsl <- split(sets$member, sets$set_id)
  coll_of <- stats::setNames(sets$collection, sets$set_id)
  genes <- contrast_gene_table(ct, inputs$features, inputs$gene_flags)
  genes <- genes[order(-genes$expression, genes$feature_id, method = "radix"), , drop = FALSE]
  symbols <- ct$gene_symbol[match(genes$feature_id, ct$feature_id)]
  uids <- vapply(params$universes, `[[`, "", "id")
  if (anyDuplicated(uids)) abort("geneset_test: duplicated universe id(s): ", toString(unique(uids[duplicated(uids)])))
  universes <- lapply(params$universes, function(u) gene_subset(genes, u$subset, context$organism %||% "human"))
  names(universes) <- uids
  # Ranking statistic of one contrast in one universe, in the gene order above (genes without lfc dropped).
  ranked <- function(u, k) {
    x <- ct[ct$contrast_id == k, c("feature_id", "lfc")]
    keep <- universes[[u]]
    lfc <- x$lfc[match(genes$feature_id[keep], x$feature_id)]
    ok <- is.finite(lfc)
    list(stat = stats::setNames(lfc[ok], genes$feature_id[keep][ok]), symbol = symbols[keep][ok])
  }
  methods <- unlist(params$methods)
  seed <- params$seed %||% context$seed %||% 1
  res <- list()
  if ("fgsea" %in% methods) {
    if (!requireNamespace("fgsea", quietly = TRUE) || !requireNamespace("BiocParallel", quietly = TRUE)) {
      abort("geneset_test: method 'fgsea' needs the fgsea and BiocParallel packages")
    }
    res$fgsea <- withr::with_seed(seed, do.call(rbind, lapply(uids, function(u) {
      do.call(rbind, lapply(contrasts, function(k) {
        r <- ranked(u, k)
        geneset_fgsea(gsl, r$stat, r$symbol, params, u, k)
      }))
    })))
  }
  if ("camera_pr" %in% methods) {
    if (!requireNamespace("limma", quietly = TRUE)) abort("geneset_test: method 'camera_pr' needs the limma package")
    res$camera_pr <- do.call(rbind, lapply(uids, function(u) {
      do.call(rbind, lapply(contrasts, function(k) {
        r <- ranked(u, k)
        idx <- lapply(gsl, function(g) which(names(r$stat) %in% g))
        idx <- idx[lengths(idx) >= params$min_size & lengths(idx) <= params$max_size]
        do.call(rbind, lapply(unlist(params$inter_gene_cor), function(igc) {
          do.call(rbind, lapply(unlist(params$use_ranks), function(rk) {
            geneset_camera(r$stat, idx, igc, rk, u, k)
          }))
        }))
      }))
    }))
  }
  results <- do.call(rbind, c(list(geneset_empty_results()), res[intersect(methods, names(res))]))
  results$collection <- unname(coll_of[results$set_id])
  results <- results[, names(geneset_empty_results())]
  rownames(results) <- NULL
  alpha <- params$alpha
  variant <- geneset_variant(results)
  comparison <- unique(results[, c("method", "inter_gene_cor", "use_ranks", "universe", "contrast_id")])
  key <- paste(variant, results$universe, results$contrast_id, sep = "\r")
  ckey <- paste(geneset_variant(comparison), comparison$universe, comparison$contrast_id, sep = "\r")
  comparison$n_sets <- as.integer(table(factor(key, levels = ckey)))
  comparison$n_padj_below_alpha <- vapply(ckey, function(z) sum(results$padj[key == z] < alpha, na.rm = TRUE), 0L)
  rownames(comparison) <- NULL
  tables <- list(geneset_results = results, geneset_method_comparison = comparison)
  pairs <- Filter(function(p) all(c(p$first, p$second) %in% contrasts), context$cfg$contrast_pairs %||% list())
  if (length(pairs) > 0) {
    conc <- geneset_pair_concordance(results, pairs, alpha, unlist(params$concordance_collections))
    if (!is.null(conc$summary)) tables$geneset_concordance <- conc$summary
    tables$geneset_concordant_sets <- conc$sets
  }
  figs <- list()
  sizes <- list()
  fg <- results[results$method == "fgsea", , drop = FALSE]
  if (nrow(fg) > 0 && length(pairs) > 0) {
    fu <- params$figure_universe %||% uids[[1]]
    if (!fu %in% uids) abort("geneset_test: figure_universe '", fu, "' is not a universe id")
    fc <- unlist(params$figure_collections)
    if (length(fc) == 0) fc <- unique(sets$collection)
    figs$nes_concordance <- plot_nes_concordance(fg[fg$universe == fu & fg$collection %in% fc, ], pairs, alpha, fu)
    sizes$nes_concordance <- c(6.5 * length(pairs), 7)
  }
  hc <- params$heatmap_collection %||% params$collections[[1]]$id
  hm <- fg[fg$collection == hc, , drop = FALSE]
  if (nrow(hm) > 0 && any(hm$padj < alpha, na.rm = TRUE)) {
    figs$nes_heatmap <- plot_nes_heatmap(hm, alpha, contrasts, uids)
    n_sets <- length(unique(hm$set_id[hm$set_id %in% hm$set_id[hm$padj < alpha]]))
    sizes$nes_heatmap <- c(4 + 1.5 * length(contrasts) * length(uids), max(4, 2 + 0.2 * n_sets))
  }
  item <- context$cfg$resources$items[[params$resource]]
  list(
    tables = tables, figures = figs, figure_sizes = sizes,
    meta = list(
      geneset_results = list(
        description = paste0(
          "Gene set tests on the log2 fold changes of '", est, "'; padj = BH over all sets of all collections per ",
          "method variant, universe and contrast. direction: sign of NES (fgsea) or cameraPR Direction."
        ),
        estimator = est, resource = params$resource, resource_version = item$version %||% NA_character_,
        db_version = sort(unique(gs$db_version)),
        match_by = if (id_type == "gene_symbol") "gene_symbol" else "ensembl_gene",
        collections = lapply(params$collections, function(cl) {
          list(
            id = cl$id, collection = cl$collection, subcollection = cl$subcollection %||% NA_character_,
            set_pattern = cl$set_pattern %||% NA_character_,
            n_sets = sum(sets$collection[!duplicated(sets$set_id)] == cl$id)
          )
        }),
        universes = lapply(params$universes, function(u) {
          list(id = u$id, subset = u$subset, n_genes = sum(universes[[u$id]]))
        }),
        seed = if ("fgsea" %in% methods) seed else NA, fgsea_tie_jitter_sd = params$fgsea_tie_jitter_sd,
        packages = geneset_package_versions(methods),
        caveats = c(
          "fgsea permutes genes: p-values rank sets within one contrast, not across samples (ranking_only).",
          "cameraPR tests genes for a given contrast (conditional_on_contrast), not samples."
        )
      ),
      geneset_method_comparison = list(
        description = paste0("Sets tested and sets with padj < ", alpha, " per method variant, universe and contrast."),
        alpha = alpha
      ),
      geneset_concordance = list(
        description = paste0(
          "Per contrast pair and collection: sets tested in both contrasts, correlation of NES (fgsea), sets with ",
          "padj < ", alpha, " in both and, among them, with the same direction (descriptive)."
        ),
        alpha = alpha
      ),
      geneset_concordant_sets = list(
        description = paste0("Sets with padj < ", alpha, " in both contrasts of a pair and the same direction."),
        alpha = alpha
      )
    )
  )
}

# Contrast ids to test: the configured list (all must exist), else every contrast of the estimator in config order
# followed by contrasts absent from the config (imported) in table order.
geneset_contrast_order <- function(wanted, ct, cfg) {
  avail <- unique(ct$contrast_id)
  if (length(wanted) > 0) {
    missing <- setdiff(wanted, avail)
    if (length(missing) > 0) abort("geneset_test: contrast(s) without estimates: ", toString(missing))
    return(wanted)
  }
  cfg_ids <- vapply(cfg$contrasts %||% list(), `[[`, "", "id")
  c(intersect(cfg_ids, avail), setdiff(avail, cfg_ids))
}

# Members of the configured collections: data.frame(set_id, member, collection). A set may belong to one
# collection only.
geneset_members <- function(gs, collections, key) {
  out <- lapply(collections, function(cl) {
    keep <- gs$gs_collection %in% cl$collection
    if (!is.null(cl$subcollection)) keep <- keep & gs$gs_subcollection %in% cl$subcollection
    if (!is.null(cl$set_pattern)) keep <- keep & grepl(cl$set_pattern, gs$gs_name)
    if (!any(keep)) abort("geneset_test: collection '", cl$id, "' selects no gene set of the export")
    x <- gs[keep, , drop = FALSE]
    data.frame(set_id = x$gs_name, member = x[[key]], collection = cl$id, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  owners <- unique(out[, c("set_id", "collection")])
  dup <- unique(owners$set_id[duplicated(owners$set_id)])
  if (length(dup) > 0) {
    abort("geneset_test: set(s) selected by more than one collection: ", toString(utils::head(dup, 5)))
  }
  out[!is.na(out$member) & out$member != "", , drop = FALSE]
}

geneset_empty_results <- function() {
  data.frame(
    method = character(), inter_gene_cor = numeric(), use_ranks = logical(), universe = character(),
    contrast_id = character(), collection = character(), set_id = character(), size = integer(), es = numeric(),
    nes = numeric(), direction = character(), pvalue = numeric(), padj = numeric(), log2err = numeric(),
    n_ranked = integer(), leading_edge = character(), inference_scope = character(), stringsAsFactors = FALSE
  )
}

geneset_fgsea <- function(gsl, stat, symbol, params, universe, contrast) {
  keep <- !duplicated(names(stat))
  stat <- stat[keep]
  symbol <- symbol[keep]
  if (params$fgsea_tie_jitter_sd > 0) stat <- stat + stats::rnorm(length(stat), sd = params$fgsea_tie_jitter_sd)
  # The backend fgsea builds itself for nproc = 1 on Unix (one worker), without its progress bar.
  bp <- BiocParallel::MulticoreParam(workers = 1, progressbar = FALSE)
  fg <- fgsea::fgsea(gsl, stat, minSize = params$min_size, maxSize = params$max_size, eps = 0, BPPARAM = bp)
  if (nrow(fg) == 0) {
    return(NULL)
  }
  le <- vapply(fg$leadingEdge, function(g) {
    paste(utils::head(symbol[match(g, names(stat))], params$leading_edge_max), collapse = ",")
  }, "")
  out <- data.frame(
    method = "fgsea", inter_gene_cor = NA_real_, use_ranks = NA, universe = universe, contrast_id = contrast,
    collection = NA_character_, set_id = fg$pathway, size = as.integer(fg$size), es = fg$ES, nes = fg$NES,
    direction = ifelse(fg$NES > 0, "Up", ifelse(fg$NES < 0, "Down", NA_character_)), pvalue = fg$pval,
    padj = fg$padj, log2err = fg$log2err, n_ranked = length(stat), leading_edge = le,
    inference_scope = "ranking_only", stringsAsFactors = FALSE
  )
  out[order(out$set_id, method = "radix"), , drop = FALSE]
}

geneset_camera <- function(stat, idx, igc, use_ranks, universe, contrast) {
  if (length(idx) == 0) {
    return(NULL)
  }
  cp <- limma::cameraPR(stat, idx, inter.gene.cor = igc, use.ranks = use_ranks, sort = FALSE)
  out <- data.frame(
    method = "camera_pr", inter_gene_cor = igc, use_ranks = use_ranks, universe = universe, contrast_id = contrast,
    collection = NA_character_, set_id = rownames(cp), size = as.integer(cp$NGenes), es = NA_real_, nes = NA_real_,
    direction = as.character(cp$Direction), pvalue = cp$PValue, padj = stats::p.adjust(cp$PValue, "BH"),
    log2err = NA_real_, n_ranked = length(stat), leading_edge = NA_character_,
    inference_scope = "conditional_on_contrast", stringsAsFactors = FALSE
  )
  out[order(out$set_id, method = "radix"), , drop = FALSE]
}

# Label of a method variant (fgsea, or cameraPR with its inter-gene correlation and use of ranks).
geneset_variant <- function(x) {
  ifelse(x$method == "camera_pr", paste0("camera_pr|", x$inter_gene_cor, "|", x$use_ranks), x$method)
}

geneset_pair_concordance <- function(results, pairs, alpha, collections) {
  variant <- geneset_variant(results)
  summary <- list()
  sets <- list()
  for (v in unique(variant)) {
    r <- results[variant == v, , drop = FALSE]
    for (u in unique(r$universe)) {
      ru <- r[r$universe == u, , drop = FALSE]
      for (p in pairs) {
        a <- ru[ru$contrast_id == p$first, , drop = FALSE]
        b <- ru[ru$contrast_id == p$second, , drop = FALSE]
        ids <- sort(union(a$set_id, b$set_id), method = "radix")
        ia <- match(ids, a$set_id)
        ib <- match(ids, b$set_id)
        coll <- ifelse(is.na(ia), b$collection[ib], a$collection[ia])
        both <- a$padj[ia] < alpha & b$padj[ib] < alpha
        both[is.na(both)] <- FALSE
        same <- both & a$direction[ia] == b$direction[ib]
        same[is.na(same)] <- FALSE
        head <- data.frame(
          method = r$method[[1]], inter_gene_cor = r$inter_gene_cor[[1]], use_ranks = r$use_ranks[[1]],
          universe = u, pair_id = p$id, first = p$first, second = p$second, stringsAsFactors = FALSE
        )
        if (any(same)) {
          w <- which(same)
          sets[[length(sets) + 1]] <- cbind(head, data.frame(
            collection = coll[w], set_id = ids[w], direction = a$direction[ia][w], nes_first = a$nes[ia][w],
            nes_second = b$nes[ib][w], padj_first = a$padj[ia][w], padj_second = b$padj[ib][w],
            inference_scope = a$inference_scope[ia][w], stringsAsFactors = FALSE
          ))
        }
        colls <- if (length(collections) > 0) collections else unique(results$collection)
        for (cl in colls) {
          i <- coll == cl
          tested <- i & !is.na(ia) & !is.na(ib)
          r_nes <- if (all(is.na(a$nes[ia][tested]))) {
            NA_real_
          } else {
            suppressWarnings(stats::cor(a$nes[ia][tested], b$nes[ib][tested], use = "complete.obs"))
          }
          summary[[length(summary) + 1]] <- cbind(head, data.frame(
            collection = cl, n_sets = sum(tested), r_nes = r_nes, n_significant_both = sum(both[i]),
            n_concordant = sum(same[i]), concordant_sets = paste(ids[i & same], collapse = ","),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  empty_sets <- data.frame(
    method = character(), inter_gene_cor = numeric(), use_ranks = logical(), universe = character(),
    pair_id = character(), first = character(), second = character(), collection = character(),
    set_id = character(), direction = character(), nes_first = numeric(), nes_second = numeric(),
    padj_first = numeric(), padj_second = numeric(), inference_scope = character(), stringsAsFactors = FALSE
  )
  list(
    summary = if (length(summary)) do.call(rbind, summary) else NULL,
    sets = do.call(rbind, c(list(empty_sets), sets))
  )
}

geneset_package_versions <- function(methods) {
  pk <- c(fgsea = "fgsea", camera_pr = "limma")[intersect(c("fgsea", "camera_pr"), methods)]
  as.list(stats::setNames(vapply(pk, function(p) as.character(utils::packageVersion(p)), ""), pk))
}

#' NES of two contrasts per contrast pair (fgsea)
#'
#' @param fg fgsea rows of `geneset_results` (one universe).
#' @param pairs `list` of contrast pairs (`id`, `first`, `second`).
#' @param alpha `numeric(1)`. Significance level of both contrasts.
#' @param universe `character(1)`. Universe id (subtitle).
#' @return A `ggplot` object faceted by pair.
#' @examples
#' fg <- data.frame(
#'   contrast_id = rep(c("a", "b"), each = 3), set_id = rep(c("S1", "S2", "S3"), 2),
#'   collection = "H", nes = c(2, -1, 0.5, 1.8, 1, 0.2), padj = c(0.01, 0.5, 0.9, 0.02, 0.4, 0.8)
#' )
#' plot_nes_concordance(fg, list(list(id = "ab", first = "a", second = "b")), 0.05, "all")
#' @export
plot_nes_concordance <- function(fg, pairs, alpha, universe) {
  d <- do.call(rbind, lapply(pairs, function(p) {
    a <- fg[fg$contrast_id == p$first, ]
    b <- fg[fg$contrast_id == p$second, ]
    ids <- intersect(a$set_id, b$set_id)
    ia <- match(ids, a$set_id)
    ib <- match(ids, b$set_id)
    sig <- a$padj[ia] < alpha & b$padj[ib] < alpha & sign(a$nes[ia]) == sign(b$nes[ib])
    data.frame(
      pair = p$id, set_id = ids, collection = a$collection[ia], x = a$nes[ia], y = b$nes[ib],
      class = ifelse(!is.na(sig) & sig, "concordant", "other"), stringsAsFactors = FALSE
    )
  }))
  d$class <- factor(d$class, levels = c("other", "concordant"))
  d <- d[order(d$class), ]
  lab <- d[d$class == "concordant", ]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#c9c8c3", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#c9c8c3", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$class, shape = .data$collection), size = 2.4) +
    ggrepel::geom_text_repel(
      data = lab, ggplot2::aes(label = .data$set_id), size = 2.3, colour = "#1b1b1b", seed = 3,
      max.overlaps = 25, min.segment.length = 0
    ) +
    ggplot2::scale_colour_manual(
      values = c(other = "#b5b4af", concordant = "#2a78d6"), drop = FALSE,
      name = paste0("padj < ", alpha, " in both, same sign")
    ) +
    ggplot2::facet_wrap(~pair) +
    ggplot2::labs(
      x = "NES (first contrast)", y = "NES (second contrast)", shape = NULL, subtitle = paste("Universe:", universe)
    ) +
    theme_bulkrnaseq() +
    ggplot2::theme(legend.position = "bottom")
}

#' Heatmap of NES per contrast (fgsea, one collection)
#'
#' Sets with padj below `alpha` in at least one contrast and universe; a dot marks padj below `alpha`.
#'
#' @param fg fgsea rows of `geneset_results` (one collection).
#' @param alpha `numeric(1)`. Significance level.
#' @param contrasts,universes `character`. Order of the columns and facets.
#' @return A `ggplot` object.
#' @examples
#' fg <- data.frame(
#'   universe = "all", contrast_id = rep(c("a", "b"), each = 2), set_id = rep(c("S1", "S2"), 2),
#'   nes = c(2, -1, 1.5, -0.5), padj = c(0.01, 0.5, 0.03, 0.6)
#' )
#' plot_nes_heatmap(fg, 0.05, c("a", "b"), "all")
#' @export
plot_nes_heatmap <- function(fg, alpha, contrasts, universes) {
  keep <- unique(fg$set_id[!is.na(fg$padj) & fg$padj < alpha])
  d <- fg[fg$set_id %in% keep, , drop = FALSE]
  ord <- stats::aggregate(d$nes, list(set_id = d$set_id), mean)
  d$set_id <- factor(d$set_id, levels = ord$set_id[order(ord$x)])
  d$contrast_id <- factor(d$contrast_id, levels = contrasts)
  d$universe <- factor(d$universe, levels = universes)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$contrast_id, y = .data$set_id, fill = .data$nes)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_point(data = d[!is.na(d$padj) & d$padj < alpha, ], shape = 16, size = 1.1, colour = "#1b1b1b") +
    ggplot2::facet_wrap(~universe) +
    ggplot2::scale_fill_gradient2(low = "#2a78d6", mid = "#f4f3ef", high = "#d1495b", midpoint = 0, name = "NES") +
    ggplot2::labs(x = NULL, y = NULL, subtitle = paste0("dot: padj < ", alpha)) +
    theme_bulkrnaseq() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
