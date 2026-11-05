#' Signatures module
#'
#' Per-sample scores of gene signatures (MSigDB sets and custom gene lists), with their differences between the
#' two sides of each configured contrast and the sign agreement of those differences within each contrast pair.
#' All results are descriptive: no p-values.
#'
#' Scores are computed on log2 CPM (`norm_method` normalisation factors, `prior_count`) of the genes kept by the
#' project's gene filter, restricted to each universe (a subset expression over the gene flags):
#' * `zscore`: per gene, z-score across the samples of the run; the score is the mean z-score of the set's genes
#'   (relative to the samples in the run).
#' * `ssgsea`: `GSVA::gsva(ssgseaParam(..., minSize = min_size))` on all sets of the universe together, so the
#'   default normalisation by the range of the scores depends on the sets scored.
#' * `singscore`: `singscore::simpleScore()` on the within-sample ranks of the universe (`TotalScore`).
#'
#' A composite `{id, up, down}` is `score(up) - score(down)` for `zscore` and `ssgsea`, and the bidirectional
#' `simpleScore(upSet = up, downSet = down)` for `singscore`. Sets with fewer than `min_size` genes in a universe
#' are not scored there (listed in the metadata).
#'
#' @param inputs `list` with `se` (counts), optional `features` and `gene_flags`.
#' @param params `list`. `modules.signatures` from the config.
#' @param context `list` with `cfg`, `organism`, `gene_lists` and `resources` (path of `params$resource`).
#' @return `list(tables, figures, meta, figure_sizes)`; tables `signature_scores` (one row per universe, method,
#'   set and sample), `signature_differences` (per configured contrast: mean score of the numerator samples minus
#'   that of the denominator samples) and `signature_method_comparison` (per contrast pair: both differences and
#'   whether their signs agree); figure `signature_scores`.
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- matrix(rnbinom(800, mu = 200, size = 10), 200, dimnames = list(paste0("G", 1:200), paste0("s", 1:4)))
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   rowData = S4Vectors::DataFrame(
#'     feature_id = rownames(m), gene_symbol = rownames(m), biotype = "protein_coding", chromosome = "1"
#'   ),
#'   colData = S4Vectors::DataFrame(sample_id = colnames(m), row.names = colnames(m))
#' )
#' p <- list(custom_sets = list(list(id = "first", genes = paste0("G", 1:20))))
#' res <- run_signatures(list(se = se), p, list(cfg = list(annotation = list(id_type = "gene_symbol"))))
#' head(res$tables$signature_scores)
#' }
#' @export
run_signatures <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "signatures"))
  se <- inputs$se
  if (is.null(se) || !"counts" %in% SummarizedExperiment::assayNames(se)) abort("signatures: needs counts")
  cfg <- context$cfg %||% list()
  organism <- context$organism %||% "human"
  methods <- unlist(params$methods)
  signatures_check_packages(methods)
  sets <- signatures_sets(params, cfg, context, inputs$features)
  composites <- params$composites %||% list()
  comp_ids <- vapply(composites, `[[`, "", "id")
  all_ids <- c(names(sets), comp_ids)
  if (anyDuplicated(all_ids)) {
    abort("signatures: duplicated set id(s): ", toString(unique(all_ids[duplicated(all_ids)])))
  }
  for (cp in composites) {
    missing <- setdiff(c(cp$up, cp$down), names(sets))
    if (length(missing) > 0) abort("signatures: composite '", cp$id, "' refers to unknown set(s): ", toString(missing))
  }
  keep <- config_filter_genes(se, cfg, context)
  counts <- SummarizedExperiment::assay(se, "counts")[keep, , drop = FALSE]
  lc <- log_cpm(counts, norm_method = params$norm_method, prior_count = params$prior_count)
  genes <- data.frame(feature_id = rownames(lc), expression = rowMeans(lc), stringsAsFactors = FALSE)
  features <- inputs$features
  genes$chromosome <- if (!is.null(features$chromosome)) {
    features$chromosome[match(genes$feature_id, features$feature_id)]
  } else {
    NA_character_
  }
  genes <- merge_gene_flags(genes, inputs$gene_flags)
  uids <- vapply(params$universes, `[[`, "", "id")
  if (anyDuplicated(uids)) abort("signatures: duplicated universe id(s): ", toString(unique(uids[duplicated(uids)])))
  scores <- list()
  dropped <- character()
  universe_sizes <- list()
  for (u in params$universes) {
    x <- lc[gene_subset(genes, u$subset, organism), , drop = FALSE]
    universe_sizes[[u$id]] <- nrow(x)
    members <- lapply(sets, intersect, rownames(x))
    small <- lengths(members) < params$min_size
    dropped <- c(dropped, paste0(u$id, ":", names(members)[small]))
    members <- members[!small]
    for (meth in methods) {
      s <- signatures_score(x, members, composites, meth, params$min_size)
      if (nrow(s) > 0) scores[[length(scores) + 1]] <- data.frame(universe = u$id, method = meth, s)
    }
  }
  long <- signatures_long(scores, colnames(lc))
  samples <- as.data.frame(SummarizedExperiment::colData(se))
  if (is.null(samples$sample_id)) samples$sample_id <- colnames(se)
  diffs <- signatures_differences(long, cfg$contrasts %||% list(), samples)
  pairs <- Filter(function(p) all(c(p$first, p$second) %in% diffs$contrast_id), cfg$contrast_pairs %||% list())
  tables <- list(signature_scores = long, signature_differences = diffs)
  if (length(pairs) > 0) tables$signature_method_comparison <- signatures_pairs(diffs, pairs)
  fu <- params$figure_universe %||% uids[[1]]
  fm <- params$figure_method %||% methods[[1]]
  if (!fu %in% uids) abort("signatures: figure_universe '", fu, "' is not a universe id")
  if (!fm %in% methods) abort("signatures: figure_method '", fm, "' is not among the methods")
  fig <- long[long$universe == fu & long$method == fm, , drop = FALSE]
  f <- sample_factors(se, cfg)
  cd <- SummarizedExperiment::colData(se)
  figures <- list()
  if (nrow(fig) > 0) {
    figures$signature_scores <- plot_signature_scores(fig, colnames(se), sample_labels(se, cfg),
      colour = if (length(f) >= 1) as.character(cd[[f[[1]]]]), shape = if (length(f) >= 2) as.character(cd[[f[[2]]]]),
      subtitle = paste0("method ", fm, ", universe ", fu), titles = f
    )
  }
  item <- if (!is.null(params$resource)) cfg$resources$items[[params$resource]]
  list(
    tables = tables, figures = figures,
    figure_sizes = list(signature_scores = c(9, max(4, 1.5 + 0.3 * length(unique(fig$set_id))))),
    meta = list(
      signature_scores = list(
        description = paste0(
          "Per-sample signature scores on log2 CPM (", params$norm_method, ", prior count ", params$prior_count,
          ") of the genes kept by gene_filter, per universe and method. zscore: mean per-gene z-score across the ",
          "samples of the run; ssgsea: GSVA ssGSEA (normalised over the sets scored together); singscore: ",
          "simpleScore TotalScore. Composites: score(up) - score(down) for zscore and ssgsea, bidirectional ",
          "simpleScore(up, down) for singscore. n_genes: genes of the set in the universe (NA for composites)."
        ),
        resource = params$resource %||% NA_character_, resource_version = item$version %||% NA_character_,
        match_by = signatures_match_by(cfg),
        sets = lapply(names(sets), function(s) list(id = s, n_members = length(sets[[s]]))),
        composites = lapply(composites, function(cp) list(id = cp$id, up = cp$up, down = cp$down)),
        universes = lapply(params$universes, function(u) {
          list(id = u$id, subset = u$subset, n_genes = universe_sizes[[u$id]])
        }),
        min_size = params$min_size, sets_below_min_size = dropped,
        packages = signatures_package_versions(methods),
        caveats = c(
          "zscore is relative to the samples in the run; ssgsea and singscore are within-sample rank scores.",
          "Scores and their differences are descriptive (no p-values)."
        )
      ),
      signature_differences = list(
        description = "Mean score of the numerator samples minus that of the denominator samples, per contrast."
      ),
      signature_method_comparison = list(
        description = "Per contrast pair: the two score differences and whether their signs agree (descriptive)."
      )
    )
  )
}

signatures_check_packages <- function(methods) {
  need <- c(ssgsea = "GSVA", singscore = "singscore")
  for (m in intersect(methods, names(need))) {
    if (!requireNamespace(need[[m]], quietly = TRUE)) {
      abort("signatures: method '", m, "' needs the ", need[[m]], " package")
    }
  }
}

signatures_package_versions <- function(methods) {
  pk <- c(ssgsea = "GSVA", singscore = "singscore")[intersect(methods, c("ssgsea", "singscore"))]
  as.list(stats::setNames(vapply(pk, function(p) as.character(utils::packageVersion(p)), ""), pk))
}

signatures_match_by <- function(cfg) {
  if (identical(cfg$annotation$id_type, "gene_symbol")) "gene_symbol" else "ensembl_gene"
}

# Named list of member feature IDs per set: MSigDB sets (by name, from the resource) then custom sets.
signatures_sets <- function(params, cfg, context, features) {
  msig <- unlist(params$msigdb_sets)
  out <- list()
  if (length(msig) > 0) {
    key <- params$resource
    if (is.null(key)) abort("signatures: msigdb_sets needs resource (an msigdb_export in resources.items)")
    path <- context$resources[[key]]
    if (is.null(path)) abort("signatures: resource '", key, "' is not available in the cache")
    id_type <- cfg$annotation$id_type %||% "ensembl_gene_id"
    if (!id_type %in% c("ensembl_gene_id", "gene_symbol")) {
      abort("signatures: MSigDB sets need Ensembl gene IDs or gene symbols as feature IDs (id_type is ", id_type, ")")
    }
    gs <- read_gene_sets(path, msig)
    member <- gs[[signatures_match_by(cfg)]]
    for (s in msig) out[[s]] <- unique(member[gs$gs_name == s & !is.na(member) & member != ""])
  }
  for (cs in params$custom_sets %||% list()) {
    symbols <- resolve_gene_list(cs$genes, context$gene_lists %||% list())
    if (length(symbols) == 0) abort("signatures: custom set '", cs$id, "' has no genes")
    out[[cs$id]] <- signatures_symbols_to_ids(symbols, features, cfg)
  }
  if (length(out) == 0) abort("signatures: no sets configured (modules.signatures.msigdb_sets or custom_sets)")
  out
}

# Custom sets are written as gene symbols; feature IDs come from the annotation when the IDs are Ensembl.
signatures_symbols_to_ids <- function(symbols, f, cfg) {
  if (identical(cfg$annotation$id_type, "gene_symbol") || is.null(f)) {
    return(unique(symbols))
  }
  unique(f$feature_id[f$gene_symbol %in% symbols])
}

# Scores of one method on matrix x (genes x samples) for sets (members within x) and composites; returns a
# data.frame with set_id, n_genes and one column per sample.
signatures_score <- function(x, members, composites, method, min_size) {
  ok_comp <- Filter(function(cp) all(c(cp$up, cp$down) %in% names(members)), composites)
  if (length(members) == 0) {
    return(data.frame())
  }
  rk <- if (method == "singscore") singscore::rankGenes(x)
  s <- switch(method,
    zscore = {
      z <- t(scale(t(x)))
      z[is.na(z)] <- 0
      t(vapply(members, function(g) colMeans(z[g, , drop = FALSE]), numeric(ncol(x))))
    },
    ssgsea = {
      r <- GSVA::gsva(GSVA::ssgseaParam(x, members, minSize = min_size), verbose = FALSE)
      r[names(members), , drop = FALSE]
    },
    singscore = {
      t(vapply(members, function(g) singscore::simpleScore(rk, upSet = g)$TotalScore, numeric(ncol(x))))
    }
  )
  s <- matrix(s, nrow = length(members), dimnames = list(names(members), colnames(x)))
  out <- data.frame(set_id = rownames(s), n_genes = lengths(members), s, check.names = FALSE, row.names = NULL)
  for (cp in ok_comp) {
    v <- if (method == "singscore") {
      singscore::simpleScore(rk, upSet = members[[cp$up]], downSet = members[[cp$down]])$TotalScore
    } else {
      s[cp$up, ] - s[cp$down, ]
    }
    row <- data.frame(set_id = cp$id, n_genes = NA_integer_, t(stats::setNames(v, colnames(x))), check.names = FALSE)
    out <- rbind(out, row)
  }
  out
}

signatures_long <- function(scores, samples) {
  empty <- data.frame(
    universe = character(), method = character(), set_id = character(), n_genes = integer(), sample_id = character(),
    score = numeric(), stringsAsFactors = FALSE
  )
  rows <- lapply(scores, function(w) {
    data.frame(
      universe = rep(w$universe, each = length(samples)), method = rep(w$method, each = length(samples)),
      set_id = rep(w$set_id, each = length(samples)), n_genes = rep(as.integer(w$n_genes), each = length(samples)),
      sample_id = rep(samples, times = nrow(w)), score = as.vector(t(as.matrix(w[, samples, drop = FALSE]))),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(list(empty), rows))
}

signatures_differences <- function(long, contrasts, samples) {
  empty <- data.frame(
    universe = character(), method = character(), set_id = character(), contrast_id = character(),
    difference = numeric(), n_numerator = integer(), n_denominator = integer(), stringsAsFactors = FALSE
  )
  if (nrow(long) == 0) {
    return(empty)
  }
  key <- paste(long$universe, long$method, long$set_id, sep = "\r")
  first <- !duplicated(key)
  rows <- lapply(contrasts, function(cn) {
    num <- contrast_samples(cn$numerator, samples)
    den <- contrast_samples(cn$denominator, samples)
    mean_of <- function(ids) {
      v <- tapply(long$score[long$sample_id %in% ids], factor(key[long$sample_id %in% ids], levels = unique(key)), mean)
      as.numeric(v)
    }
    data.frame(
      universe = long$universe[first], method = long$method[first], set_id = long$set_id[first],
      contrast_id = cn$id, difference = mean_of(num) - mean_of(den), n_numerator = length(num),
      n_denominator = length(den), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, c(list(empty), rows))
  rownames(out) <- NULL
  out
}

signatures_pairs <- function(diffs, pairs) {
  rows <- lapply(pairs, function(p) {
    a <- diffs[diffs$contrast_id == p$first, , drop = FALSE]
    b <- diffs[diffs$contrast_id == p$second, , drop = FALSE]
    i <- match(paste(a$universe, a$method, a$set_id), paste(b$universe, b$method, b$set_id))
    data.frame(
      universe = a$universe, method = a$method, set_id = a$set_id, pair_id = p$id, first = p$first,
      second = p$second, first_difference = a$difference, second_difference = b$difference[i],
      concordant = sign(a$difference) == sign(b$difference[i]), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Dot plot of signature scores per sample
#'
#' @param scores Rows of `signature_scores` for one universe and method.
#' @param samples,labels `character`. Sample ids and their labels.
#' @param colour,shape `character` or `NULL`. Sample groupings (in the order of `samples`) used as point colour
#'   and shape.
#' @param subtitle `character(1)`. Plot subtitle.
#' @param titles `character`. Legend titles of `colour` and `shape` (e.g. the factor names).
#' @return A `ggplot` object: one row per set, one point per sample.
#' @examples
#' s <- data.frame(set_id = rep(c("A", "B"), each = 2), sample_id = c("s1", "s2"), score = c(-1, 1, 0.5, -0.5))
#' plot_signature_scores(s, c("s1", "s2"), c("S1", "S2"), colour = c("x", "y"))
#' @export
plot_signature_scores <- function(scores, samples, labels, colour = NULL, shape = NULL, subtitle = NULL,
                                  titles = character()) {
  d <- scores
  i <- match(d$sample_id, samples)
  d$label <- labels[i]
  d$colour <- if (is.null(colour)) "all" else colour[i]
  d$shape <- if (is.null(shape)) "all" else shape[i]
  d$set_id <- factor(d$set_id, levels = rev(unique(d$set_id)))
  # Explicit legend order: without it ggplot2 may order the colour and shape legends differently between runs.
  point_aes <- if (is.null(shape)) {
    ggplot2::aes(colour = .data$colour)
  } else {
    ggplot2::aes(colour = .data$colour, shape = .data$shape)
  }
  ggplot2::ggplot(d, ggplot2::aes(x = .data$score, y = .data$set_id)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#c9c8c3", linewidth = 0.4) +
    ggplot2::geom_point(point_aes, size = 3, alpha = 0.9) +
    ggplot2::scale_colour_manual(
      values = palette_bulkrnaseq(unique(d$colour)), name = if (length(titles) >= 1) titles[[1]]
    ) +
    ggplot2::scale_shape_discrete(name = if (length(titles) >= 2) titles[[2]]) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1), shape = ggplot2::guide_legend(order = 2)) +
    ggplot2::labs(title = "Signature scores per sample", subtitle = subtitle, x = "Score", y = NULL) +
    theme_bulkrnaseq(base_size = 10) +
    ggplot2::theme(legend.position = if (is.null(colour) && is.null(shape)) "none" else "right")
}
