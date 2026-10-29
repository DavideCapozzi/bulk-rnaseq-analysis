#' Contamination module
#'
#' Signals of contaminating cells per sample and a mask of contaminant genes. Profiles are gene lists (e.g. the
#' genes driving PC1, or pan-T genes) scored on CPM of all genes (library size only): mean log2(CPM + 1) and CPM
#' sum, optionally relative to the sample with the largest sum (`reference_profile`). Ratio checks give per-gene
#' log2 CPM ratios between two samples over a profile (similar ratios are compatible with one source at different
#' doses); excess checks relate the CPM of other genes above a baseline to a profile's CPM sum. The mask joins
#' gene sets of contaminating cell types (minus genes in sets of the target cells) and the genes of chosen
#' profiles; downstream modules exclude it with the subset token `no_masked`.
#'
#' @param inputs `list` with `se` (counts).
#' @param params `list`. `modules.contamination` from the config.
#' @param context `list` with `cfg`, `gene_lists`, `resources` (paths by key).
#' @return `list(tables, figures, meta, gene_flags)`; tables `contamination_scores`,
#'   `contamination_profile_genes`, `profile_ratio_consistency`, `profile_excess`, `contaminant_mask`; figure
#'   `contamination_profiles`; `gene_flags` with `masked`, `masked_by_sets`, `masked_by_profiles` for every feature.
#' @examples
#' \donttest{
#' # called by the workflow on the project's SummarizedExperiment
#' }
#' @export
run_contamination <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "contamination"))
  if (length(params$profiles) == 0 && length(params$mask$include) == 0) {
    abort("contamination: no profiles and no mask gene sets configured; set modules.contamination.enabled: false")
  }
  se <- inputs$se
  counts <- SummarizedExperiment::assay(se, "counts")
  cpm <- edgeR::cpm(counts)
  rd <- SummarizedExperiment::rowData(se)
  profiles <- stats::setNames(
    lapply(params$profiles, function(p) resolve_gene_list(p$genes, context$gene_lists)),
    vapply(params$profiles, `[[`, "", "id")
  )
  check_profile <- function(id, field) {
    if (!id %in% names(profiles)) {
      abort("contamination: ", field, " '", id, "' is not among the profile ids (", toString(names(profiles)), ")")
    }
  }
  check_samples <- function(ids, field) {
    bad <- setdiff(ids, colnames(se))
    if (length(bad) > 0) abort("contamination: ", field, " not among the samples with counts: ", toString(bad))
  }
  for (field in c("profiles", "ratio_checks", "excess_checks")) {
    ids <- vapply(params[[field]], `[[`, "", "id")
    if (anyDuplicated(ids)) {
      abort("contamination: duplicated ", field, " id(s): ", toString(unique(ids[duplicated(ids)])))
    }
  }
  if (!is.null(params$reference_profile)) check_profile(params$reference_profile, "reference_profile")
  # First feature with each symbol, as a symbol lookup in a gene table.
  first_id <- function(symbols) rownames(se)[match(symbols, rd$gene_symbol)]
  tables <- list()
  meta <- list()
  if (length(profiles) > 0) {
    tables$contamination_scores <- profile_scores(profiles, cpm, first_id, se, context$cfg, params$reference_profile)
    tables$contamination_profile_genes <- profile_genes(profiles, cpm, first_id)
    meta$contamination_scores <- list(
      description = paste(
        "Per profile and sample: number of genes found, mean log2(CPM + 1) and CPM sum over the profile's genes",
        "(CPM on all genes, library size only); rel_to_max = CPM sum / largest CPM sum, reference profile only."
      ),
      profiles = lapply(profiles, unname)
    )
    meta$contamination_profile_genes <- list(description = "CPM (all genes, library size only) of each profile gene.")
  }
  if (length(params$ratio_checks) > 0) {
    for (rc in params$ratio_checks) {
      check_profile(rc$profile, "ratio_checks profile")
      check_samples(c(rc$numerator_sample, rc$denominator_sample), paste0("ratio_checks '", rc$id, "' samples"))
    }
    tables$profile_ratio_consistency <- ratio_consistency(
      params$ratio_checks, profiles, cpm, first_id, params$ratio_pseudocount
    )
    meta$profile_ratio_consistency <- list(
      description = paste0(
        "log2((CPM numerator + ", params$ratio_pseudocount, ") / (CPM denominator + ", params$ratio_pseudocount,
        ")) per gene of a profile."
      ),
      interpretation = paste(
        "Similar ratios across the genes are compatible with one source at different doses, not proof of it."
      )
    )
  }
  if (length(params$excess_checks) > 0) {
    for (ec in params$excess_checks) {
      check_profile(ec$profile, "excess_checks profile")
      check_samples(c(unlist(ec$baseline_samples), unlist(ec$samples)), paste0("excess_checks '", ec$id, "' samples"))
    }
    tables$profile_excess <- profile_excess(params$excess_checks, profiles, cpm, first_id, context$gene_lists)
    meta$profile_excess <- list(
      description = paste(
        "Per gene: CPM above the mean of the baseline samples, divided by the profile's CPM sum, in two samples;",
        "ratio = excess in the first sample / excess in the second. The profile row gives the same ratio for the",
        "profile's CPM sum."
      ),
      interpretation = paste(
        "A gene ratio similar to the profile ratio is compatible with the gene coming from the profile's cells."
      )
    )
  }
  for (id in unlist(params$mask$add_profiles)) check_profile(id, "mask.add_profiles")
  mask <- contaminant_mask(se, params$mask, profiles, context)
  tables$contaminant_mask <- mask$table
  meta$contaminant_mask <- mask$meta
  figures <- list()
  if (!is.null(tables$contamination_scores)) {
    fct <- sample_factors(se, context$cfg)
    fill <- if (length(fct)) as.character(SummarizedExperiment::colData(se)[[fct[[1]]]])
    figures$contamination_profiles <- plot_contamination(
      tables$contamination_scores, colnames(se), sample_labels(se, context$cfg), fill
    )
  }
  list(
    tables = tables, figures = figures, meta = meta, gene_flags = mask$flags,
    figure_sizes = list(contamination_profiles = c(10, 4.2))
  )
}

profile_scores <- function(profiles, cpm, first_id, se, cfg, reference_profile) {
  labels <- sample_labels(se, cfg)
  out <- do.call(rbind, lapply(names(profiles), function(p) {
    ids <- first_id(profiles[[p]])
    ids <- ids[!is.na(ids)]
    x <- cpm[ids, , drop = FALSE]
    data.frame(
      profile_id = p, sample_id = colnames(cpm), label = labels, n_genes = length(ids),
      mean_log2_cpm_p1 = colMeans(log2(x + 1)), sum_cpm = colSums(x), rel_to_max = NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  if (!is.null(reference_profile)) {
    i <- out$profile_id == reference_profile
    out$rel_to_max[i] <- out$sum_cpm[i] / max(out$sum_cpm[i])
  }
  rownames(out) <- NULL
  out
}

profile_genes <- function(profiles, cpm, first_id) {
  out <- do.call(rbind, lapply(names(profiles), function(p) {
    ids <- first_id(profiles[[p]])
    g <- expand.grid(i = seq_along(ids), sample_id = colnames(cpm), stringsAsFactors = FALSE)
    found <- !is.na(ids[g$i])
    v <- rep(NA_real_, nrow(g))
    v[found] <- cpm[cbind(ids[g$i][found], g$sample_id[found])]
    data.frame(
      profile_id = p, gene_symbol = profiles[[p]][g$i], feature_id = ids[g$i], sample_id = g$sample_id, cpm = v,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

ratio_consistency <- function(checks, profiles, cpm, first_id, pseudocount) {
  out <- do.call(rbind, lapply(checks, function(rc) {
    sym <- profiles[[rc$profile]]
    ids <- first_id(sym)
    num <- cpm[ids, rc$numerator_sample]
    den <- cpm[ids, rc$denominator_sample]
    data.frame(
      check_id = rc$id, profile_id = rc$profile, gene_symbol = sym, feature_id = ids,
      numerator_sample = rc$numerator_sample, denominator_sample = rc$denominator_sample,
      log2_ratio = log2((num + pseudocount) / (den + pseudocount)), stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

profile_excess <- function(checks, profiles, cpm, first_id, gene_lists) {
  out <- lapply(checks, function(ec) {
    base_s <- unlist(ec$baseline_samples)
    s <- unlist(ec$samples)
    prof_ids <- first_id(profiles[[ec$profile]])
    prof_sum <- colSums(cpm[prof_ids[!is.na(prof_ids)], , drop = FALSE])
    sym <- resolve_gene_list(ec$genes, gene_lists)
    ids <- first_id(sym)
    v <- cpm[ids, , drop = FALSE]
    excess <- sweep(v, 1, rowMeans(v[, base_s, drop = FALSE]))
    genes <- data.frame(
      check_id = ec$id, row_type = "gene", gene_symbol = sym, feature_id = ids, profile_id = ec$profile,
      stringsAsFactors = FALSE
    )
    prof_excess <- prof_sum - mean(prof_sum[base_s])
    prof <- data.frame(
      check_id = ec$id, row_type = "profile", gene_symbol = NA_character_, feature_id = NA_character_,
      profile_id = ec$profile, stringsAsFactors = FALSE
    )
    wide <- rbind(v, prof_sum)
    colnames(wide) <- paste0("cpm_", colnames(cpm))
    per <- rbind(
      sweep(excess[, s, drop = FALSE], 2, prof_sum[s], "/"),
      matrix(NA_real_, 1, length(s))
    )
    colnames(per) <- paste0("excess_per_profile_", s)
    ratio <- c(excess[, s[[1]]] / excess[, s[[2]]], prof_excess[[s[[1]]]] / prof_excess[[s[[2]]]])
    data.frame(rbind(genes, prof), wide, per, ratio = unname(ratio), check.names = FALSE, row.names = NULL)
  })
  # Checks on different samples have different columns: bind on the union, missing columns as NA.
  cols <- unique(unlist(lapply(out, names)))
  out <- do.call(rbind, lapply(out, function(x) {
    x[setdiff(cols, names(x))] <- NA
    x[cols]
  }))
  rownames(out) <- NULL
  out
}

contaminant_mask <- function(se, mask, profiles, context) {
  rd <- SummarizedExperiment::rowData(se)
  ids <- rownames(se)
  by_symbol <- identical(context$cfg$annotation$id_type, "gene_symbol")
  groups <- c(mask$include, mask$exclude)
  flags <- data.frame(feature_id = ids, gene_symbol = rd$gene_symbol, stringsAsFactors = FALSE)
  set_meta <- list()
  if (length(groups) > 0) {
    if (is.null(mask$resource)) abort("contamination: mask.include/exclude need mask.resource (an msigdb_export)")
    path <- context$resources[[mask$resource]]
    if (is.null(path)) abort("contamination: resource '", mask$resource, "' is not available in the cache")
    gs <- read_gene_sets(path, unique(unlist(lapply(groups, `[[`, "sets"))))
    key <- if (by_symbol) gs$gene_symbol else gs$ensembl_gene
    target <- if (by_symbol) rd$gene_symbol else ids
    for (g in groups) flags[[paste0("in_", g$id)]] <- target %in% key[gs$gs_name %in% unlist(g$sets)]
    item <- context$cfg$resources$items[[mask$resource]]
    set_meta <- list(
      resource = mask$resource, resource_version = item$version %||% NA_character_,
      db_version = sort(unique(gs$db_version)), match_by = if (by_symbol) "gene_symbol" else "ensembl_gene",
      include = lapply(mask$include, function(g) list(id = g$id, sets = unlist(g$sets))),
      exclude = lapply(mask$exclude, function(g) list(id = g$id, sets = unlist(g$sets)))
    )
  }
  any_of <- function(gs) {
    if (length(gs) == 0) {
      return(rep(FALSE, length(ids)))
    }
    Reduce(`|`, lapply(gs, function(g) flags[[paste0("in_", g$id)]]))
  }
  inc <- any_of(mask$include)
  exc <- any_of(mask$exclude)
  add <- unique(unlist(profiles[unlist(mask$add_profiles)]))
  flags$observed_marker <- !is.na(rd$gene_symbol) & rd$gene_symbol %in% add
  flags$masked <- (inc & !exc) | flags$observed_marker
  flag_cols <- setdiff(names(flags), c("feature_id", "gene_symbol"))
  table <- flags[Reduce(`|`, flags[flag_cols]), , drop = FALSE]
  rownames(table) <- NULL
  meta <- c(
    list(
      description = paste(
        "Genes with any flag: masked = (in an include group and in no exclude group) or observed_marker",
        "(a gene of mask.add_profiles); genes not listed have every flag FALSE."
      ),
      add_profiles = unlist(mask$add_profiles), n_masked = sum(flags$masked),
      n_masked_by_sets = sum(inc & !exc), n_kept_by_exclude = sum(inc & exc),
      n_observed_markers = sum(flags$observed_marker),
      n_observed_markers_covered_by_sets = sum(flags$observed_marker & inc & !exc),
      caveat = paste(
        "The mask removes known contaminant markers only; it cannot remove contaminant signal from other genes",
        "or show that the remaining signal is intrinsic to the target cells. Profiles chosen from these samples",
        "make the mask partly data-derived."
      )
    ),
    set_meta
  )
  gene_flags <- data.frame(
    feature_id = ids, masked = flags$masked, masked_by_sets = inc & !exc, masked_by_profiles = flags$observed_marker,
    stringsAsFactors = FALSE
  )
  list(table = table, meta = meta, flags = gene_flags)
}

plot_contamination <- function(scores, samples, labels, fill) {
  d <- scores
  n <- tapply(d$n_genes, d$profile_id, `[`, 1)
  d$panel <- factor(
    sprintf("%s (%d genes)", d$profile_id, n[d$profile_id]),
    unique(sprintf("%s (%d genes)", d$profile_id, n[d$profile_id]))
  )
  lab <- stats::setNames(labels, samples)
  d$label <- factor(lab[d$sample_id], unique(labels))
  d$fill <- if (!is.null(fill)) fill[match(d$sample_id, samples)] else "all"
  ggplot2::ggplot(d, ggplot2::aes(.data$label, .data$sum_cpm, fill = .data$fill)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f", .data$sum_cpm)), vjust = -0.3, size = 3) +
    ggplot2::facet_wrap(~panel, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(x = NULL, y = "CPM (sum)", fill = NULL, title = "Profile signals per sample") +
    theme_bulkrnaseq() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
