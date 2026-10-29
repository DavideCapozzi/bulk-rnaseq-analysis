#' Pairing module
#'
#' Expression clues on whether samples come from the same donors (or pools): sex genes (XIST, chrY genes) and a
#' "fingerprint" of genes whose expression depends strongly on genotype or common deletions, compared between
#' every pair of samples as the mean |log2(CPM + 1)| difference. Clues, not genotype evidence: the expression of
#' these genes also depends on cell composition and depth, and the definitive test needs SNPs from the reads.
#'
#' @param inputs `list` with `se` (counts; `tpm` optional).
#' @param params `list`. `modules.pairing` from the config.
#' @param context `list` with `cfg`, `gene_lists`, `organism`.
#' @return `list(tables, figures, meta)`; tables `sex_genes`, `sex_summary`, `fingerprint_genes`,
#'   `fingerprint_distances`, `pairing_summary` (with `group_column` and `across_factor`), `male_fraction` and
#'   `male_fraction_genes` (with `male_reference_samples`); figure `pairing_sex_fingerprint`.
#' @examples
#' \donttest{
#' # called by the workflow on the project's SummarizedExperiment
#' }
#' @export
run_pairing <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "pairing"))
  se <- inputs$se
  counts <- SummarizedExperiment::assay(se, "counts")
  cpm <- edgeR::cpm(counts)
  rd <- SummarizedExperiment::rowData(se)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  org <- organism_info(context$organism %||% "human")
  lists <- context$gene_lists
  # First feature with each symbol, as a symbol lookup in a gene table.
  first_id <- function(symbols) rownames(se)[match(symbols, rd$gene_symbol)]
  sex <- sex_gene_tables(
    se, cpm, resolve_gene_list(params$sex_genes, lists),
    resolve_gene_list(params$xist_genes, lists), org
  )
  fp <- fingerprint_distances(se, counts, cpm, params, context, first_id)
  if (is.null(sex) && is.null(fp$distances)) {
    abort(
      "pairing: neither the sex genes (XIST and chrY genes) nor the genes of any fingerprint variant are among ",
      "the features; set modules.pairing.enabled: false"
    )
  }
  # Notes on what could not be computed go to the meta of every table written.
  notes <- c(
    if (is.null(sex)) "sex tables omitted: no XIST or chrY sex gene among the features",
    if (length(fp$skipped)) paste0("fingerprint variant(s) without genes among the features: ", toString(fp$skipped))
  )
  tables <- list(
    sex_genes = sex$sex_genes, sex_summary = sex$sex_summary, fingerprint_genes = fp$genes,
    fingerprint_distances = fp$distances
  )
  meta <- list(
    sex_genes = list(description = "Sex genes on the X and Y chromosomes: counts, CPM (all genes) and TPM."),
    sex_summary = list(
      description = "XIST and chrY sex genes per sample; all_chr_y_genes_cpm sums every chrY gene.",
      units = list(xist_cpm = "CPM", chr_y_cpm_sum = "CPM (sum over the chrY sex genes)")
    ),
    fingerprint_distances = list(
      description = "Mean |log2(CPM + 1) difference| over the fingerprint genes, per variant and pair of samples.",
      variants = fp$variants,
      caveat = "Expression clues, not genotype evidence: cell composition and depth also move these genes."
    )
  )
  if (!is.null(params$group_column) && !is.null(fp$distances)) {
    tables$pairing_summary <- pairing_summary(fp$distances, cd, params$group_column, params$across_factor)
    meta$pairing_summary <- list(
      description = paste0(
        "Per variant: are the pairs of the same '", params$group_column, "' the nearest pairs across '",
        params$across_factor, "'? chance_at_most bounds the probability of that outcome when the grouping is ",
        "unrelated to the samples (1/k! for k groups)."
      ),
      interpretation = paste(
        "This orients the question to whoever holds the sample records (are these the same donors?); it does",
        "not resolve it. A concordant result is weak evidence (see chance_at_most); a sex mismatch within a",
        "hypothesised pair (XIST against chrY genes in sex_summary) is the strongest expression clue against it."
      )
    )
  }
  refs <- unlist(params$male_reference_samples)
  if (length(refs) > 0) {
    bad <- setdiff(refs, colnames(se))
    if (length(bad) > 0) abort("pairing: male_reference_samples not among the samples with counts: ", toString(bad))
    if (is.null(sex)) abort("pairing: male_reference_samples is set but no XIST or chrY sex gene is among the features")
    mf <- male_fraction(cpm, sex, params$sex_denominators, refs, lists, first_id)
    tables <- c(tables, mf)
    meta$male_fraction <- list(
      description = paste(
        "Share of male RNA: ratio (chrY genes, optionally over a denominator) relative to its mean in the",
        "male reference samples, capped at 1; assumes similar chrY expression across samples and donors."
      ),
      male_reference_samples = refs
    )
    meta$male_fraction_genes <- list(
      description = "Per chrY gene: CPM and ratio to the male reference mean (uncapped)."
    )
  } else {
    notes <- c(notes, "male fraction not estimated: modules.pairing.male_reference_samples is empty")
  }
  tables <- Filter(Negate(is.null), tables)
  for (nm in names(tables)) meta[[nm]]$notes <- notes
  fill <- if (!is.null(params$across_factor)) cd[[params$across_factor]]
  first <- if (length(fp$variants)) fp$genes[fp$genes$variant_id == fp$variants[[1]]$id, ]
  figures <- list(pairing_sex_fingerprint = plot_pairing(
    colnames(se), sex$sex_summary, first, sample_labels(se, context$cfg), fill
  ))
  list(tables = tables, figures = figures, meta = meta, figure_sizes = list(pairing_sex_fingerprint = c(11, 7)))
}

sex_gene_tables <- function(se, cpm, sex_symbols, xist_symbols, org) {
  rd <- SummarizedExperiment::rowData(se)
  counts <- SummarizedExperiment::assay(se, "counts")
  tpm <- if ("tpm" %in% SummarizedExperiment::assayNames(se)) SummarizedExperiment::assay(se, "tpm")
  on_xy <- rd$chromosome %in% c(org$chromosome_x, org$chromosome_y)
  ids <- rownames(se)[rd$gene_symbol %in% sex_symbols & on_xy]
  ids <- ids[order(rd[ids, "chromosome"], rd[ids, "gene_symbol"])]
  y_ids <- ids[rd[ids, "chromosome"] == org$chromosome_y]
  x_ids <- ids[rd[ids, "gene_symbol"] %in% xist_symbols]
  if (length(y_ids) == 0 || length(x_ids) == 0) {
    return(NULL)
  }
  long <- expand.grid(feature_id = ids, sample_id = colnames(se), stringsAsFactors = FALSE)
  long <- data.frame(
    long[, "feature_id", drop = FALSE],
    gene_symbol = rd[long$feature_id, "gene_symbol"], chromosome = rd[long$feature_id, "chromosome"],
    sample_id = long$sample_id, count = counts[cbind(long$feature_id, long$sample_id)],
    cpm = cpm[cbind(long$feature_id, long$sample_id)],
    tpm = if (!is.null(tpm)) tpm[cbind(long$feature_id, long$sample_id)] else NA_real_
  )
  sums <- function(m, i) colSums(m[i, , drop = FALSE])
  summ <- data.frame(
    sample_id = colnames(se), xist_count = sums(counts, x_ids), chr_y_count_sum = sums(counts, y_ids),
    xist_cpm = sums(cpm, x_ids), chr_y_cpm_sum = sums(cpm, y_ids)
  )
  if (!is.null(tpm)) {
    summ$xist_tpm <- sums(tpm, x_ids)
    summ$chr_y_tpm_sum <- sums(tpm, y_ids)
  }
  summ$all_chr_y_genes_cpm <- sums(cpm, rd$chromosome %in% org$chromosome_y)
  summ$log2_xist_over_chr_y_cpm <- log2((summ$xist_cpm + 1) / (summ$chr_y_cpm_sum + 1))
  if (!is.null(tpm)) summ$log2_xist_over_chr_y_tpm <- log2((summ$xist_tpm + 1) / (summ$chr_y_tpm_sum + 1))
  rownames(summ) <- NULL
  y_symbols <- sex_symbols[sex_symbols %in% rd$gene_symbol[rd$chromosome %in% org$chromosome_y]]
  list(sex_genes = long, sex_summary = summ, y_ids = y_ids, y_symbols = y_symbols)
}

fingerprint_distances <- function(se, counts, cpm, params, context, first_id) {
  pairs <- utils::combn(colnames(se), 2)
  rd <- SummarizedExperiment::rowData(se)
  genes <- list()
  dists <- list()
  variants <- list()
  skipped <- character()
  tmm <- NULL
  for (v in params$fingerprint_variants) {
    symbols <- resolve_gene_list(v$genes, context$gene_lists)
    ids <- first_id(symbols)
    ids <- ids[!is.na(ids)]
    ids <- ids[rowSums(counts[ids, , drop = FALSE]) >= v$min_total_count]
    vals <- cpm
    if (v$normalization == "cpm_tmm") {
      if (is.null(tmm)) {
        keep <- config_filter_genes(se, context$cfg, context, "keep")
        tmm <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts[keep, , drop = FALSE])))
      }
      vals <- tmm
      ids <- intersect(ids, rownames(tmm))
    }
    if (length(ids) == 0) {
      skipped <- c(skipped, v$id)
      next
    }
    variants[[v$id]] <- list(id = v$id, normalization = v$normalization, genes = unname(rd[ids, "gene_symbol"]))
    l <- log2(vals[ids, , drop = FALSE] + 1)
    d <- apply(pairs, 2, function(p) mean(abs(l[, p[[1]]] - l[, p[[2]]])))
    dists[[v$id]] <- data.frame(
      variant_id = v$id, n_genes = length(ids), sample_a = pairs[1, ], sample_b = pairs[2, ], distance = d,
      rank_within_variant = rank(d)
    )
    genes[[v$id]] <- do.call(rbind, lapply(ids, function(i) {
      data.frame(
        variant_id = v$id, feature_id = i, gene_symbol = rd[i, "gene_symbol"], sample_id = colnames(se),
        count = counts[i, ], cpm = vals[i, ]
      )
    }))
  }
  if (length(dists) == 0) {
    return(list(skipped = skipped))
  }
  distances <- do.call(rbind, unname(dists))
  group <- params$group_column
  if (!is.null(group)) {
    cd <- SummarizedExperiment::colData(se)
    if (!group %in% names(cd)) abort("pairing: group_column '", group, "' is not in the sample sheet")
    g <- stats::setNames(as.character(cd[[group]]), colnames(se))
    distances$same_group <- g[distances$sample_a] == g[distances$sample_b]
    distances$pair_type <- ifelse(distances$same_group, "same_group", "different_group")
    if (!is.null(params$across_factor)) {
      f <- stats::setNames(as.character(cd[[params$across_factor]]), colnames(se))
      same_level <- f[distances$sample_a] == f[distances$sample_b]
      distances$pair_type[!distances$same_group] <- ifelse(same_level, "different_group_same_level",
        "different_group_different_level"
      )[!distances$same_group]
    }
  }
  rownames(distances) <- NULL
  out_genes <- do.call(rbind, unname(genes))
  rownames(out_genes) <- NULL
  list(distances = distances, genes = out_genes, variants = unname(variants), skipped = skipped)
}

# Per variant: are the same-group pairs the k nearest pairs across the two levels of `across`? With k groups,
# each with one sample per level, at most 1 of the k! matchings of the two levels has that property, so a
# grouping unrelated to the samples shows it with probability at most 1/k!.
pairing_summary <- function(distances, cd, group, across) {
  if (is.null(across)) abort("pairing: group_column needs across_factor (the levels a donor spans)")
  if (!across %in% names(cd)) abort("pairing: across_factor '", across, "' is not in the sample sheet")
  tab <- table(cd[[group]], cd[[across]])
  balanced <- ncol(tab) == 2 && all(tab == 1)
  k <- nrow(tab)
  rows <- lapply(split(distances, distances$variant_id), function(d) {
    across_pairs <- d$pair_type %in% c("same_group", "different_group_different_level")
    within <- d$distance[d$same_group & across_pairs]
    other <- d$distance[!d$same_group & across_pairs]
    data.frame(
      variant_id = d$variant_id[[1]], n_groups = k, max_same_group_distance = max(within),
      min_other_group_distance_across = min(other),
      same_group_pairs_nearest_across = if (balanced) max(within) < min(other) else NA,
      chance_at_most = if (balanced) 1 / factorial(k) else NA_real_
    )
  })
  out <- do.call(rbind, rows[unique(distances$variant_id)])
  rownames(out) <- NULL
  out
}

male_fraction <- function(cpm, sex, denominators, refs, lists, first_id) {
  summ <- sex$sex_summary
  out <- data.frame(
    sample_id = summ$sample_id, xist_cpm = summ$xist_cpm, chr_y_cpm_sum = summ$chr_y_cpm_sum,
    reference_chr_y_cpm = mean(summ$chr_y_cpm_sum[summ$sample_id %in% refs])
  )
  rel <- function(r) pmin(1, r / mean(r[colnames(cpm) %in% refs]))
  sum_cpm <- function(ref) {
    ids <- first_id(resolve_gene_list(ref, lists))
    if (anyNA(ids)) abort("pairing: sex_denominators gene(s) not among the features: ", toString(unlist(ref)))
    colSums(cpm[ids, , drop = FALSE])
  }
  for (dn in denominators) {
    num <- if (is.null(dn$numerator)) colSums(cpm[sex$y_ids, , drop = FALSE]) else sum_cpm(dn$numerator)
    if (is.null(dn$denominator)) {
      out[[paste0("male_fraction_", dn$id)]] <- unname(rel(num))
    } else {
      den <- sum_cpm(dn$denominator)
      out[[paste0("male_fraction_", dn$id)]] <- unname(rel(num / den))
      out[[paste0("xist_over_", dn$id)]] <- unname(summ$xist_cpm / den)
    }
  }
  genes <- do.call(rbind, lapply(sex$y_symbols, function(g) {
    v <- cpm[first_id(g), ]
    data.frame(gene_symbol = g, sample_id = colnames(cpm), cpm = unname(v), male_fraction = unname(v / mean(v[refs])))
  }))
  list(male_fraction = out, male_fraction_genes = genes)
}

plot_pairing <- function(samples, summ, fp_genes, labels, fill) {
  lab <- stats::setNames(labels, samples)
  d <- rbind(
    if (!is.null(summ)) data.frame(sample_id = summ$sample_id, panel = "XIST (CPM)", value = summ$xist_cpm),
    if (!is.null(summ)) {
      data.frame(sample_id = summ$sample_id, panel = "chrY sex genes (CPM sum)", value = summ$chr_y_cpm_sum)
    },
    if (!is.null(fp_genes)) {
      data.frame(
        sample_id = fp_genes$sample_id, panel = paste0(fp_genes$gene_symbol, " (log2 CPM+1)"),
        value = log2(fp_genes$cpm + 1)
      )
    }
  )
  d$panel <- factor(d$panel, unique(d$panel))
  d$label <- factor(lab[d$sample_id], unique(labels))
  d$fill <- if (!is.null(fill)) fill[match(d$sample_id, samples)] else "all"
  ggplot2::ggplot(d, ggplot2::aes(.data$label, .data$value, fill = .data$fill)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::facet_wrap(~panel, scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL, title = "Sex genes and genotype-sensitive genes") +
    theme_bulkrnaseq() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
