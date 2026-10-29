#' Reference module
#'
#' External reference profiles (Human Protein Atlas: Monaco immune cells in nTPM, single cell types in nCPM) for
#' every gene, and a gene universe defined by them alone: `reference_universe` = the target cells express the
#' gene (`target >= min_target_expr`) and no contaminating cell type dominates it
#' (`max(contaminants) >= contam_min` and `>= contam_fold x (target + 1)` makes a gene `contam_sensitive`).
#' The universe does not include the contamination mask: downstream subsets combine the two explicitly
#' (`in_reference_universe+no_masked`). When the mask is available it is used only to report why candidate
#' genes are excluded (`universe_exclusions`).
#'
#' @param inputs `list` with `se` (counts) and optional `gene_flags` (from the contamination module).
#' @param params `list`. `modules.reference` from the config.
#' @param context `list` with `cfg`, `gene_lists`, `organism`, `resources` (paths by key).
#' @return `list(tables, figures, meta, gene_flags)`; tables `reference_profiles`, `reference_universe` (gene
#'   counts per criterion), `universe_exclusions` (with `exclusion_candidates`); `gene_flags` with
#'   `target_expressed`, `contam_sensitive`, `reference_universe` for every feature.
#' @examples
#' \donttest{
#' # called by the workflow on the project's SummarizedExperiment
#' }
#' @export
run_reference <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "reference"))
  se <- inputs$se
  rd <- SummarizedExperiment::rowData(se)
  if (length(params$profiles) == 0) abort("reference: modules.reference.profiles is empty")
  ids <- vapply(params$profiles, `[[`, "", "id")
  if (anyDuplicated(ids)) abort("reference: duplicated profile id(s): ", toString(unique(ids[duplicated(ids)])))
  if (is.null(params$target) || !params$target %in% ids) {
    abort("reference: target '", params$target %||% "", "' is not among the profile ids (", toString(ids), ")")
  }
  bad <- setdiff(unlist(params$contaminants), ids)
  if (length(bad) > 0) abort("reference: contaminants not among the profile ids: ", toString(bad))
  prof <- reference_matrix(params$profiles, rownames(se), context)
  target <- prof[, params$target]
  contam_max <- if (length(params$contaminants) > 0) {
    do.call(pmax, c(as.data.frame(prof[, unlist(params$contaminants), drop = FALSE]), na.rm = TRUE))
  } else {
    rep(NA_real_, nrow(prof))
  }
  target_expressed <- !is.na(target) & target >= params$min_target_expr
  contam_sensitive <- !is.na(contam_max) & contam_max >= params$contam_min &
    contam_max >= params$contam_fold * (ifelse(is.na(target), 0, target) + 1)
  universe <- target_expressed & !contam_sensitive
  profiles <- data.frame(
    feature_id = rownames(se), gene_symbol = rd$gene_symbol, prof, target_expressed = target_expressed,
    contam_max = contam_max, contam_sensitive = contam_sensitive, reference_universe = universe,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  rownames(profiles) <- NULL
  filtered <- config_filter_genes(se, context$cfg, context)
  mask <- inputs$gene_flags
  masked <- if (!is.null(mask) && "masked" %in% names(mask)) {
    v <- mask$masked[match(rownames(se), mask$feature_id)]
    v[is.na(v)] <- FALSE
    v
  }
  count_row <- function(set, keep) {
    data.frame(
      genes = set, n_genes = sum(keep), n_with_target_reference = sum(keep & !is.na(target)),
      n_target_expressed = sum(keep & target_expressed), n_contam_sensitive = sum(keep & contam_sensitive),
      n_reference_universe = sum(keep & universe),
      n_reference_universe_unmasked = if (is.null(masked)) NA_integer_ else sum(keep & universe & !masked)
    )
  }
  tables <- list(
    reference_profiles = profiles,
    reference_universe = rbind(count_row("all_features", rep(TRUE, nrow(se))), count_row("filtered", filtered))
  )
  units <- vapply(params$profiles, function(p) {
    c(hpa_immune_cell = "nTPM", hpa_single_cell = "nCPM")[[context$cfg$resources$items[[p$resource]]$kind]]
  }, "")
  meta <- list(
    reference_profiles = list(
      description = paste(
        "Reference expression per gene (NA: gene absent from the reference). target_expressed: target >=",
        params$min_target_expr, "; contam_sensitive: max(contaminants) >=", params$contam_min, "and >=",
        params$contam_fold, "x (target + 1); reference_universe: target_expressed and not contam_sensitive."
      ),
      profiles = lapply(params$profiles, function(p) {
        list(
          id = p$id, resource = p$resource, cell = p$cell, unit = units[[match(p$id, ids)]],
          resource_version = context$cfg$resources$items[[p$resource]]$version %||% NA_character_
        )
      }),
      target = params$target, contaminants = unlist(params$contaminants),
      caveat = "Profiles come from other platforms, tissues and units: thresholds are coarse filters."
    ),
    reference_universe = list(
      description = paste(
        "Genes per criterion over all features and over the genes kept by gene_filter;",
        "n_reference_universe_unmasked also excludes the contamination mask (NA when the mask is not available)."
      )
    )
  )
  cand <- resolve_gene_list(params$exclusion_candidates, context$gene_lists)
  if (length(cand) > 0) {
    tables$universe_exclusions <- universe_exclusions(cand, se, profiles, params$target, filtered, mask)
    meta$universe_exclusions <- list(
      description = paste(
        "Candidate genes (protein-coding features with the symbol) and the first criterion that excludes them:",
        "not_filtered, mask_gene_sets, mask_profiles, absent_in_reference, target_below_min, contam_sensitive;",
        "included = filtered, in reference_universe and not masked."
      ),
      mask_applied = !is.null(masked)
    )
  }
  list(
    tables = tables, figures = list(), meta = meta,
    gene_flags = data.frame(
      feature_id = rownames(se), target_expressed = target_expressed, contam_sensitive = contam_sensitive,
      reference_universe = universe, stringsAsFactors = FALSE
    )
  )
}

# Reference values, features x profile ids (NA where a gene is absent from the reference). Each resource is read
# once, for all its cell types.
reference_matrix <- function(profiles, feature_ids, context) {
  ids <- vapply(profiles, `[[`, "", "id")
  out <- matrix(NA_real_, length(feature_ids), length(profiles), dimnames = list(NULL, ids))
  for (key in unique(vapply(profiles, `[[`, "", "resource"))) {
    path <- context$resources[[key]]
    if (is.null(path)) abort("reference: resource '", key, "' is not available in the cache")
    kind <- context$cfg$resources$items[[key]]$kind %||% ""
    if (!kind %in% c("hpa_immune_cell", "hpa_single_cell")) {
      abort("reference: resource '", key, "' needs kind hpa_immune_cell or hpa_single_cell (resources.items)")
    }
    mine <- Filter(function(p) identical(p$resource, key), profiles)
    long <- read_reference_profiles(path, kind, unique(vapply(mine, `[[`, "", "cell")))
    for (p in mine) {
      x <- long[long$cell == p$cell, , drop = FALSE]
      out[, p$id] <- x$value[match(feature_ids, x$feature_id)]
    }
  }
  out
}

universe_exclusions <- function(cand, se, profiles, target, filtered, mask) {
  rd <- SummarizedExperiment::rowData(se)
  pc <- which(rd$biotype %in% "protein_coding")
  rows <- lapply(cand, function(s) {
    i <- pc[rd$gene_symbol[pc] %in% s]
    if (length(i) == 0) i <- NA_integer_
    data.frame(gene_symbol = s, row = i)
  })
  x <- do.call(rbind, rows)
  i <- x$row
  flag <- function(col) if (is.null(mask)) rep(NA, length(i)) else mask[[col]][match(rownames(se)[i], mask$feature_id)]
  out <- data.frame(
    feature_id = rownames(se)[i], gene_symbol = x$gene_symbol, in_filtered = !is.na(i) & filtered[i],
    masked = flag("masked"), masked_by_sets = flag("masked_by_sets"), masked_by_profiles = flag("masked_by_profiles"),
    target = profiles[[target]][i], contam_max = profiles$contam_max[i],
    target_expressed = profiles$target_expressed[i], contam_sensitive = profiles$contam_sensitive[i],
    reference_universe = profiles$reference_universe[i], stringsAsFactors = FALSE
  )
  names(out)[names(out) == "target"] <- target
  has <- function(v) !is.na(v) & v
  out$exclusion_reason <- ifelse(!out$in_filtered, "not_filtered",
    ifelse(has(out$reference_universe) & !has(out$masked), "included",
      ifelse(has(out$masked_by_sets), "mask_gene_sets",
        ifelse(has(out$masked), "mask_profiles",
          ifelse(is.na(out[[target]]), "absent_in_reference",
            ifelse(!has(out$target_expressed), "target_below_min",
              ifelse(has(out$contam_sensitive), "contam_sensitive", "other")
            )
          )
        )
      )
    )
  )
  out <- out[order(out$exclusion_reason, out$gene_symbol), , drop = FALSE]
  rownames(out) <- NULL
  out
}
