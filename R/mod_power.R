#' Power module
#'
#' Indicative sample sizes for a two-group test on negative-binomial counts (Hart et al. 2013, J Comput Biol,
#' PMID 23961961): `n = 2 (z[1 - alpha/2] + z[power])^2 (1/mu + CV^2) / log(FC)^2` per group, over a grid of
#' biological CV (BCV), mean counts, fold changes and alpha. The BCV of the data is estimated with edgeR on a
#' design of the config; with few residual degrees of freedom it is only a reference point, and the table says
#' how wide a single gene's interval is.
#'
#' @param inputs `list` with `se` (counts).
#' @param params `list`. `modules.power` from the config.
#' @param context `list` with `cfg`, `gene_lists`, `organism`, `seed`.
#' @return `list(tables, figures, meta)`; tables `bcv`, `sample_size_grid`.
#' @examples
#' n_hart(cv = 0.3, mu = 100, fold_change = 2, alpha = 0.01)
#' @export
run_power <- function(inputs, params, context) {
  params <- with_defaults(params, c("modules", "power"))
  cfg <- context$cfg
  se <- inputs$se
  designs <- cfg$designs %||% list()
  d <- if (is.null(params$design)) designs[[1]] else Find(function(x) identical(x$id, params$design), designs)
  if (is.null(d)) abort("power module: design '", params$design %||% "", "' not found in designs")
  counts <- SummarizedExperiment::assay(se, "counts")
  des <- stats::model.matrix(stats::as.formula(d$formula), data = as.data.frame(SummarizedExperiment::colData(se)))
  df <- ncol(counts) - qr(des)$rank
  if (df < 1) abort("power module: design '", d$id, "' leaves no residual degree of freedom")
  keep <- config_filter_genes(se, cfg, context)
  universes <- list(all = rownames(se)[keep])
  rs <- params$reference_scenario
  grid <- params$grid
  # One seed for the whole bootstrap: universes run in sequence on the same random stream.
  bcv <- withr::with_seed(context$seed %||% 1, do.call(rbind, lapply(names(universes), function(u) {
    bcv_estimate(counts[universes[[u]], , drop = FALSE], des, params$bootstrap_replicates, u)
  })))
  bcv$residual_df <- df
  bcv$gene_bcv_ci95_factor_low <- sqrt(df / stats::qchisq(0.975, df))
  bcv$gene_bcv_ci95_factor_high <- sqrt(df / stats::qchisq(0.025, df))
  bcv$n_hart_common_bcv <- n_hart(bcv$common_bcv, rs$mu, rs$fold_change, rs$alpha, grid$power)
  bcv$n_hart_q90_raw <- n_hart(bcv$tagwise_raw_prior0_q90, rs$mu, rs$fold_change, rs$alpha, grid$power)
  sizes <- expand.grid(
    alpha = unlist(grid$alpha), fold_change = unlist(grid$fold_change), mu = unlist(grid$mu), cv = unlist(grid$cv)
  )[, c("cv", "mu", "fold_change", "alpha")]
  sizes <- sizes[do.call(order, sizes), ]
  sizes$power <- grid$power
  sizes$n_per_group <- n_hart(sizes$cv, sizes$mu, sizes$fold_change, sizes$alpha, grid$power)
  labels <- unlist(grid$cv_labels)
  if (length(labels) > 0 && length(labels) != length(unlist(grid$cv))) {
    abort("power module: grid.cv_labels needs one label per grid.cv value")
  }
  sizes$scenario <- if (length(labels)) labels[match(sizes$cv, unlist(grid$cv))] else NA_character_
  rownames(sizes) <- NULL
  width <- sprintf("x%.3g to x%.3g", bcv$gene_bcv_ci95_factor_low[[1]], bcv$gene_bcv_ci95_factor_high[[1]])
  list(
    tables = list(bcv = bcv, sample_size_grid = sizes), figures = list(),
    meta = list(
      bcv = list(
        description = paste0(
          "BCV (sqrt of edgeR dispersion) on design '", d$id, "' (", d$formula, "), genes of gene_filter. ",
          "common_bcv_boot_*: 2.5% and 97.5% quantiles over ", params$bootstrap_replicates, " bootstrap samples of ",
          "genes; tagwise_raw_prior0_*: per-gene BCV without shrinkage."
        ),
        residual_df = df, seed = context$seed %||% 1, reference_scenario = rs,
        caveats = c(
          paste0(
            "With ", df, " residual df the 95% interval of a single gene's BCV spans ", width,
            " its estimate (gene_bcv_ci95_factor_*)."
          ),
          "The bootstrap resamples genes: it measures stability to the choice of genes, not biological variability.",
          "Sample sizes are indicative (Hart et al. 2013 approximation), not a power analysis of a planned design."
        )
      ),
      sample_size_grid = list(
        description = "n per group from the Hart et al. 2013 formula over the grid of modules.power.grid.",
        caveat = "cv values and their labels are assumptions of the config, not estimates."
      )
    )
  )
}

#' Sample size per group for a two-group negative-binomial test (Hart et al. 2013)
#'
#' @param cv Biological coefficient of variation (BCV).
#' @param mu Mean count of the gene.
#' @param fold_change Fold change to detect.
#' @param alpha Two-sided significance level.
#' @param power Power.
#' @return Integer sample size per group (rounded up).
#' @examples
#' n_hart(0.4, 20, 1.5, 0.001)
#' @export
n_hart <- function(cv, mu, fold_change, alpha, power = 0.8) {
  ceiling(2 * (stats::qnorm(1 - alpha / 2) + stats::qnorm(power))^2 * (1 / mu + cv^2) / log(fold_change)^2)
}

bcv_estimate <- function(counts, des, replicates, universe) {
  y <- edgeR::calcNormFactors(edgeR::DGEList(counts))
  y <- edgeR::estimateDisp(y, des)
  boot <- replicate(replicates, {
    i <- sample(nrow(y), replace = TRUE)
    sqrt(edgeR::estimateGLMCommonDisp(y[i, , keep.lib.sizes = TRUE], des)$common.dispersion)
  })
  raw <- sqrt(edgeR::estimateGLMTagwiseDisp(y, des, prior.df = 0)$tagwise.dispersion)
  data.frame(
    universe = universe, n_genes = nrow(y), common_bcv = sqrt(y$common.dispersion),
    common_bcv_boot_q025 = unname(stats::quantile(boot, 0.025)),
    common_bcv_boot_q975 = unname(stats::quantile(boot, 0.975)),
    median_tagwise_bcv = stats::median(sqrt(y$tagwise.dispersion)),
    tagwise_raw_prior0_q10 = unname(stats::quantile(raw, 0.1)), tagwise_raw_prior0_median = stats::median(raw),
    tagwise_raw_prior0_q90 = unname(stats::quantile(raw, 0.9)), prior_df = y$prior.df[[1]]
  )
}
