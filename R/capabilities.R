#' Data levels provided by the configured inputs
#'
#' @param cfg Config from [read_config()].
#' @return Named `list`: one entry per available level with a short description of its source; attribute
#'   `contrasts` is a `data.frame` (`contrast_id`, `significance_metric`, `source`) of imported contrasts.
#' @examples
#' input_levels(list(inputs = list(counts = list(adapter = "matrix", file = "x.xlsx"))))
#' @export
input_levels <- function(cfg) {
  inp <- cfg$inputs %||% list()
  lv <- list()
  if (!is.null(inp$counts)) lv$counts <- paste0(inp$counts$adapter, " (", inp$counts$files %||% inp$counts$file, ")")
  if (!is.null(inp$abundance)) {
    lv$abundance <- paste0(inp$abundance$adapter, " (", inp$abundance$files %||% inp$abundance$file, ")")
  }
  if (!is.null(inp$normalized)) lv$normalized <- paste0("matrix (", inp$normalized$file, ")")
  if (!is.null(inp$coverage)) lv$coverage <- paste0("bedgraph (", inp$coverage$files, ")")
  qc <- c(
    if (identical(inp$counts$adapter, "star_reads_per_gene")) "STAR N_* rows",
    if (!is.null(inp$qc$multiqc)) {
      paste0("MultiQC (", paste(vapply(inp$qc$multiqc, `[[`, "", "id"), collapse = ", "), ")")
    }
  )
  if (length(qc) > 0) lv$qc <- paste(qc, collapse = "; ")
  tabs <- inp$contrast_tables %||% list()
  contrasts <- data.frame(contrast_id = character(), significance_metric = character(), source = character())
  for (t in tabs) {
    metric <- switch(t$adapter,
      deseq2 = "padj",
      noiseq = "posterior_prob",
      generic = if (!is.null(t$columns$padj)) {
        "padj"
      } else if (!is.null(t$columns$posterior_prob)) {
        "posterior_prob"
      } else {
        "none"
      }
    )
    row <- data.frame(contrast_id = t$id, significance_metric = metric, source = basename(t$file))
    contrasts <- rbind(contrasts, row)
  }
  if (nrow(contrasts) > 0) {
    lv$contrast_stats <- paste0(nrow(contrasts), " imported table(s): ", paste(contrasts$source, collapse = ", "))
    lv$ranked <- "derived from imported contrast tables"
  }
  attr(lv, "contrasts") <- contrasts
  lv
}

#' Decide which modules run
#'
#' Starts from the levels provided by the inputs, adds the levels provided by modules whose requirements are met,
#' and repeats until nothing changes. Every module ends as `run`, `skipped` (with the reason) or `disabled`
#' (`modules.<name>.enabled: false`); nothing is dropped silently.
#'
#' @param cfg Config from [read_config()].
#' @param samples Sample sheet (`data.frame`) or `NULL`.
#' @param registry Module registry from [read_registry()].
#' @param gene_lists Named `list` of available gene lists (for list-type resources).
#' @return `data.frame` with `module`, `number`, `dir`, `status`, `reason`, `depends_on` (comma-separated
#'   step names: `ingest` and/or modules).
#' @examples
#' cfg <- list(
#'   inputs = list(contrast_tables = list(list(id = "a", adapter = "deseq2", file = "a.xlsx"))),
#'   resources = list(items = list())
#' )
#' resolve_capabilities(cfg)[, c("module", "status", "reason")]
#' @export
resolve_capabilities <- function(cfg, samples = NULL, registry = read_registry(), gene_lists = list()) {
  lv <- input_levels(cfg)
  imported <- attr(lv, "contrasts")
  status <- stats::setNames(rep(NA_character_, length(registry)), names(registry))
  reason <- stats::setNames(rep("", length(registry)), names(registry))
  for (nm in names(registry)) {
    if (isFALSE(cfg$modules[[nm]]$enabled)) {
      status[[nm]] <- "disabled"
      reason[[nm]] <- paste0(nm, " disabled: modules.", nm, ".enabled is false in the config")
    }
  }
  repeat {
    provided <- names(lv)
    for (nm in names(registry)[status %in% "run"]) provided <- union(provided, unlist(registry[[nm]][["provides"]]))
    contrasts <- available_contrasts(cfg, imported, status, samples)
    changed <- FALSE
    for (nm in names(registry)[is.na(status) | status == "skipped"]) {
      why <- module_unmet(registry[[nm]], cfg, provided, contrasts, samples, lv, gene_lists)
      new_status <- if (is.null(why)) "run" else "skipped"
      if (!identical(status[[nm]], new_status)) changed <- TRUE
      status[[nm]] <- new_status
      reason[[nm]] <- if (is.null(why)) "" else paste0(nm, " skipped: ", why)
    }
    if (!changed) break
  }
  depends <- vapply(names(registry), function(nm) {
    m <- registry[[nm]]
    used <- unique(c(unlist(m[["requires"]]), unlist(m[["requires_any"]]), unlist(m[["optional"]])))
    deps <- if (any(used %in% names(lv)) || length(used) == 0) "ingest" else character()
    for (other in names(registry)[status %in% "run"]) {
      if (other != nm && any(unlist(registry[[other]][["provides"]]) %in% used)) deps <- c(deps, other)
    }
    paste(unique(c("ingest", deps)), collapse = ",")
  }, "")
  data.frame(
    module = names(registry), number = vapply(registry, function(m) as.integer(m$number), 1L),
    dir = vapply(names(registry), module_dir_name, "", registry = registry), status = unname(status),
    reason = unname(reason), depends_on = ifelse(status == "run", depends, ""), stringsAsFactors = FALSE,
    row.names = NULL
  )
}

available_contrasts <- function(cfg, imported, status, samples) {
  out <- imported
  if (identical(status[["contrasts"]], "run")) {
    for (cn in cfg$contrasts %||% list()) {
      out <- rbind(out, data.frame(contrast_id = cn$id, significance_metric = "none", source = "contrasts module"))
    }
  }
  if (identical(status[["de_test"]], "run")) {
    for (d in cfg$designs %||% list()) {
      if (!is.null(d$contrast)) {
        out <- rbind(out, data.frame(contrast_id = d$id, significance_metric = "padj", source = "de_test module"))
      }
    }
  }
  out
}

module_unmet <- function(m, cfg, provided, contrasts, samples, lv, gene_lists) {
  if (!exists(m[["function"]], envir = asNamespace("bulkrnaseq"), inherits = FALSE)) {
    return(paste0("not implemented in bulkrnaseq ", utils::packageVersion("bulkrnaseq"), " (", m[["function"]], ")"))
  }
  req <- unlist(m[["requires"]])
  missing <- setdiff(req, provided)
  if (length(missing) > 0) {
    return(paste0(
      "needs data level(s) ", paste(missing, collapse = ", "), "; input provides ",
      if (length(names(lv))) paste(names(lv), collapse = ", ") else "none"
    ))
  }
  any_of <- unlist(m[["requires_any"]])
  if (length(any_of) > 0 && !any(any_of %in% provided)) {
    return(paste0(
      "needs one of ", paste(any_of, collapse = ", "), "; input provides ",
      if (length(names(lv))) paste(names(lv), collapse = ", ") else "none"
    ))
  }
  for (cond in m[["conditions"]] %||% list()) {
    why <- condition_unmet(cond, cfg, contrasts, samples, m$name)
    if (!is.null(why)) {
      return(why)
    }
  }
  for (r in unlist(m[["resources"]])) {
    key <- if (r == "annotation") cfg$annotation$gtf %||% "gtf" else r
    if (!key %in% names(cfg$resources$items %||% list()) && !r %in% names(gene_lists)) {
      return(paste0("needs resource '", key, "' (resources.items) or gene list '", r, "'"))
    }
  }
  keys <- param_resource_keys(with_defaults(cfg$modules[[m$name]], c("modules", m$name)))
  missing <- setdiff(keys, names(cfg$resources$items %||% list()))
  if (length(missing) > 0) {
    return(paste0(
      "needs resource(s) ", paste0("'", missing, "'", collapse = ", "), " named in modules.", m$name,
      " (resources.items)"
    ))
  }
  NULL
}

# Resource keys named by a module's parameters: the values of every field called `resource`, at any depth.
param_resource_keys <- function(params) {
  if (!is.list(params)) {
    return(character())
  }
  nms <- names(params) %||% rep("", length(params))
  own <- unlist(params[nms == "resource"])
  unique(as.character(c(own, unlist(lapply(params[nms != "resource"], param_resource_keys)))))
}

condition_unmet <- function(cond, cfg, contrasts, samples, module = NULL) {
  key <- names(cond)[[1]]
  val <- cond[[1]]
  n_samples <- if (is.null(samples)) 0L else nrow(samples_with_level(samples, "counts"))
  switch(key,
    min_samples = if (n_samples < val) paste0("needs >= ", val, " samples with counts; input provides ", n_samples),
    contrasts_defined = if (length(cfg$contrasts %||% list()) == 0) {
      "needs contrasts defined in the config (contrasts:)"
    },
    contrast_pairs_defined = if (length(cfg$contrast_pairs %||% list()) == 0) {
      "needs contrast_pairs defined in the config"
    },
    min_contrasts = {
      n <- nrow(contrasts)
      if (n < val) {
        paste0(
          "needs >= ", val, " contrasts with per-gene statistics; input provides ", n,
          if (n > 0) paste0(" (", paste(unique(contrasts$source), collapse = ", "), ")") else ""
        )
      }
    },
    significance_available = if (!any(contrasts$significance_metric != "none")) {
      paste0(
        "needs a contrast with a significance metric (padj or posterior_prob); available contrasts are ",
        if (nrow(contrasts)) "descriptive only (significance_metric = none)" else "none"
      )
    },
    design_df_min = {
      ok <- FALSE
      if (!is.null(samples)) {
        ss <- samples_with_level(samples, "counts")
        for (d in cfg$designs %||% list()) {
          if (residual_df(stats::as.formula(d$formula), ss) >= val) ok <- TRUE
        }
      }
      if (!ok) {
        paste0("needs a design with >= ", val, " residual degree(s) of freedom (n samples - rank of the design matrix)")
      }
    },
    configured = {
      params <- with_defaults(cfg$modules[[module]], c("modules", module))
      set <- vapply(unlist(val), function(path) {
        x <- params
        for (k in strsplit(path, ".", fixed = TRUE)[[1]]) x <- x[[k]]
        length(unlist(x)) > 0
      }, TRUE)
      if (!any(set)) paste0("needs ", paste0("modules.", module, ".", unlist(val), collapse = " or "), " in the config")
    },
    abort("modules.yaml: unknown condition '", key, "'")
  )
}
