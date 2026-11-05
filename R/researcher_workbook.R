#' Write the researcher workbook
#'
#' One Excel workbook that collects the main results of a run for the person who asked for the analysis, instead
#' of one workbook per module: a README sheet (project, design, thresholds, how to read each sheet, a glossary of
#' every column), the methods text, then one or more sheets per result type, built from the tables the modules
#' wrote. Sheets whose module did not run are left out. Every column must have an entry in
#' `inst/extdata/column_glossary.yaml` (columns named after samples are described once in the README).
#'
#' Sheets (`export.researcher_workbook.sheets`): `de` (per contrast with a significance metric: all tested genes,
#' then the UP and DOWN genes), `ora`, `geneset_test`, `markers`, `signatures`, `qc`, `pca`.
#'
#' @param cfg Config from [read_config()] (defaults filled).
#' @param out `character(1)`. Project output folder.
#' @return Path of the workbook, invisibly.
#' @examples
#' \donttest{
#' # called by the report step when export.researcher_workbook.enabled is true
#' }
#' @export
write_researcher_workbook <- function(cfg, out) {
  wb_cfg <- cfg$export$researcher_workbook
  caps <- read_capabilities(out)
  ran <- caps$module[caps$status == "run"]
  dir_of <- stats::setNames(caps$dir, caps$module)
  read_tab <- function(module, table) {
    if (!module %in% ran) {
      return(NULL)
    }
    f <- file.path(out, dir_of[[module]], "tables", paste0(table, ".tsv"))
    if (!file.exists(f)) f <- paste0(f, ".gz")
    if (!file.exists(f)) {
      return(NULL)
    }
    as.data.frame(data.table::fread(f, sep = "\t", na.strings = ""))
  }
  samples <- workbook_samples(cfg, out)
  acc <- new.env()
  acc$sheets <- list()
  acc$guide <- list()
  add <- function(name, df, what) {
    if (is.null(df) || nrow(df) == 0) {
      return(invisible())
    }
    name <- workbook_sheet_name(name, names(acc$sheets))
    acc$sheets[[name]] <- df
    acc$guide[[name]] <- what
  }
  wanted <- unlist(wb_cfg$sheets)
  thresholds <- with_defaults(cfg$thresholds, "thresholds")
  if ("de" %in% wanted) {
    ct <- workbook_contrast_stats(out, caps)
    for (key in unique(paste(ct$contrast_id, ct$method, sep = "\r"))) {
      x <- ct[paste(ct$contrast_id, ct$method, sep = "\r") == key, , drop = FALSE]
      tab <- workbook_de_table(x, thresholds)
      id <- x$contrast_id[[1]]
      src <- paste0(x$method[[1]], ", ", x$source[[1]])
      add(paste0("DE ", id), tab, paste0("All genes tested in contrast ", id, " (", src, "), most significant first."))
      add(paste0("UP ", id), tab[tab$class == "UP", , drop = FALSE], paste0("Genes UP in ", id, "."))
      add(paste0("DOWN ", id), tab[tab$class == "DOWN", , drop = FALSE], paste0("Genes DOWN in ", id, "."))
    }
  }
  if ("ora" %in% wanted) {
    tabs <- c("ora_go_bp", "ora_go_cc", "ora_go_mf", "ora_kegg", "ora_reactome")
    ora <- do.call(rbind, lapply(tabs, function(t) {
      x <- read_tab("ora", t)
      if (!is.null(x) && nrow(x) > 0) x
    }))
    if (!is.null(ora)) ora <- ora[order(ora$contrast_id, ora$database, ora$padj, method = "radix"), , drop = FALSE]
    add("ORA", workbook_select(ora, c(
      "contrast_id", "method", "database", "term_id", "description", "gene_ratio", "bg_ratio", "pvalue", "padj",
      "count", "gene_symbols", "inference_scope"
    )), "Over-represented terms among the significant genes of each contrast, all databases, by padj.")
  }
  if ("geneset_test" %in% wanted) {
    gs <- read_tab("geneset_test", "geneset_results")
    if (!is.null(gs)) {
      alpha <- cfg$modules$geneset_test$alpha %||% 0.05
      gs <- gs[order(gs$contrast_id, gs$method, gs$padj, method = "radix"), , drop = FALSE]
      add("Gene sets", workbook_select(gs, c(
        "contrast_id", "method", "inter_gene_cor", "use_ranks", "universe", "collection", "set_id", "size",
        "nes", "direction", "pvalue", "padj", "leading_edge", "inference_scope"
      )), paste0("Gene set tests, all sets (padj < ", alpha, " is significant), sorted by padj per contrast."))
    }
  }
  if ("markers" %in% wanted) {
    me <- read_tab("markers", "marker_expression")
    if (!is.null(me)) {
      add(
        "Markers CPM", workbook_wide(me, c("panel_id", "gene_symbol", "feature_id"), "cpm", samples),
        "Marker genes: counts per million (library size) per sample, one column per sample."
      )
    }
    ps <- read_tab("markers", "marker_panel_scores")
    if (!is.null(ps)) {
      add(
        "Marker panels", workbook_wide(ps, c("panel_id", "n_genes"), "cpm_sum", samples),
        "Sum of the CPM of each panel's genes, per sample."
      )
    }
    add("Marker ratios", read_tab("markers", "marker_ratios"), "log2 ratios between chosen genes and samples.")
  }
  if ("signatures" %in% wanted) {
    sc <- read_tab("signatures", "signature_scores")
    if (!is.null(sc)) {
      add(
        "Signature scores", workbook_wide(sc, c("universe", "method", "set_id", "n_genes"), "score", samples),
        "Signature scores per sample (one column per sample), per universe and method."
      )
    }
    add(
      "Signature differences", read_tab("signatures", "signature_differences"),
      "Mean score of the numerator samples minus the denominator samples, per contrast (descriptive)."
    )
  }
  if ("qc" %in% wanted) {
    add("QC", workbook_transpose(read_tab("qc", "qc_summary"), samples), "QC metrics (rows) per sample.")
  }
  if ("pca" %in% wanted) add("PCA", read_tab("pca", "pca_scores"), "PCA scores of the samples (main gene filter).")
  sheets <- acc$sheets
  guide <- acc$guide
  glossary <- workbook_glossary()
  sample_cols <- unique(c(samples$sample_id, samples$label, names(cfg$samples$factors), cfg$samples$label_column))
  used <- unique(unlist(lapply(sheets, names)))
  undocumented <- setdiff(used, c(glossary$column, sample_cols))
  if (length(undocumented) > 0) {
    abort("researcher workbook: column(s) without a glossary entry: ", toString(undocumented))
  }
  methods_path <- file.path(out, "report", "methods.md")
  methods <- if (file.exists(methods_path)) data.frame(methods = readLines(methods_path, warn = FALSE))
  readme <- workbook_readme(cfg, thresholds, guide, glossary[glossary$column %in% used, , drop = FALSE], samples)
  all_sheets <- c(list(README = readme), if (!is.null(methods)) list(Methods = methods), sheets)
  path <- file.path(out, "report", wb_cfg$file)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  workbook_write(all_sheets, path)
  invisible(path)
}

# Samples of the sample sheet (all of them, including QC-only samples) with their labels.
workbook_samples <- function(cfg, out) {
  sheet <- cfg$samples$sheet
  if (is.null(sheet)) {
    return(data.frame(sample_id = character(), label = character()))
  }
  ss <- read_sample_sheet(config_path(cfg, sheet), cfg$samples$factors %||% list())
  col <- cfg$samples$label_column
  label <- if (!is.null(col) && col %in% names(ss)) as.character(ss[[col]]) else ss$sample_id
  data.frame(sample_id = ss$sample_id, label = label, stringsAsFactors = FALSE)
}

# Contrast tables with a significance metric: imported (ingest) and computed by modules that provide them.
workbook_contrast_stats <- function(out, caps) {
  files <- c(
    file.path(out, "01_ingest", "objects", "contrast_stats.rds"),
    file.path(out, caps$dir[caps$status == "run"], "objects", "contrast_stats.rds")
  )
  # Fixed columns plus posterior_prob, the significance metric of NOISeq-type tables (an extension column).
  cols <- c(contrast_table_columns(), "posterior_prob")
  pieces <- lapply(files[file.exists(files)], function(f) {
    x <- readRDS(f)
    if (is.null(x$posterior_prob)) x$posterior_prob <- rep(NA_real_, nrow(x))
    x[, cols]
  })
  if (length(pieces) == 0) {
    return(NULL)
  }
  ct <- do.call(rbind, pieces)
  ct[ct$significance_metric != "none", , drop = FALSE]
}

workbook_de_table <- function(x, thresholds) {
  sig <- significance_values(x)
  out <- data.frame(
    gene_symbol = x$gene_symbol, feature_id = x$feature_id, base_mean = x$base_mean, lfc = x$lfc,
    fold_change = ifelse(is.na(x$lfc), NA_real_, sign(x$lfc) * 2^abs(x$lfc)), lfc_se = x$lfc_se, stat = x$stat,
    pvalue = x$pvalue, significance_metric = sig$metric, significance = sig$values,
    class = classify_de(x, thresholds), inference_scope = x$inference_scope, stringsAsFactors = FALSE
  )
  o <- if (sig$metric == "padj") order(out$significance, -abs(out$lfc)) else order(-out$significance, -abs(out$lfc))
  out <- out[o, , drop = FALSE]
  drop <- vapply(out, function(v) all(is.na(v)), TRUE)
  drop[c("gene_symbol", "feature_id", "lfc", "class")] <- FALSE
  out[, !drop, drop = FALSE]
}

# QC summary (one row per sample) as one row per metric and one column per sample label.
workbook_transpose <- function(df, samples) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  metrics <- setdiff(names(df), "sample_id")
  out <- data.frame(metric = metrics, stringsAsFactors = FALSE)
  labels <- samples$label[match(df$sample_id, samples$sample_id)]
  labels[is.na(labels)] <- df$sample_id[is.na(labels)]
  for (i in seq_len(nrow(df))) out[[labels[[i]]]] <- unlist(lapply(df[i, metrics], as.character))
  out
}

workbook_select <- function(df, cols) {
  if (is.null(df)) {
    return(NULL)
  }
  df[, intersect(cols, names(df)), drop = FALSE]
}

# Long table (one row per sample) to wide (one column per sample, named by the sample label).
workbook_wide <- function(df, id_cols, value, samples) {
  key <- do.call(paste, c(df[id_cols], sep = "\r"))
  first <- !duplicated(key)
  out <- df[first, id_cols, drop = FALSE]
  labels <- samples$label[match(unique(df$sample_id), samples$sample_id)]
  ids <- unique(df$sample_id)
  for (i in seq_along(ids)) {
    v <- df[df$sample_id == ids[[i]], , drop = FALSE]
    out[[labels[[i]] %||% ids[[i]]]] <- v[[value]][match(key[first], do.call(paste, c(v[id_cols], sep = "\r")))]
  }
  rownames(out) <- NULL
  out
}

workbook_sheet_name <- function(name, taken) {
  for (ch in c("[", "]", ":", "*", "?", "/", "\\")) name <- gsub(ch, "_", name, fixed = TRUE)
  base <- substr(name, 1, 31)
  i <- 1
  while (base %in% taken) {
    i <- i + 1
    suffix <- paste0("~", i)
    base <- paste0(substr(name, 1, 31 - nchar(suffix)), suffix)
  }
  base
}

workbook_glossary <- function() {
  g <- yaml::read_yaml(package_file("extdata", "column_glossary.yaml"))
  data.frame(column = names(g), description = unlist(g, use.names = FALSE), stringsAsFactors = FALSE)
}

workbook_sample_names <- function(samples) {
  same <- samples$label == samples$sample_id
  ifelse(same, samples$sample_id, paste0(samples$label, " (", samples$sample_id, ")"))
}

workbook_readme <- function(cfg, thresholds, guide, glossary, samples) {
  ph <- methods_phrases(cfg)
  rows <- list(
    c("Project", cfg$project$name),
    c("Generated by", paste0("bulkrnaseq ", utils::packageVersion("bulkrnaseq"))),
    c("Samples", paste(workbook_sample_names(samples), collapse = ", ")),
    if (!is.null(ph$factors)) c("Sample factors", ph$factors),
    if (!is.null(ph$designs)) c("Designs (DESeq2)", ph$designs),
    if (!is.null(ph$contrast_tables)) c("Imported DE tables", ph$contrast_tables),
    if (!is.null(ph$contrasts)) c("Contrasts", ph$contrasts),
    c("DE thresholds", ph$de_thresholds %||% ""),
    c(
      "inference_scope", paste(
        "sample: p-values from a test whose units are samples;",
        "conditional_on_contrast: genes are the units, for one given contrast;",
        "ranking_only: gene permutations, p-values only rank sets within one contrast."
      )
    ),
    c("Sample columns", "Columns named after a sample hold that sample's value (label as in the sample sheet)."),
    c("", ""), c("Sheet", "Content"),
    if (!is.null(cfg$export$methods) && isTRUE(cfg$export$methods)) c("Methods", "Methods text of the run."),
    NULL
  )
  rows <- c(Filter(Negate(is.null), rows), lapply(names(guide), function(n) c(n, guide[[n]])))
  rows <- c(rows, list(c("", ""), c("Column", "Meaning")), lapply(seq_len(nrow(glossary)), function(i) {
    c(glossary$column[[i]], glossary$description[[i]])
  }))
  data.frame(item = vapply(rows, `[[`, "", 1), value = vapply(rows, `[[`, "", 2), stringsAsFactors = FALSE)
}

# Deterministic workbook with a bold, frozen header row and filters on every data sheet.
workbook_write <- function(sheets, path) {
  wb <- openxlsx::createWorkbook(creator = "bulkrnaseq")
  bold <- openxlsx::createStyle(textDecoration = "bold")
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    df <- sheets[[nm]]
    header <- !nm %in% c("README", "Methods")
    openxlsx::writeData(wb, nm, df, colNames = header || nm == "README", headerStyle = bold)
    if (header && nrow(df) > 0) {
      openxlsx::freezePane(wb, nm, firstRow = TRUE)
      openxlsx::addFilter(wb, nm, rows = 1, cols = seq_len(ncol(df)))
    }
    openxlsx::setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = if (header) 14 else c(28, 120)[seq_len(ncol(df))])
  }
  tmp <- tempfile(fileext = ".xlsx")
  openxlsx::saveWorkbook(wb, tmp, overwrite = TRUE)
  normalise_zip(tmp, path)
  unlink(tmp)
  invisible(path)
}
