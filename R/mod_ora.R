#' Map gene identifiers to Entrez IDs
#'
#' @param ids `character`. Gene symbols or Ensembl IDs.
#' @param id_type `character(1)`. `"gene_symbol"`, `"ensembl_gene_id"` or `"entrez_id"`.
#' @param organism `character(1)`. Organism key (the OrgDb package comes from the organism table).
#' @return `data.frame` with `feature_id` and `entrez_id` (first Entrez ID per input; `NA` if unmapped).
#' @examples
#' \donttest{
#' if (requireNamespace(organism_info("human")$orgdb, quietly = TRUE)) {
#'   map_to_entrez(c("TP53", "NOTAGENE"), "gene_symbol")
#' }
#' }
#' @export
map_to_entrez <- function(ids, id_type = c("gene_symbol", "ensembl_gene_id", "entrez_id"), organism = "human") {
  id_type <- match.arg(id_type)
  if (id_type == "entrez_id") {
    return(data.frame(feature_id = ids, entrez_id = ids, stringsAsFactors = FALSE))
  }
  orgdb <- load_orgdb(organism)
  key <- c(gene_symbol = "SYMBOL", ensembl_gene_id = "ENSEMBL")[[id_type]]
  entrez <- suppressMessages(AnnotationDbi::mapIds(orgdb,
    keys = unique(ids), keytype = key, column = "ENTREZID",
    multiVals = "first"
  ))
  data.frame(feature_id = ids, entrez_id = unname(entrez[ids]), stringsAsFactors = FALSE)
}

load_orgdb <- function(organism) {
  pkg <- organism_info(organism)$orgdb
  if (!requireNamespace(pkg, quietly = TRUE)) abort("OrgDb package '", pkg, "' is not installed")
  getExportedValue(pkg, pkg)
}

#' Over-representation analysis module
#'
#' Selects significant genes of each contrast with the contrast's own significance metric and a linear
#' fold-change threshold from the config, then tests GO (BP, CC, MF), KEGG (pathway table from the resource
#' cache) and Reactome against the universe of tested genes. Hypergeometric p-values treat genes as
#' independent units for a given contrast, so every result has `inference_scope = "conditional_on_contrast"`.
#'
#' @param inputs `list` with `contrast_stats`.
#' @param params `list`. `modules.ora` (`databases`, `p_adjust_method`, `pvalue_cutoff`, `show_category`).
#' @param context `list` with `cfg`, `organism`, and `resources` (paths; `kegg_pathways` needed for KEGG).
#' @return `list(tables, figures, meta)`; tables `ora_input_genes` and one `ora_<database>` per database.
#' @examples
#' \donttest{
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'   requireNamespace(organism_info("human")$orgdb, quietly = TRUE)) {
#'   ids <- c("TP53", "MDM2", "CDKN1A", "BAX", "GAPDH", "ACTB")
#'   ct <- new_contrast_table("c", ids, c(2, 2, 2, 2, 0, 0),
#'     padj = c(0.01, 0.01, 0.01, 0.01, 0.9, 0.9),
#'     significance_metric = "padj", method = "m", source = "imported:x",
#'     inference_scope = "sample"
#'   )
#'   res <- run_ora(
#'     list(contrast_stats = ct), list(databases = "GO_BP"),
#'     list(cfg = list(annotation = list(id_type = "gene_symbol")), organism = "human")
#'   )
#'   names(res$tables)
#' }
#' }
#' @export
run_ora <- function(inputs, params, context) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) abort("ora module needs the 'clusterProfiler' package")
  ct <- inputs$contrast_stats
  ct <- ct[ct$significance_metric != "none", , drop = FALSE]
  if (nrow(ct) == 0) abort("ora module: no contrast with a significance metric")
  params <- with_defaults(params, c("modules", "ora"))
  cfg <- context$cfg
  thresholds <- with_defaults(cfg$thresholds, "thresholds")
  id_type <- with_defaults(cfg$annotation, "annotation")$id_type
  dbs <- unlist(params$databases)
  cutoff <- params$pvalue_cutoff
  adj <- params$p_adjust_method
  org <- organism_info(context$organism)
  orgdb <- load_orgdb(context$organism)
  kegg <- NULL
  if ("KEGG" %in% dbs) {
    kp <- context$resources$kegg_pathways
    if (is.null(kp)) abort("ora module: KEGG needs the 'kegg_pathways' resource (run fetch_resources)")
    kegg <- as.data.frame(data.table::fread(kp, colClasses = "character"))
  }
  input_rows <- list()
  results <- stats::setNames(vector("list", length(dbs)), dbs)
  for (key in unique(paste(ct$contrast_id, ct$method, sep = "\r"))) {
    x <- ct[paste(ct$contrast_id, ct$method, sep = "\r") == key, , drop = FALSE]
    sig <- significance_values(x)
    tested <- !is.na(sig$values) & !is.na(x$lfc)
    x <- x[tested, , drop = FALSE]
    cls <- classify_de(x, thresholds)
    map <- map_to_entrez(x$feature_id, id_type, context$organism)
    input_rows[[key]] <- data.frame(
      contrast_id = x$contrast_id, method = x$method, feature_id = x$feature_id,
      entrez_id = map$entrez_id, class = cls, stringsAsFactors = FALSE
    )
    genes <- unique(stats::na.omit(map$entrez_id[cls != "NS"]))
    universe <- unique(stats::na.omit(map$entrez_id))
    for (db in dbs) {
      res <- ora_one(db, genes, universe, orgdb, org, kegg, cutoff, adj)
      if (!is.null(res) && nrow(res) > 0) {
        res <- data.frame(
          contrast_id = x$contrast_id[[1]], method = x$method[[1]], direction = "up_and_down",
          database = db, res, inference_scope = "conditional_on_contrast",
          exploratory = contrast_exploratory(x), stringsAsFactors = FALSE
        )
        results[[db]] <- rbind(results[[db]], res)
      }
    }
  }
  inp <- do.call(rbind, unname(input_rows))
  rownames(inp) <- NULL
  tables <- list(ora_input_genes = inp)
  figs <- list()
  for (db in dbs) {
    nm <- paste0("ora_", tolower(db))
    tables[[nm]] <- results[[db]] %||% empty_ora_table()
    if (nrow(tables[[nm]]) > 0) figs[[nm]] <- plot_ora_dots(tables[[nm]], params$show_category)
  }
  meta <- list(ora_input_genes = list(
    description = "Genes tested by ORA, with Entrez mapping and UP/DOWN/NS class.",
    mapping = list(n_tested = nrow(inp), n_unmapped = sum(is.na(inp$entrez_id))), thresholds = thresholds
  ))
  for (db in dbs) {
    meta[[paste0("ora_", tolower(db))]] <- list(
      description = paste("Over-representation of", db, "terms among significant genes (UP and DOWN together)."),
      inference_scope = "conditional_on_contrast", p_adjust_method = adj, pvalue_cutoff = cutoff,
      units = list(pvalue = "p-value", padj = "adjusted p-value")
    )
  }
  list(tables = tables, figures = figs, meta = meta)
}

ora_one <- function(db, genes, universe, orgdb, org, kegg, cutoff, adj) {
  if (length(genes) == 0) {
    return(NULL)
  }
  res <- switch(db,
    GO_BP = ,
    GO_CC = ,
    GO_MF = clusterProfiler::enrichGO(genes,
      OrgDb = orgdb, ont = sub("GO_", "", db),
      universe = universe, pAdjustMethod = adj,
      pvalueCutoff = cutoff, readable = TRUE
    ),
    KEGG = {
      r <- clusterProfiler::enricher(genes,
        universe = universe, pAdjustMethod = adj, pvalueCutoff = cutoff,
        TERM2GENE = kegg[, c("pathway_id", "entrez_id")],
        TERM2NAME = unique(kegg[, c("pathway_id", "pathway_name")])
      )
      if (is.null(r)) NULL else clusterProfiler::setReadable(r, OrgDb = orgdb, keyType = "ENTREZID")
    },
    Reactome = {
      if (!requireNamespace("ReactomePA", quietly = TRUE)) abort("ora module: Reactome needs the 'ReactomePA' package")
      ReactomePA::enrichPathway(genes,
        organism = org$reactome_organism, universe = universe, pAdjustMethod = adj,
        pvalueCutoff = cutoff, readable = TRUE
      )
    },
    abort("ora module: unknown database '", db, "'")
  )
  if (is.null(res)) {
    return(NULL)
  }
  d <- as.data.frame(res)
  if (nrow(d) == 0) {
    return(NULL)
  }
  data.frame(
    term_id = d$ID, description = d$Description, gene_ratio = d$GeneRatio, bg_ratio = d$BgRatio,
    pvalue = d$pvalue, padj = d$p.adjust, qvalue = d$qvalue, count = d$Count, gene_symbols = d$geneID,
    stringsAsFactors = FALSE
  )
}

empty_ora_table <- function() {
  data.frame(
    contrast_id = character(), method = character(), direction = character(), database = character(),
    term_id = character(), description = character(), gene_ratio = character(), bg_ratio = character(),
    pvalue = numeric(), padj = numeric(), qvalue = numeric(), count = integer(), gene_symbols = character(),
    inference_scope = character(), exploratory = logical(), stringsAsFactors = FALSE
  )
}

#' Dot plot of ORA results
#'
#' @param d `data.frame` from [run_ora()] for one database.
#' @param show_category `integer(1)`. Terms shown per contrast (lowest adjusted p-value first).
#' @return A `ggplot` object.
#' @examples
#' d <- data.frame(
#'   contrast_id = "c", description = c("t1", "t2"), gene_ratio = c("2/10", "3/10"),
#'   padj = c(0.01, 0.02), count = c(2, 3)
#' )
#' plot_ora_dots(d, 2)
#' @export
plot_ora_dots <- function(d, show_category = 15) {
  d$ratio <- vapply(strsplit(d$gene_ratio, "/", fixed = TRUE), function(z) as.numeric(z[[1]]) / as.numeric(z[[2]]), 1)
  top <- function(z) utils::head(z[order(z$padj, z$description), ], show_category)
  d <- do.call(rbind, lapply(split(d, d$contrast_id), top))
  d$label <- factor(d$description, levels = rev(unique(d$description[order(d$padj)])))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$ratio, y = .data$label, size = .data$count, colour = .data$padj)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = "#e34948", high = "#2a78d6", name = "adjusted p") +
    ggplot2::facet_wrap(~contrast_id) +
    ggplot2::labs(x = "gene ratio", y = NULL, size = "genes") +
    theme_bulkrnaseq()
}
