#' Read named gene lists
#'
#' Gene lists are data, not code: the package ships generic lists in `inst/extdata/gene_lists/`; a project can
#' add its own files (`resources.gene_lists`). Each file has columns `list`, `gene_symbol`, `source`.
#'
#' @param files `character`. Extra list files; package lists are always read first, and a project list with
#'   the same name replaces the package one.
#' @return Named `list` of character vectors, with attribute `sources` (list name -> source text).
#' @examples
#' names(read_gene_lists())
#' @export
read_gene_lists <- function(files = character()) {
  pkg <- list.files(package_file("extdata", "gene_lists"), pattern = "\\.tsv$", full.names = TRUE)
  out <- list()
  sources <- list()
  for (f in c(pkg, files)) {
    x <- as.data.frame(data.table::fread(f, sep = "\t", colClasses = "character"))
    assert_columns(x, c("list", "gene_symbol", "source"), paste0("gene list file ", f))
    # Files are read in order, so a project list replaces a package list with the same name.
    for (nm in unique(x$list)) {
      out[[nm]] <- x$gene_symbol[x$list == nm]
      sources[[nm]] <- paste(unique(x$source[x$list == nm]), collapse = "; ")
    }
  }
  attr(out, "sources") <- sources
  out
}

# A gene_list_ref in the config is either the name of a list or an inline vector of symbols.
resolve_gene_list <- function(ref, gene_lists) {
  if (is.null(ref) || length(ref) == 0) {
    return(character())
  }
  if (length(ref) == 1 && ref %in% names(gene_lists)) {
    return(gene_lists[[ref]])
  }
  as.character(unlist(ref))
}

#' Select analysable genes
#'
#' Keeps genes of the chosen biotypes, outside the mitochondrial chromosome, with at least `min_cpm` counts
#' per million in at least `min_samples` samples (CPM on the full library), and optionally removes
#' sex-chromosome genes so that donor sex does not drive downstream comparisons.
#'
#' @param se `SummarizedExperiment` with a `counts` assay and rowData `biotype`, `chromosome`, `gene_symbol`.
#' @param biotypes `character`. Biotypes to keep (`NULL` keeps all).
#' @param exclude_mitochondrial `logical(1)`.
#' @param min_cpm `numeric(1)`. CPM threshold.
#' @param min_samples `integer(1)`. Number of samples that must pass `min_cpm`.
#' @param sex_genes `character(1)`. `"keep"`, `"drop_xist_y"` (Y chromosome and X-inactivation transcripts)
#'   or `"drop_xy"` (X and Y chromosomes).
#' @param organism `character(1)`. Organism key (chromosome names come from the organism table).
#' @param x_inactivation_genes `character`. Symbols removed with `"drop_xist_y"`.
#' @return `logical` vector, one value per gene.
#' @examples
#' m <- matrix(c(100L, 0L, 50L, 80L, 1L, 60L), 3, dimnames = list(c("a", "b", "c"), c("s1", "s2")))
#' se <- SummarizedExperiment::SummarizedExperiment(list(counts = m),
#'   rowData = S4Vectors::DataFrame(
#'     feature_id = rownames(m), gene_symbol = c("A", "B", "XIST"),
#'     biotype = "protein_coding", chromosome = c("1", "2", "X")
#'   )
#' )
#' filter_genes(se,
#'   min_cpm = 2, min_samples = 2, sex_genes = "drop_xist_y", x_inactivation_genes = "XIST"
#' )
#' @export
filter_genes <- function(se, biotypes = c("protein_coding", "lncRNA"), exclude_mitochondrial = TRUE, min_cpm = 2,
                         min_samples = 2, sex_genes = c("keep", "drop_xist_y", "drop_xy"), organism = "human",
                         x_inactivation_genes = character()) {
  sex_genes <- match.arg(sex_genes)
  org <- organism_info(organism)
  rd <- SummarizedExperiment::rowData(se)
  m <- SummarizedExperiment::assay(se, "counts")
  keep_bt <- if (is.null(biotypes)) rep(TRUE, nrow(m)) else rd$biotype %in% biotypes
  if (exclude_mitochondrial) keep_bt <- keep_bt & rd$chromosome != org$chromosome_mito
  keep_ex <- rowSums(edgeR::cpm(m) >= min_cpm) >= min_samples
  keep_sx <- switch(sex_genes,
    keep = rep(TRUE, nrow(m)),
    drop_xist_y = !(rd$chromosome == org$chromosome_y | rd$gene_symbol %in% x_inactivation_genes),
    drop_xy = !(rd$chromosome %in% c(org$chromosome_x, org$chromosome_y))
  )
  out <- keep_bt & keep_ex & keep_sx
  out[is.na(out)] <- FALSE
  unname(out)
}

# filter_genes() with the project's gene_filter settings; `sex_genes` overrides gene_filter.sex_genes (variants).
config_filter_genes <- function(se, cfg, context, sex_genes = NULL) {
  gf <- with_defaults(cfg$gene_filter, "gene_filter")
  filter_genes(se,
    biotypes = unlist(gf$biotypes), exclude_mitochondrial = gf$exclude_mitochondrial, min_cpm = gf$min_cpm,
    min_samples = gf$min_samples, sex_genes = sex_genes %||% gf$sex_genes,
    organism = context$organism %||% "human",
    x_inactivation_genes = resolve_gene_list("x_inactivation", context$gene_lists)
  )
}

#' Log2 counts per million after library normalisation
#'
#' TMM (default) assumes that most genes are not differentially expressed. When one sample carries a large
#' share of RNA from another cell type this fails, so callers can pass normalisation factors computed on
#' reference genes instead (see [reference_norm_factors()]).
#'
#' @param counts Integer matrix, genes x samples.
#' @param norm_method `character(1)`. Method of `edgeR::calcNormFactors()` (`"TMM"`, `"RLE"`,
#'   `"upperquartile"`, `"none"`), ignored when `norm_factors` is given.
#' @param prior_count `numeric(1)`. Average count added before the log, scaled by library size.
#' @param norm_factors `numeric` or `NULL`. Normalisation factors, one per sample.
#' @param lib_size `numeric` or `NULL`. Library sizes (default: column sums).
#' @return Numeric matrix of log2 CPM (unit: `log2_cpm`).
#' @examples
#' m <- matrix(c(10L, 20L, 30L, 5L, 40L, 25L), 3)
#' log_cpm(m)
#' @export
log_cpm <- function(counts, norm_method = "TMM", prior_count = 2, norm_factors = NULL, lib_size = NULL) {
  y <- if (is.null(lib_size)) edgeR::DGEList(counts) else edgeR::DGEList(counts, lib.size = lib_size)
  if (is.null(norm_factors)) {
    y <- edgeR::calcNormFactors(y, method = norm_method)
  } else {
    y$samples$norm.factors <- norm_factors
  }
  edgeR::cpm(y, log = TRUE, prior.count = prior_count)
}

#' Normalisation factors from reference genes
#'
#' Scales each library so that the summed counts of reference genes (e.g. genes expressed by the cell type of
#' interest) are equal across samples after normalisation; factors are centred on their geometric mean, as
#' edgeR expects.
#'
#' @param counts Integer matrix, genes x samples (the library size is its column sum).
#' @param reference_rows `logical` or `integer` index of the reference genes.
#' @return `numeric` normalisation factors, one per sample.
#' @examples
#' m <- matrix(c(10L, 20L, 70L, 20L, 20L, 160L), 3)
#' reference_norm_factors(m, 1:2)
#' @export
reference_norm_factors <- function(counts, reference_rows) {
  lib <- colSums(counts)
  ref <- colSums(counts[reference_rows, , drop = FALSE])
  share <- ref / lib
  share / exp(mean(log(share)))
}
