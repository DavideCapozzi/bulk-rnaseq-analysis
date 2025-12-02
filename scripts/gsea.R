# 1. CONFIGURATION -------------------------------------------------------------

# Analysis Configuration parameters
CONFIG <- list(
  # Input/Output settings
  input_file      = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/P_vs_M/Comparison_P_vs_M.xlsx",
  sheet_name      = "DESeq2_results",
  output_dir      = "/home/davidec/projects/bulk-rnaseq-analysis/Results_Enrichment",
  output_pdf_name = "Enrichment_Analysis.pdf",
  output_xls_name = "Enrichment_Results.xlsx",
  
  # Statistical Thresholds
  padj_cutoff     = 0.05,
  log2fc_cutoff   = 1.0, 
  
  # Visualization parameters
  show_category   = 15,            # Number of top terms to plot
  wrap_width      = 55,            # Max characters before wrapping text
  organism_db     = "org.Hs.eg.db",
  kegg_organism   = "hsa"
)

# 2. LIBRARY LOADING -----------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse) 
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
  library(ggrepel)
  library(openxlsx)
  # 'viridis' removed to use standard Red/Blue scales
})

# Create output directory
if (!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)

# 3. DATA IMPORT AND PROCESSING ------------------------------------------------

message(">>> Importing data...")
raw_data <- read_excel(CONFIG$input_file, sheet = CONFIG$sheet_name)

if (!"Gene_ID" %in% colnames(raw_data)) stop("Error: 'Gene_ID' column missing.")

# Filter valid data and assign significance
clean_data <- raw_data %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  mutate(Significance = case_when(
    padj <= CONFIG$padj_cutoff & log2FoldChange >= CONFIG$log2fc_cutoff ~ "UP",
    padj <= CONFIG$padj_cutoff & log2FoldChange <= -CONFIG$log2fc_cutoff ~ "DOWN",
    TRUE ~ "NS"
  ))

degs_data <- clean_data %>% filter(Significance %in% c("UP", "DOWN"))

message(">>> DEGs Summary:")
print(table(clean_data$Significance))

# 4. ID MAPPING ----------------------------------------------------------------

message(">>> Mapping IDs...")
# Map DEGs
mapped_degs <- bitr(degs_data$Gene_ID, fromType="SYMBOL", toType="ENTREZID", OrgDb=CONFIG$organism_db)
target_genes <- mapped_degs$ENTREZID

# Map Universe (Background)
mapped_universe <- bitr(clean_data$Gene_ID, fromType="SYMBOL", toType="ENTREZID", OrgDb=CONFIG$organism_db)
universe_genes <- mapped_universe$ENTREZID

message(paste(">>> Mapped", length(target_genes), "DEGs and", length(universe_genes), "background genes."))

# 5. FUNCTIONAL ENRICHMENT -----------------------------------------------------

message(">>> Running Enrichment Analysis...")

# We store results in a named list for both plotting and Excel export
analysis_results <- list()

# 5.1 GO Enrichment (Looping BP, CC, MF)
for (ont in c("BP", "CC", "MF")) {
  message(paste("   ... GO:", ont))
  ego <- enrichGO(gene = target_genes, universe = universe_genes, OrgDb = CONFIG$organism_db,
                  ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego)) analysis_results[[paste0("GO_", ont)]] <- ego
}

# 5.2 KEGG
message("   ... KEGG")
ekegg <- enrichKEGG(gene = target_genes, universe = universe_genes, organism = CONFIG$kegg_organism, pvalueCutoff = 0.05)
if (!is.null(ekegg)) {
  ekegg <- setReadable(ekegg, OrgDb = CONFIG$organism_db, keyType = "ENTREZID")
  analysis_results[["KEGG"]] <- ekegg
}

# 5.3 Reactome
message("   ... Reactome")
ereact <- enrichPathway(gene = target_genes, universe = universe_genes, pvalueCutoff = 0.05, readable = TRUE)
if (!is.null(ereact)) analysis_results[["Reactome"]] <- ereact

# 6. VISUALIZATION -------------------------------------------------------------

message(">>> Generating Plots...")
pdf(file.path(CONFIG$output_dir, CONFIG$output_pdf_name), width = 10, height = 7)

# 6.1 Volcano Plot
# Standard Red/Blue for Fold Change (UP=Red, DOWN=Blue)
top_genes <- clean_data %>% filter(Significance != "NS") %>% arrange(padj) %>% head(15)

p_volc <- ggplot(clean_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c("DOWN" = "#377EB8", "NS" = "grey90", "UP" = "#E41A1C")) +
  geom_vline(xintercept = c(-CONFIG$log2fc_cutoff, CONFIG$log2fc_cutoff), linetype="dashed", linewidth=0.4) +
  geom_hline(yintercept = -log10(CONFIG$padj_cutoff), linetype="dashed", linewidth=0.4) +
  geom_text_repel(data = top_genes, aes(label = Gene_ID), box.padding = 0.4, max.overlaps = 20, fontface="italic") +
  theme_bw(base_size = 14) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P.adj") +
  theme(legend.position = "top", plot.title = element_text(face="bold", hjust = 0.5))

print(p_volc)

# 6.2 Enrichment Plots Helper 
plot_enrichment <- function(res_obj, title) {
  if (is.null(res_obj) || nrow(as.data.frame(res_obj)) == 0) return(NULL)
  
  # Funzione helper per wrappare il testo
  wrapper <- function(x) stringr::str_wrap(x, width = CONFIG$wrap_width)
  
  # A. BARPLOT
  # Rimosso scale_fill_gradient manuale per usare il default di enrichplot (come nello script di riferimento)
  p_bar <- barplot(res_obj, showCategory = CONFIG$show_category) +
    ggtitle(paste(title, "- Barplot")) +
    scale_y_discrete(labels = wrapper) +
    theme(plot.title = element_text(face="bold"))
  
  print(p_bar)
  
  # B. DOTPLOT
  # Rimosso scale_color_gradient manuale per usare il default di enrichplot (come nello script di riferimento)
  p_dot <- dotplot(res_obj, showCategory = CONFIG$show_category) +
    ggtitle(paste(title, "- Dotplot")) +
    scale_y_discrete(labels = wrapper) +
    theme(plot.title = element_text(face="bold"))
  
  print(p_dot)
}

# Iterate and plot
names(analysis_results) %>% walk(~ plot_enrichment(analysis_results[[.x]], .x))

dev.off()

# 7. EXPORT EXCEL --------------------------------------------------------------

message(">>> Exporting Excel...")
# Extract dataframes from the result objects
excel_list <- lapply(analysis_results, as.data.frame)
excel_list[["DEGs_Input"]] <- degs_data

write.xlsx(excel_list, file.path(CONFIG$output_dir, CONFIG$output_xls_name))

message(">>> Pipeline Finished. Saved in: ", CONFIG$output_dir)