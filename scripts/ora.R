# ==============================================================================
# UNIVERSAL ENRICHMENT ANALYSIS PIPELINE
# Supports: DESeq2 and NOISeq outputs
# ==============================================================================

# 1. CONFIGURATION -------------------------------------------------------------

CONFIG <- list(
  # Analysis Mode: "DESeq2" or "NOISeq"
  analysis_mode   = "DESeq2",  
  
  # Input/Output settings
  input_file      = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/P_vs_M/Comparison_P_vs_M_deseq.xlsx",
  sheet_name      = "DESeq2_results",  # Update sheet name for NOISeq if needed
  output_dir      = "/home/davidec/projects/bulk-rnaseq-analysis/Results_Enrichment",
  output_pdf_name = "Enrichment_Analysis_deseq.pdf",
  output_xls_name = "Enrichment_Results_deseq.xlsx",
  
  # Statistical Thresholds
  # DESeq2 specific
  padj_cutoff     = 0.05,
  # NOISeq specific
  prob_cutoff     = 0.8,
  # Common threshold
  log2fc_cutoff   = 1.0, 
  
  # Visualization parameters
  show_category   = 15,
  wrap_width      = 55,
  organism_db     = "org.Hs.eg.db",
  kegg_organism   = "hsa"
)

# Validate analysis mode
if (!CONFIG$analysis_mode %in% c("DESeq2", "NOISeq")) {
  stop("ERROR: analysis_mode must be either 'DESeq2' or 'NOISeq'")
}

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
})

# Create output directory
if (!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)

# 3. DATA IMPORT AND PROCESSING ------------------------------------------------

message(paste(">>> Running in", CONFIG$analysis_mode, "mode..."))
message(">>> Importing data...")
raw_data <- read_excel(CONFIG$input_file, sheet = CONFIG$sheet_name)

# ==============================================================================
# MODE-SPECIFIC DATA PROCESSING
# ==============================================================================

if (CONFIG$analysis_mode == "DESeq2") {
  # ---- DESeq2 MODE ----
  message(">>> Processing DESeq2 data...")
  
  # Check required columns
  required_cols <- c("Gene_ID", "padj", "log2FoldChange")
  missing_cols <- setdiff(required_cols, colnames(raw_data))
  if (length(missing_cols) > 0) {
    stop(paste("ERROR: Missing required columns for DESeq2:", paste(missing_cols, collapse = ", ")))
  }
  
  # Filter and assign significance
  clean_data <- raw_data %>%
    filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
    mutate(
      Significance = case_when(
        padj <= CONFIG$padj_cutoff & log2FoldChange >= CONFIG$log2fc_cutoff ~ "UP",
        padj <= CONFIG$padj_cutoff & log2FoldChange <= -CONFIG$log2fc_cutoff ~ "DOWN",
        TRUE ~ "NS"
      ),
      plot_metric = padj  # For volcano plot
    )
  
  y_axis_label <- "-Log10 P.adj"
  y_axis_transform <- function(x) -log10(x)
  y_threshold <- -log10(CONFIG$padj_cutoff)
  
} else {
  # ---- NOISeq MODE ----
  message(">>> Processing NOISeq data...")
  
  # Check required columns (NOISeq typically has: GeneID, M, prob)
  if (!"GeneID" %in% colnames(raw_data)) {
    stop("ERROR: 'GeneID' column missing in NOISeq data")
  }
  if (!"M" %in% colnames(raw_data)) {
    stop("ERROR: 'M' column (log2FC) missing in NOISeq data")
  }
  if (!"prob" %in% colnames(raw_data)) {
    stop("ERROR: 'prob' column missing in NOISeq data")
  }
  
  # Rename and process NOISeq columns
  clean_data <- raw_data %>%
    dplyr::rename(
      Gene_ID = GeneID,
      log2FoldChange = M
    ) %>%
    dplyr::filter(!is.na(prob) & !is.na(log2FoldChange)) %>%
    dplyr::mutate(
      # NOISeq significance based on probability
      Significance = case_when(
        prob >= CONFIG$prob_cutoff & log2FoldChange >= CONFIG$log2fc_cutoff ~ "UP",
        prob >= CONFIG$prob_cutoff & log2FoldChange <= -CONFIG$log2fc_cutoff ~ "DOWN",
        TRUE ~ "NS"
      ),
      # Create synthetic padj-like metric for plotting
      # High prob (e.g., 0.99) -> low plot_metric (e.g., 0.01)
      # This makes the volcano plot work correctly
      plot_metric = 1 - prob + 1e-10,
      # Keep original prob for export
      padj = plot_metric
    )
  
  y_axis_label <- "-Log10 (1 - Probability)"
  y_axis_transform <- function(x) -log10(1 - x + 1e-10)
  y_threshold <- -log10(1 - CONFIG$prob_cutoff + 1e-10)
}

# ==============================================================================
# COMMON PROCESSING (WORKS FOR BOTH MODES)
# ==============================================================================

degs_data <- clean_data %>% filter(Significance %in% c("UP", "DOWN"))

message(">>> DEGs Summary:")
print(table(clean_data$Significance))
message(paste("   UP-regulated:", sum(clean_data$Significance == "UP")))
message(paste("   DOWN-regulated:", sum(clean_data$Significance == "DOWN")))
message(paste("   Not Significant:", sum(clean_data$Significance == "NS")))

# 4. ID MAPPING ----------------------------------------------------------------

message(">>> Mapping Gene IDs to ENTREZ...")

# Map DEGs
mapped_degs <- bitr(
  degs_data$Gene_ID, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = CONFIG$organism_db
)
target_genes <- mapped_degs$ENTREZID

# Map Universe (Background)
mapped_universe <- bitr(
  clean_data$Gene_ID, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = CONFIG$organism_db
)
universe_genes <- mapped_universe$ENTREZID

message(paste(">>> Mapped", length(target_genes), "DEGs and", length(universe_genes), "background genes."))

if (length(target_genes) == 0) {
  stop("ERROR: No genes mapped successfully. Check gene symbols.")
}

# 5. FUNCTIONAL ENRICHMENT -----------------------------------------------------

message(">>> Running Enrichment Analysis...")
analysis_results <- list()

# 5.1 GO Enrichment (BP, CC, MF)
for (ont in c("BP", "CC", "MF")) {
  message(paste("   ... GO:", ont))
  ego <- tryCatch({
    enrichGO(
      gene = target_genes, 
      universe = universe_genes, 
      OrgDb = CONFIG$organism_db,
      ont = ont, 
      pAdjustMethod = "BH", 
      pvalueCutoff = 0.05, 
      readable = TRUE
    )
  }, error = function(e) {
    message(paste("   WARNING: GO", ont, "enrichment failed:", e$message))
    return(NULL)
  })
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    analysis_results[[paste0("GO_", ont)]] <- ego
  }
}

# 5.2 KEGG
message("   ... KEGG")
ekegg <- tryCatch({
  enrichKEGG(
    gene = target_genes, 
    universe = universe_genes, 
    organism = CONFIG$kegg_organism, 
    pvalueCutoff = 0.05
  )
}, error = function(e) {
  message(paste("   WARNING: KEGG enrichment failed:", e$message))
  return(NULL)
})

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  ekegg <- setReadable(ekegg, OrgDb = CONFIG$organism_db, keyType = "ENTREZID")
  analysis_results[["KEGG"]] <- ekegg
}

# 5.3 Reactome
message("   ... Reactome")
ereact <- tryCatch({
  enrichPathway(
    gene = target_genes, 
    universe = universe_genes, 
    pvalueCutoff = 0.05, 
    readable = TRUE
  )
}, error = function(e) {
  message(paste("   WARNING: Reactome enrichment failed:", e$message))
  return(NULL)
})

if (!is.null(ereact) && nrow(as.data.frame(ereact)) > 0) {
  analysis_results[["Reactome"]] <- ereact
}

if (length(analysis_results) == 0) {
  message("WARNING: No significant enrichment results found.")
}

# 6. VISUALIZATION -------------------------------------------------------------

message(">>> Generating Plots...")
pdf(file.path(CONFIG$output_dir, CONFIG$output_pdf_name), width = 10, height = 7)

# 6.1 VOLCANO PLOT (Mode-adaptive)
top_genes <- clean_data %>% 
  filter(Significance != "NS") %>% 
  arrange(plot_metric) %>% 
  head(15)

# Calculate y-axis values based on mode
clean_data_plot <- clean_data %>%
  mutate(y_value = -log10(plot_metric))

top_genes_plot <- top_genes %>%
  mutate(y_value = -log10(plot_metric))

p_volc <- ggplot(clean_data_plot, aes(x = log2FoldChange, y = y_value)) +
  geom_point(aes(color = Significance), alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c("DOWN" = "#377EB8", "NS" = "grey90", "UP" = "#E41A1C")) +
  geom_vline(
    xintercept = c(-CONFIG$log2fc_cutoff, CONFIG$log2fc_cutoff), 
    linetype = "dashed", 
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = y_threshold, 
    linetype = "dashed", 
    linewidth = 0.4
  ) +
  geom_text_repel(
    data = top_genes_plot, 
    aes(label = Gene_ID), 
    box.padding = 0.4, 
    max.overlaps = 20, 
    fontface = "italic"
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = paste("Volcano Plot -", CONFIG$analysis_mode, "Analysis"), 
    x = "Log2 Fold Change", 
    y = y_axis_label
  ) +
  theme(
    legend.position = "top", 
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_volc)

# 6.2 ENRICHMENT PLOTS
plot_enrichment <- function(res_obj, title) {
  if (is.null(res_obj) || nrow(as.data.frame(res_obj)) == 0) return(NULL)
  
  wrapper <- function(x) stringr::str_wrap(x, width = CONFIG$wrap_width)
  
  # Barplot
  p_bar <- tryCatch({
    barplot(res_obj, showCategory = CONFIG$show_category) +
      ggtitle(paste(title, "- Barplot")) +
      scale_y_discrete(labels = wrapper) +
      theme(plot.title = element_text(face = "bold"))
  }, error = function(e) {
    message(paste("   WARNING: Barplot for", title, "failed"))
    return(NULL)
  })
  
  if (!is.null(p_bar)) print(p_bar)
  
  # Dotplot
  p_dot <- tryCatch({
    dotplot(res_obj, showCategory = CONFIG$show_category) +
      ggtitle(paste(title, "- Dotplot")) +
      scale_y_discrete(labels = wrapper) +
      theme(plot.title = element_text(face = "bold"))
  }, error = function(e) {
    message(paste("   WARNING: Dotplot for", title, "failed"))
    return(NULL)
  })
  
  if (!is.null(p_dot)) print(p_dot)
}

# Plot all enrichment results
names(analysis_results) %>% walk(~ plot_enrichment(analysis_results[[.x]], .x))

dev.off()
message(paste(">>> Plots saved to:", file.path(CONFIG$output_dir, CONFIG$output_pdf_name)))

# 7. EXPORT EXCEL --------------------------------------------------------------

message(">>> Exporting Results to Excel...")

# Prepare Excel export
excel_list <- list()

# Add metadata sheet
excel_list[["Analysis_Info"]] <- data.frame(
  Parameter = c("Analysis_Mode", "Input_File", "Sheet_Name", 
                if(CONFIG$analysis_mode == "DESeq2") "Padj_Cutoff" else "Prob_Cutoff",
                "Log2FC_Cutoff", "Total_Genes", "UP_Regulated", "DOWN_Regulated"),
  Value = c(CONFIG$analysis_mode, CONFIG$input_file, CONFIG$sheet_name,
            if(CONFIG$analysis_mode == "DESeq2") CONFIG$padj_cutoff else CONFIG$prob_cutoff,
            CONFIG$log2fc_cutoff, nrow(clean_data),
            sum(clean_data$Significance == "UP"),
            sum(clean_data$Significance == "DOWN"))
)

# Add DEGs
excel_list[["DEGs_Input"]] <- degs_data

# Add enrichment results
for (name in names(analysis_results)) {
  excel_list[[name]] <- as.data.frame(analysis_results[[name]])
}

# Write Excel file
write.xlsx(excel_list, file.path(CONFIG$output_dir, CONFIG$output_xls_name))
message(paste(">>> Excel file saved to:", file.path(CONFIG$output_dir, CONFIG$output_xls_name)))

# 8. SUMMARY -------------------------------------------------------------------

message("\n=== PIPELINE COMPLETED SUCCESSFULLY ===")
message(paste("Mode:", CONFIG$analysis_mode))
message(paste("Total DEGs:", nrow(degs_data)))
message(paste("Enrichment analyses performed:", length(analysis_results)))
message(paste("Output directory:", CONFIG$output_dir))
message("========================================\n")

