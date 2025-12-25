#!/usr/bin/env Rscript
################################################################################
# Script 2: Expression Analysis & Functional Enrichment
#
# Purpose: Perform condition-specific annotation, differential expression
#          analysis, GO enrichment, and visualization
#
# Required Input Files:
#   - annotation_table_with_GO.txt: Annotation table from Script 1
#   - sample2condition.txt: Sample metadata (columns: sample, path, group)
#   - salmon quant.sf files: In directories specified in sample2condition.txt
#
# Output:
#   - annotation_table_final.txt: Complete annotation with TPM & DEG results
#   - PCA_plot.png: PCA visualization
#   - topGO results and plots for each comparison
#
# Usage:
#   Rscript 02_expression_analysis.R
#
################################################################################

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(topGO)
  library(stringr)
})

message("================================================================================")
message("Fungal Annotation Workflow - Script 2: Expression Analysis & Enrichment")
message("================================================================================\n")

################################################################################
# CONFIGURATION - Modify these as needed
################################################################################

# Input files
annotation_file <- "annotation_table_with_GO.txt"
sample_metadata <- "sample2condition.txt"

# Output files
final_output <- "annotation_table_final.txt"

# Analysis parameters
tpm_threshold <- 1          # Minimum mean TPM for condition-specific annotation
cv_threshold <- 1           # Maximum CV for condition-specific annotation
min_mean_tpm <- 2          # Minimum mean TPM for final selection
deseq_padj_threshold <- 1e-9  # Adjusted p-value threshold for DEGs

# Comparison groups (modify based on your experiment)
# Format: list(contrast_vector, output_filename)
comparisons <- list(
  list(c("condition", "group2", "group1"), "deseq2-group2-vs-group1.txt"),
  list(c("condition", "group3", "group2"), "deseq2-group3-vs-group2.txt")
)

################################################################################
# PART 1: Load Data
################################################################################

message("Step 1: Loading data...")

# Load annotation table
if (!file.exists(annotation_file)) {
  stop("Error: Annotation file not found: ", annotation_file)
}
df <- read_tsv(annotation_file, show_col_types = FALSE)
message("  Loaded annotation table: ", nrow(df), " entries")

# Load sample metadata
if (!file.exists(sample_metadata)) {
  stop("Error: Sample metadata not found: ", sample_metadata)
}
s2c <- read.table(sample_metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
s2c$group <- gsub(" ", "_", s2c$group)
message("  Loaded sample metadata: ", nrow(s2c), " samples\n")

# Auto-generate colors for groups if not specified

message("Step 2: Importing Salmon quantification data...")

files <- s2c$path
names(files) <- s2c$sample

# Check if all salmon files exist
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Error: Missing salmon files:\n", paste(missing_files, collapse = "\n"))
}

txi <- tximport(files, type = "salmon", txOut = TRUE)

# Avoid division by zero
txi$length[txi$length == 0] <- 1

message("  Imported quantification for ", ncol(txi$counts), " samples")
message("  Total transcripts: ", nrow(txi$counts), "\n")

################################################################################
# PART 3: Condition-Specific Annotation
################################################################################

message("Step 3: Performing condition-specific annotation...")

# Convert to TPM data frame
tpm <- as.data.frame(txi$abundance)
tpm <- tibble::rownames_to_column(tpm, var = "Gene")

# Reshape to long format and join with metadata
tpm_long <- pivot_longer(tpm, -Gene, names_to = "sample", values_to = "TPM")
tpm_annotated <- left_join(tpm_long, s2c %>% select(sample, group), by = "sample")

# Calculate statistics per gene and condition
selected_genes <- tpm_annotated %>%
  group_by(Gene, group) %>%
  summarise(
    meanTPM = mean(TPM, na.rm = TRUE),
    sdTPM   = sd(TPM, na.rm = TRUE),
    CV      = ifelse(meanTPM > 0, sdTPM / meanTPM, NA),
    .groups = "drop"
  ) %>%
  group_by(Gene) %>%
  summarise(
    high_expr_groups = sum(meanTPM >= tpm_threshold),
    expressed_group = paste(group[meanTPM >= tpm_threshold], collapse = ";"),
    all_CV_below_threshold = all(CV < cv_threshold, na.rm = TRUE),
    meanTPM = max(meanTPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    high_expr_groups == 1,
    all_CV_below_threshold == TRUE,
    meanTPM >= min_mean_tpm
  ) %>%
  select(Gene, expressed_group, meanTPM)

message("  Identified ", nrow(selected_genes), " condition-specific genes")

# Add condition-specific annotation to main table
df <- left_join(df, selected_genes, by = c("Node" = "Gene"))

# Add TPM values for all samples
tpm_wide <- tpm %>% 
  rename_with(~ paste0("TPM_", .x), -Gene)

df <- left_join(df, tpm_wide, by = c("Node" = "Gene"))

message("  Added TPM values for all samples\n")

################################################################################
# PART 4: Differential Expression Analysis
################################################################################

message("Step 4: Performing differential expression analysis...")

# Create DESeq2 object
sampleTable <- data.frame(condition = s2c$group)
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

# Run DESeq2
dds <- DESeq(dds)

message("  DESeq2 analysis completed")

# Perform all comparisons
for (i in seq_along(comparisons)) {
  comp <- comparisons[[i]]
  contrast <- comp[[1]]
  outfile <- comp[[2]]
  
  comp_name <- paste(contrast[2], "vs", contrast[3], sep = "_")
  message("  Analyzing: ", comp_name)
  
  res <- results(dds, contrast = contrast)
  res_naomit <- na.omit(res)
  res_sorted <- res_naomit[order(res_naomit$padj), ]
  
  # Save full results
  write.table(res_sorted, outfile, sep = "\t", quote = FALSE)
  
  # Count significant DEGs
  n_degs <- sum(res_sorted$padj <= deseq_padj_threshold, na.rm = TRUE)
  message("    Found ", n_degs, " DEGs (padj <= ", deseq_padj_threshold, ")")
  
  # Add to main table
  res_df <- as.data.frame(res_sorted)
  res_df$Gene <- rownames(res_df)
  res_df <- res_df[, c("Gene", "log2FoldChange", "padj")]
  colnames(res_df) <- c("Gene",
                        paste0("log2FC_", comp_name),
                        paste0("padj_", comp_name))
  res_df <- distinct(res_df, Gene, .keep_all = TRUE)
  
  df <- left_join(df, res_df, by = c("Node" = "Gene"))
}

message("\n")

# Extract DEG sets for each comparison (for topGO later)
deg_sets <- list()
for (i in seq_along(comparisons)) {
  comp <- comparisons[[i]]
  contrast <- comp[[1]]
  comp_name <- paste(contrast[2], "vs", contrast[3], sep = "_")
  padj_col <- paste0("padj_", comp_name)
  
  if (padj_col %in% colnames(df)) {
    deg_set <- df %>%
      filter(!is.na(.data[[padj_col]]),
             .data[[padj_col]] <= deseq_padj_threshold) %>%
      pull(Node) %>%
      unique()
    
    deg_sets[[comp_name]] <- deg_set
    message("  DEGs in ", comp_name, ": ", length(deg_set))
  }
}

message("\n")

################################################################################
# PART 5: PCA Plot
################################################################################

message("Step 5: Creating PCA plot...")

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Use DESeq2's built-in plotPCA function with default settings
pca_plot <- plotPCA(vsd, intgroup = "condition")

ggsave("PCA_plot.png", plot = pca_plot, width = 8, height = 6, dpi = 300, bg = "white")
message("  PCA plot saved to: PCA_plot.png\n")

################################################################################
# PART 6: GO Enrichment Analysis with topGO
################################################################################

message("Step 5: Performing GO enrichment analysis...")

# Collect all GO columns
go_columns <- grep("^GO_", colnames(df), value = TRUE)

if (length(go_columns) == 0) {
  warning("  No GO columns found. Skipping enrichment analysis.")
} else {
  # Prepare gene2GO mapping
  gene2go_long <- df %>%
    pivot_longer(cols = all_of(go_columns), names_to = "source", values_to = "GO") %>%
    filter(!is.na(GO), GO != "") %>%
    mutate(GO_list = strsplit(GO, ";")) %>%
    unnest(GO_list) %>%
    mutate(GO_list = str_trim(GO_list)) %>%
    filter(GO_list != "", GO_list != "-") %>%
    distinct(Node, GO_list)
  
  all_genes <- unique(gene2go_long$Node)
  gene2GO_map <- split(gene2go_long$GO_list, gene2go_long$Node)
  gene2GO_map <- lapply(gene2GO_map, unique)
  
  message("  Prepared GO mapping for ", length(all_genes), " genes")
  
  # Function to perform topGO analysis
  perform_topgo <- function(sig_genes, all_genes, gene2GO_map, comparison_name) {
    message("  Analyzing: ", comparison_name)
    
    geneList <- factor(as.integer(all_genes %in% sig_genes))
    names(geneList) <- all_genes
    
    GOdata <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2GO_map)
    
    # Run different algorithms
    res_classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    res_elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
    res_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    
    # Save results for each algorithm
    for (algo in c("classic", "elim", "weight01")) {
      res_obj <- get(paste0("res_", algo))
      
      allRes <- GenTable(
        GOdata,
        Pvalue = res_obj,
        orderBy = "Pvalue",
        topNodes = 200,
        numChar = 1000
      )
      
      outfile <- paste0("topGO_", comparison_name, "_", algo, ".txt")
      write.table(allRes, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
      
      # Create dot plot
      create_topgo_dotplot(GOdata, res_obj, algo, comparison_name)
    }
    
    message("    Completed topGO analysis for ", comparison_name)
  }
  
  # Function to create topGO dot plot
  create_topgo_dotplot <- function(GOdata, result_obj, method, comparison_name,
                                   top_n = 15, label_wrap = 40) {
    tab <- GenTable(
      GOdata,
      Pvalue = result_obj,
      topNodes = max(200, top_n),
      numChar = 1000
    )
    
    plot_df <- tab %>%
      mutate(
        Pvalue_num = suppressWarnings(as.numeric(Pvalue)),
        Annotated = as.numeric(Annotated),
        Significant = as.numeric(Significant),
        GeneRatio = Significant / pmax(Annotated, 1),
        Term_wrapped = str_wrap(Term, width = label_wrap)
      ) %>%
      filter(is.finite(Pvalue_num)) %>%
      arrange(Pvalue_num) %>%
      slice_head(n = top_n) %>%
      mutate(Term_wrapped = factor(Term_wrapped, levels = rev(Term_wrapped)))
    
    if (nrow(plot_df) == 0) {
      message("      No significant terms found for ", method)
      return(invisible(NULL))
    }
    
    p <- ggplot(plot_df, aes(x = GeneRatio, y = Term_wrapped)) +
      geom_point(aes(size = Significant, color = -log10(Pvalue_num))) +
      scale_size_continuous(name = "Significant") +
      scale_color_viridis_c(name = expression(-log[10](p)), option = "D") +
      labs(
        x = "GeneRatio (Significant / Annotated)",
        y = NULL,
        title = paste0("GO Enrichment: ", comparison_name, " (", method, ")")
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_text(size = 11))
    
    outfile <- paste0("topGO_dotplot_", comparison_name, "_", method, ".png")
    ggsave(outfile, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
  }
  
  # Perform topGO for each DEG set
  for (set_name in names(deg_sets)) {
    if (length(deg_sets[[set_name]]) > 0) {
      sig_genes <- deg_sets[[set_name]][deg_sets[[set_name]] %in% all_genes]
      if (length(sig_genes) >= 5) {
        perform_topgo(sig_genes, all_genes, gene2GO_map, set_name)
      } else {
        message("  Skipping ", set_name, ": too few genes with GO annotations")
      }
    }
  }
}

message("\n")

################################################################################
# PART 7: Save Final Output
################################################################################

message("Step 6: Saving final annotation table...")

write_tsv(df, final_output, na = "")

message("  Final table saved to: ", final_output)
message("  Total entries: ", nrow(df))
message("  Total columns: ", ncol(df))

# Summary statistics
message("\nSummary Statistics:")
message("  Condition-specific genes: ", sum(!is.na(df$expressed_group)))
for (set_name in names(deg_sets)) {
  message("  DEGs in ", set_name, ": ", length(deg_sets[[set_name]]))
}

message("\n================================================================================")
message("Expression analysis and functional enrichment completed successfully!")
message("================================================================================\n")
