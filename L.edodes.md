## Computational Environment

All analyses were conducted on a **MacBook Pro (Apple M2 Max)** running **macOS 15.5**.  
A **Docker-based Ubuntu environment** was used unless otherwise specified for individual steps.

## Prepare fastq files

```bash
# Download 57 RNA-seq samples of Lentinula edodes (H600) from SRA
prefetch --option-file SRR.txt

# Convert .sra files to .fastq (SRR21185407 to SRR21185463)
for i in $(seq 21185407 21185463); do
  fasterq-dump "SRR${i}"
done

# Merge all forward and reverse reads into single files
cat *_1.fastq > cat-L.edodes_1.fastq
cat *_2.fastq > cat-L.edodes_2.fastq

# Trimming reads
docker attach trimgalore-0.6.10
trim_galore --paired cat-L.edodes_1.fastq cat-L.edodes_2.fastq
```

## Transcriptome Assembly

```bash
# RNA-Seq assembly
SPAdes-3.15.5-Linux/bin/rnaspades.py \
  -1 cat-shiitake_1_val_1.fq \
  -2 cat-shiitake_2_val_2.fq \
  -o shiitake_rnaspades_trimmed

# Check number of assembled transcripts (92,304)
grep -c "^>" shiitake_rnaspades_trimmed/transcripts.fasta

# BUSCO assessment
busco \
  -i shiitake_rnaspades_trimmed/transcripts.fasta \
  -o busco-results \
  -l eukaryota_odb10 \
  -m transcriptome \
  -c 8

**Results:**
- C: 97.3% [S: 5.5%, D: 91.8%]
- F: 2.4%, M: 0.4%, n: 255
- 248 Complete BUSCOs  
- 14 Complete and single-copy  
- 234 Complete and duplicated  
- 6 Fragmented  
- 1 Missing

# rnaQUAST assessment
docker attach rnaquast
python /usr/local/share/rnaquast-2.2.3-0/rnaQUAST.py \
  --transcripts shiitake_rnaspades_trimmed/transcripts.fasta \
  --reference GCF_021015755.1_Lenedo1_genomic.fna \
  --gtf genomic.gtf \
  --busco fungi_odb10
```

## ORF Prediction with TransDecoder

```bash
# Run TransDecoder (executed on host system due to Docker error)
TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs \
  -t shiitake_rnaspades_trimmed/transcripts.fasta

TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict \
  -t shiitake_rnaspades_trimmed/transcripts.fasta

# Check number of predicted ORFs (227,580)
grep -c "^>" transcripts.fasta.transdecoder_dir/longest_orfs.pep
```
## ggsearch & InterProScan

```bash
# ggsearch
time ggsearch36 -Q -T 96 -d 1 -m 10 -E 0.1 transcripts.fasta.transdecoder_dir/longest_orfs.pep Homo_sapiens.GRCh38.pep.all.fa > ggsearch-shiitake-human.txt & #96,665
time ggsearch36 -Q -T 96 -d 1 -m 10 -E 0.1 transcripts.fasta.transdecoder_dir/longest_orfs.pep Mus_musculus.GRCm39.pep.all.fa > ggsearch-shiitake-mouse.txt & #87,169
time ggsearch36 -Q -T 96 -d 1 -m 10 -E 0.1 transcripts.fasta.transdecoder_dir/longest_orfs.pep Saccharomyces_cerevisiae.R64-1-1.pep.all.fa > ggsearch-shiitake-S.cerevisiae.txt & #68,747
time ggsearch36 -Q -T 96 -d 1 -m 10 -E 0.1 transcripts.fasta.transdecoder_dir/longest_orfs.pep fungidb-68.fasta > ggsearch-shiitake-fungidb.txt & #182,249
time ggsearch36 -Q -T 96 -d 1 -m 10 -E 0.1 transcripts.fasta.transdecoder_dir/longest_orfs.pep uniprot_sprot.fasta > ggsearch-shiitake-uniprot.txt & #124,213

# InterProScan
# Remove asterisks (*) at the end of amino acid sequences from the FASTA file
sed 's/\*//g' transcripts.fasta.transdecoder_dir/longest_orfs.pep > longest_orfs_cleaned.pep
# Run InterProScan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/interproscan-5.67-99.0-64-bit.tar.gz
tar -pxvzf interproscan-5.67-99.0-*-bit.tar.gz
python3 setup.py -f interproscan.properties
./interproscan.sh -i ../output/longest_orfs_cleaned.pep -f tsv #130,603
./interproscan.sh -i longest_orfs_cleaned.pep -f tsv --goterms
```
## annotation table
```
# parse
perl -nle 'print $1 if(/^>(\S+)/)' transcripts.fasta.transdecoder_dir/longest_orfs.pep > longest_orf-pid.txt
for file in ggsearch-*.txt; do cat "$file" | perl 15parseggsearch.pl > "${file%.txt}-parsing.txt"; done
# make annotation table
cat longest_orf-pid.txt | perl 15mkannotbl-for-fungi.pl ggsearch-shiitake-human-parsing.txt | perl 15mkannotbl.pl ggsearch-shiitake-mouse-parsing.txt | perl 15mkannotbl.pl ggsearch-shiitake-S.cerevisiae-parsing.txt | perl 15mkannotbl.pl ggsearch-shiitake-uniprot-parsing.txt | perl 15mkannotbl.pl ggsearch-shiitake-fungidb-parsing.txt > ggsearch_annotbl.txt

# add InterProScan
R
library(dplyr)
library(stringr)

interpro_raw <- read.table(
  "interpro_longest_orfs.tsv",
  header = FALSE, sep = "\t", stringsAsFactors = FALSE,
  fill = TRUE, quote = ""   
)

colnames(interpro_raw) <- c(
  "Node", "ID", "Length", "Database", "Accession", "Description",
  "Start", "End", "Evalue", "Status", "Date", "InterProID", "InterProDescription",
  "Extra1", "Extra2"
)

interpro_top1 <- interpro_raw %>%
  dplyr::filter(stringr::str_detect(Evalue, "^\\d+\\.?\\d*([eE][-+]?\\d+)?$")) %>%
  dplyr::mutate(Evalue = as.numeric(Evalue)) %>%
  dplyr::group_by(Node) %>%
  dplyr::slice_min(order_by = Evalue, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(Node, Database, Accession, Description, InterProID, InterProDescription)

ggsearch_data <- read.table(
  "ggsearch_annotbl.txt",
  sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = ""
)

colnames(ggsearch_data) <- c(
  "Node",
  "human","human.1","human.2",
  "mouse","mouse.1","mouse.2",
  "yeast","yeast.1","yeast.2",
  "uniprot","uniprot.1","uniprot.2",
  "fungidb","fungidb.1","fungidb.2"
)

# UniProt: "sp|P50090|KEL2_YEAST" → split into ID and NAME
valid_uniprot_rows <- grepl("^sp\\|[^|]+\\|[^|]+$", ggsearch_data$uniprot)
if (any(valid_uniprot_rows)) {
  split_matrix <- do.call(rbind, strsplit(ggsearch_data$uniprot[valid_uniprot_rows], "\\|"))
  ggsearch_data$uniprot.id[valid_uniprot_rows]   <- split_matrix[, 2]  # P50090
  ggsearch_data$uniprot.name[valid_uniprot_rows] <- split_matrix[, 3]  # KEL2_YEAST
  ggsearch_data$uniprot   <- ggsearch_data$uniprot.id
  ggsearch_data$uniprot.1 <- ggsearch_data$uniprot.name
  ggsearch_data$uniprot.id <- NULL
  ggsearch_data$uniprot.name <- NULL
}

ggsearch_data$fungidb.1 <- ggsearch_data$fungidb.2
ggsearch_data$fungidb.2 <- NULL

df <- ggsearch_data %>%
  left_join(interpro_top1, by = "Node")

# write.table(df, "ggsearch-interpro_annotbl.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
```
## add GO term
```
library(biomaRt)
library(readr)
# human / mouse
add_go_terms_ensembl_pep <- function(df, id_col, mart_dataset, species_colname) {
  mart <- useMart("ensembl", dataset = mart_dataset)
  df[[paste0("clean_", id_col)]] <- gsub("\\..*", "", df[[id_col]])
  ids <- unique(na.omit(df[[paste0("clean_", id_col)]]))
  
  if (length(ids) == 0) return(df)
  
  go_raw <- getBM(
    attributes = c("ensembl_peptide_id", "go_id"),
    filters = "ensembl_peptide_id",
    values = ids,
    mart = mart
  )
  
  go_summarized <- go_raw %>%
    dplyr::group_by(ensembl_peptide_id) %>%
    dplyr::summarise(!!paste0("GO_", species_colname) := paste(unique(go_id), collapse = ";")) %>%
    dplyr::rename(!!paste0("clean_", id_col) := ensembl_peptide_id)
  
  left_join(df, go_summarized, by = paste0("clean_", id_col))
}

df <- add_go_terms_ensembl_pep(df, id_col = "human", mart_dataset = "hsapiens_gene_ensembl", species_colname = "human")
df <- add_go_terms_ensembl_pep(df, id_col = "mouse", mart_dataset = "mmusculus_gene_ensembl", species_colname = "mouse")
df <- dplyr::select(df, -dplyr::starts_with("clean_"))

# yeast
add_go_terms_ensembl_gene <- function(df, id_col, mart_dataset, species_colname) {
  mart <- useMart("ensembl", dataset = mart_dataset)
  ids <- unique(na.omit(df[[id_col]]))
  
  if (length(ids) == 0) return(df)
  
  go_raw <- getBM(
    attributes = c("ensembl_gene_id", "go_id"),
    filters = "ensembl_gene_id",
    values = ids,
    mart = mart
  )
  
  go_summarized <- go_raw %>%
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::summarise(!!paste0("GO_", species_colname) := paste(unique(go_id), collapse = ";")) %>%
    dplyr::rename(!!id_col := ensembl_gene_id)
  
  left_join(df, go_summarized, by = id_col)
}

df <- add_go_terms_ensembl_gene(df, id_col = "yeast", mart_dataset = "scerevisiae_gene_ensembl", species_colname = "yeast")

# UniProt ID
add_go_terms_uniprot <- function(df, id_col, uniprot_go_file, species_colname = "uniprot") {
  go_df <- read.delim(uniprot_go_file, header = FALSE, col.names = c("uniprot_id", "go_id"))
  go_summarized <- go_df %>%
    dplyr::group_by(uniprot_id) %>%
    dplyr::summarise(!!paste0("GO_", species_colname) := paste(unique(go_id), collapse = ";"))
  
  df <- left_join(df, go_summarized, by = setNames("uniprot_id", id_col))
  return(df)
}

df <- add_go_terms_uniprot(df, id_col = "uniprot", uniprot_go_file = "uniprot_go_annotations.tsv")

# InterPro
interpro2go_url <- "https://current.geneontology.org/ontology/external2go/interpro2go"
raw  <- readr::read_lines(interpro2go_url)
lines <- raw[!startsWith(raw, "!")] 
pairs_list <- lapply(lines, function(x) {
  ipr <- str_extract_all(x, "InterPro:(IPR\\d+)")[[1]]
  go  <- str_extract_all(x, "(GO:\\d+)")[[1]]
  if (length(ipr) == 0 || length(go) == 0) return(NULL)
   expand.grid(IPR = unique(sub("^InterPro:", "", ipr)),
              GO  = unique(go),
              stringsAsFactors = FALSE)
})

map_df <- dplyr::bind_rows(pairs_list) %>%
  dplyr::distinct() %>%
  dplyr::rename(InterProID = IPR, GO_ID = GO)

split_interpro <- function(x) {
  x %>%
    str_replace_all("InterPro:", "") %>%        
    str_split("\\s*[,;|]\\s*|\\s+")             
}

df_go <- df %>%
  dplyr::mutate(.rowid = dplyr::row_number(),
                .IPR_list = split_interpro(InterProID)) %>%
  tidyr::unnest_longer(.IPR_list, values_to = "InterProID_expanded") %>%
  dplyr::mutate(InterProID_expanded = dplyr::if_else(is.na(InterProID_expanded), "", InterProID_expanded),
                InterProID_expanded = stringr::str_trim(InterProID_expanded)) %>%
  dplyr::filter(InterProID_expanded != "", InterProID_expanded != "-") %>%
  dplyr::mutate(InterProID_norm = stringr::str_extract(InterProID_expanded, "IPR\\d+")) %>%
  dplyr::filter(!is.na(InterProID_norm)) %>%
  dplyr::left_join(map_df, by = c("InterProID_norm" = "InterProID")) %>%
  dplyr::group_by(.rowid) %>%
  dplyr::summarize(
    GO_from_InterPro_IDs   = paste(unique(stats::na.omit(GO_ID)), collapse = ";"),
    .groups = "drop"
  )

df <- df %>%
  dplyr::mutate(.rowid = dplyr::row_number()) %>%
  dplyr::left_join(df_go, by = ".rowid") %>%
  dplyr::mutate(
    GO_from_InterPro_IDs = tidyr::replace_na(GO_from_InterPro_IDs, "")
  ) %>%
  dplyr::select(-.rowid)
# write_tsv(df, "ggsearch-interpro-GO_annotbl.txt", na = "")
```
## Quantification
```bash
# salmon (Since amino acid sequences are not accepted by Salmon, CDS were used.)
salmon index -p 64 --keepDuplicates -t transcripts.fasta.transdecoder_dir/longest_orfs.cds -i SPAdes_index

for i in $(seq 21185407 21185463); do
    salmon quant -i SPAdes_index -l A \
                 -1 "SRR${i}"_1.fastq \
                 -2 "SRR${i}"_2.fastq \
                 -p 64 \
                 -o salmon-"SRR${i}"
done
```
## Condition Specific annotation
```R
library(tximport)
library(dplyr)
library(tidyr)
library(readr)
library(tximport)

s2c <- read.table("sample2condition.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
s2c$group <- gsub(" ", "_", s2c$group)
files <- s2c$path
names(files) <- s2c$sample
txi <- tximport(files, type = "salmon", txOut = TRUE)
txi$length[txi$length == 0] <- 1
tpm <- as.data.frame(txi$abundance)
tpm <- tibble::rownames_to_column(tpm, var = "Gene")
meta <- read_csv("meta_data.csv") 
tpm_long <- pivot_longer(tpm, -Gene, names_to = "Run", values_to = "TPM")
tpm_annotated <- left_join(tpm_long, meta, by = "Run")

selected_genes <- tpm_annotated %>%
  group_by(Gene, LargeGroup) %>%
  summarise(
    meanTPM = mean(TPM, na.rm = TRUE),
    sdTPM   = sd(TPM, na.rm = TRUE),
    CV      = ifelse(meanTPM > 0, sdTPM / meanTPM, NA),
    .groups = "drop"
  ) %>%
  group_by(Gene) %>%
  summarise(
    high_expr_groups = sum(meanTPM >= 1),
    expressed_group = paste(LargeGroup[meanTPM >= 1], collapse = ";"),
    all_CV_below_threshold = all(CV < 1, na.rm = TRUE),
    meanTPM = max(meanTPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    high_expr_groups == 1,
    all_CV_below_threshold == TRUE,
    meanTPM >= 2
  ) %>%
  dplyr::select(Gene, expressed_group, meanTPM)

go_annotated_df <- read.table("go_annotated_results.tsv",
                              sep = "\t", header = TRUE,
                              stringsAsFactors = FALSE, fill = TRUE)

colnames(go_annotated_df)[1] <- "Gene"
final_df <- left_join(go_annotated_df, selected_genes, by = "Gene")
# write.table(final_df, file = "ggsearch-interpro-GO-CS_annotbl.txt", na = "", sep = "\t", quote = FALSE, row.names = FALSE)
```
## DEGs
```R
library(DESeq2)
library(tidyverse)
library(ggplot2)

sampleTable <- data.frame(condition=s2c$group)
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
# dds_5 <- dds[rowSums(counts(dds) >= 5) >= 5, ] # Optional:Low-expression genes will also be filtered in the next step
dds_wt <- DESeq(dds)
# dds_wt_5 <- DESeq(dds_5) # for PCA

comparisons <- list(
  list(c("condition", "primordia", "mycelia"), "deseq2-primordia-vs-mycelia.txt"),
  list(c("condition", "fruiting_body", "primordia"), "deseq2-fruiting_body-vs-primordia.txt"),
  list(c("condition", "mycelia", "fruiting_body"), "deseq2-mycelia-vs-fruiting_body.txt")
)

for (comp in comparisons) {
  contrast <- comp[[1]]         
  outfile <- comp[[2]]         
  res <- results(dds_wt, contrast = contrast)
  res_naomit <- na.omit(res)
  res_sorted <- res_naomit[order(res_naomit$padj), ]
  write.table(res_sorted, outfile, sep = "\t", quote = FALSE)
}

#deseq2-primordia_body-vs-mycelia-e-10 1,926 transcripts
#deseq2-fruiting_body-vs-primordia-e-10 1,739 transcripts
#deseq2-mycelia-vs-fruiting_body-e-10 3,801 transcripts

annot_expr <- read.table("annotbl.txt",
                         sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)

de_files <- list(
  primordia_vs_mycelia     = "deseq2-primordia-vs-mycelia.txt",
  fruiting_vs_primordia    = "deseq2-fruiting_body-vs-primordia.txt",
  mycelia_vs_fruiting      = "deseq2-mycelia-vs-fruiting_body.txt"
)

for (name in names(de_files)) {
  res <- read.table(de_files[[name]], sep="\t", header=TRUE, stringsAsFactors=FALSE)
  res$Gene <- rownames(res)
  df_sub <- res[, c("Gene", "log2FoldChange", "padj")]
  colnames(df_sub) <- c("Gene",
                        paste0("log2FC_", name),
                        paste0("padj_", name))
  annot_expr <- left_join(annot_expr, df_sub, by="Gene")
}

write.table(annot_expr,
            file="annotbl_with_DEGs.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
```
## Benn Diagram
```
library(ggVennDiagram)
library(grid)      

alpha_thr <- 1e-10
set_primordia_vs_mycelia <- annot_expr %>%
  dplyr::filter(!is.na(padj_primordia_vs_mycelia),
                padj_primordia_vs_mycelia <= alpha_thr) %>%
  dplyr::pull(Gene) %>% unique()

set_fruiting_vs_primordia <- annot_expr %>%
  dplyr::filter(!is.na(padj_fruiting_vs_primordia),
                padj_fruiting_vs_primordia <= alpha_thr) %>%
  dplyr::pull(Gene) %>% unique()

set_mycelia_vs_fruiting <- annot_expr %>%
  dplyr::filter(!is.na(padj_mycelia_vs_fruiting),
                padj_mycelia_vs_fruiting <= alpha_thr) %>%
  dplyr::pull(Gene) %>% unique()

x <- list(
  "Primordia vs Mycelia"   = set_primordia_vs_mycelia,
  "Fruiting body vs Primordia"  = set_fruiting_vs_primordia,
  "Mycelia vs Fruiting body"    = set_mycelia_vs_fruiting
)

venn_plot <- ggVennDiagram(
  x[1:3], label_alpha = 0
) +
  scale_fill_gradient(low = "white", high = "turquoise") +
  theme_void() +
  theme(
    plot.margin = unit(c(2, 2, 2, 2), "cm"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  ) +
  coord_cartesian(clip = "off") +                                  
  scale_x_continuous(expand = expansion(mult = 0.08)) +             
  scale_y_continuous(expand = expansion(mult = 0.08))

ggsave("venn_padj1e-10.png", venn_plot, width = 7, height = 6.5, dpi = 300)

annot_df <- read.delim("annotbl_with_DEGs.txt", stringsAsFactors = FALSE)
if (!"Gene" %in% colnames(annot_df) && !is.null(rownames(annot_df))) {
  annot_df$Gene <- rownames(annot_df)
}

in_A <- annot_df$Gene %in% x$primordia_vs_mycelia
in_B <- annot_df$Gene %in% x$fruiting_body_vs_primordia
in_C <- annot_df$Gene %in% x$mycelia_vs_fruiting_body

annot_df$only_A          <-  in_A & !in_B & !in_C
annot_df$only_B          <- !in_A &  in_B & !in_C
annot_df$only_C          <- !in_A & !in_B &  in_C
annot_df$A_and_B         <-  in_A &  in_B & !in_C
annot_df$A_and_C         <-  in_A & !in_B &  in_C
annot_df$B_and_C         <- !in_A &  in_B &  in_C
annot_df$A_and_B_and_C   <-  in_A &  in_B &  in_C

annot_df$in_A <- in_A
annot_df$in_B <- in_B
annot_df$in_C <- in_C

write.table(annot_df,
            file = "annotbl_with_DEGs.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
```
## PCA plot
```
vsd <- vst(dds_wt, blind = FALSE)
vsd_5 <- vst(dds_wt_5, blind = FALSE)

# Use this when utilizing the built-in functions of DESeq2
png("DESeq2_PCA_plot.png", width = 1200, height = 1000, res = 150)
plotPCA(vsd, intgroup = "condition")
dev.off()

# The following is the customized version used in this analysis.
vsd_mat <- assay(vsd_5)              
vsd_df <- t(vsd_mat) %>% as.data.frame()
vsd_df <- vsd_df %>% mutate(Sample = rownames(vsd_df))
pca_input <- vsd_df %>% left_join(meta, by = c("Sample" = "Run"))
pca_result <- prcomp(
  pca_input %>% dplyr::select(-Sample, -Group, -LargeGroup),
  center = TRUE, scale. = TRUE
)
pca_scores <- as.data.frame(pca_result$x) %>%
  mutate(Sample = pca_input$Sample,
         Group = pca_input$Group,
         LargeGroup = pca_input$LargeGroup)
explained_var <- pca_result$sdev^2
explained_var_percent <- explained_var / sum(explained_var) * 100

large_group_colors <- c(
  "mycelia" = "#009E73",      # green
  "primordia" = "#0072B2",    # blue
  "fruiting_body" = "#D55E00" # orange
)

group_shapes <- c(
  "1M_sawdust_media" = 16, "2M_sawdust_media" = 17, "2Md_1Ml_sawdust_media" = 15,
  "2Md_1Ml_sawdust_media_1dc" = 18, "2Md_1Ml_sawdust_media_2dc" = 19,
  "3M_sawdust_media" = 8, "mycelia_on_agar" = 7,
  "mycelia_with_primordia_stage0" = 6, "mycelia_with_primordia_stage1" = 4,
  "primordia_stage0" = 3, "primordia_stage1" = 2, "young_fruiting_body_1" = 0,
  "young_fruiting_body_2" = 1, "cap_1" = 9, "cap_2" = 10, 
  "stipe_1" = 5, "stipe_2" = 14, "gil_D0" = 11,
  "gil_D2" = 12, "gill_D4" = 13
)

pca_scores$Group <- factor(pca_scores$Group, levels = names(group_shapes))
pca_scores$LargeGroup <- factor(pca_scores$LargeGroup, levels = names(large_group_colors))
p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = LargeGroup, shape = Group)) +
  geom_point(size = 4.2, stroke = 0.8) +                     
  theme_minimal(base_size = 12) +                             
  labs(
    x = sprintf("PC1 (%.2f%%)", explained_var_percent[1]),
    y = sprintf("PC2 (%.2f%%)", explained_var_percent[2])
  ) +
  scale_color_manual(values = large_group_colors) +
  scale_shape_manual(values = group_shapes) +
  theme(
    axis.title  = element_text(size = 14),
    axis.text   = element_text(size = 12),
    legend.title= element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.7, "cm")                        
  ) +
  guides(color = guide_legend(override.aes = list(size = 4.2)),
         shape = guide_legend(override.aes = list(size = 4.2)))

ggsave("PCA_plot.png", plot = p, width = 11, height = 8.5, units = "in", dpi = 300, bg = "white")
```
## TopGO
```
##topGO
library(topGO)
library(tidyverse)
df <- read.table("annotbl_with_DEGs.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
go_columns <- c("human_GO", "mouse_GO", "yeast_GO", "uniprot_GO", "interproscan_GO")
for (go_col in go_columns) {
  message("Running topGO for: ", go_col)
  df_bg <- df %>% filter(!is.na(.data[[go_col]]) & .data[[go_col]] != "")
  gene2go <- df_bg %>%
    dplyr::select(query_id, GO = all_of(go_col)) %>%
    dplyr::mutate(GO_list = strsplit(GO, ";")) %>%
    tidyr::unnest(GO_list)
  gene2GO_map <- split(gene2go$GO_list, gene2go$query_id)
  all_genes <- unique(df_bg$query_id)
  sig_genes <- df %>%
    filter(onlyA == TRUE & query_id %in% all_genes) %>%
    pull(query_id)
  geneList <- factor(as.integer(all_genes %in% sig_genes))
  names(geneList) <- all_genes
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = gene2GO_map,
                nodeSize = 5) # run it without this option as well
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 20)
  print(allRes)
  output_file <- paste0("topGO_Fisher_result_", go_col, "_onlyA.tsv")
  write.table(allRes, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
```
