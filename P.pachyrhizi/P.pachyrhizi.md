## Computational Environment

All analyses were conducted on a **MacBook Pro (Apple M2 Max)** running **macOS 15.5**.  
A **Docker-based Ubuntu environment** was used unless otherwise specified for individual steps.

## Iso-seq sequences
- GenBank accessions: `GHWK00000000`, `GHWL00000000`

```bash
cat GHWK01.1.fsa_nt GHWL01.1.fsa_nt > P.pachyrhizi-IsoSeq.fasta
seqkit rmdup -s P.pachyrhizi-IsoSeq.fasta -o P.pachyrhizi-IsoSeq-remdup.fasta
[INFO] 74 duplicated records removed

# BUSCO assessment
#busco
mamba install -c conda-forge -c bioconda busco=5.8.0
busco -i P.pachyrhizi-IsoSeq-remdup.fasta -o busco-results -l eukaryota_odb10 -m transcriptome -c 8
    ----------------------------------------------------
    |Results from dataset eukaryota_odb10               |
    ----------------------------------------------------
    |C:33.3%[S:22.4%,D:11.0%],F:11.8%,M:54.9%,n:255     |
    |85    Complete BUSCOs (C)                          |
    |57    Complete and single-copy BUSCOs (S)          |
    |28    Complete and duplicated BUSCOs (D)           |
    |30    Fragmented BUSCOs (F)                        |
    |140    Missing BUSCOs (M)                          |
    |255    Total BUSCO groups searched                 |
    ----------------------------------------------------
```

## Prepare fastq files

```bash
# Download 57 RNA-seq samples of Lentinula edodes (H600) from SRA
prefetch --option-file SRR.txt
for i in $(seq 10130097 10130116); do
  fasterq-dump "SRR${i}"
done

# Trimming reads
docker attach trimgalore-0.6.10
for i in $(seq 10130097 10130116); do
  trim_galore --paired "SRR${i}_1.fastq" "SRR${i}_2.fastq"
done
```
## ORF Prediction with TransDecoder

```bash
# Run TransDecoder (executed on host system due to Docker error)
#TransDecoder
sudo TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t P.pachyrhizi-IsoSeq-remdup.fasta
sudo TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t P.pachyrhizi-IsoSeq-remdup.fasta

# Check number of predicted ORFs (227,580)
grep -c "^>" transcripts.fasta.transdecoder_dir/longest_orfs.pep
```
## ggsearch & InterProScan

```bash
# ggsearch
time ggsearch36 -Q -T 8 -d 1 -m10 -E 0.1 P.pachyrhizi-IsoSeq-remdup.fasta.transdecoder_dir/longest_orfs.pep mac/Homo_sapiens.GRCh38.pep.all.fa > ggsearch-P.pachyrhizi-human.txt; ggsearch36 -Q -T 8 -d 1 -m10 -E 0.1 P.pachyrhizi-IsoSeq-remdup.fasta.transdecoder_dir/longest_orfs.pep mac/Mus_musculus.GRCm39.pep.all.fa > ggsearch-P.pachyrhizi-mouse.txt; ggsearch36 -Q -T 8 -d 1 -m10 -E 0.1 P.pachyrhizi-IsoSeq-remdup.fasta.transdecoder_dir/longest_orfs.pep mac/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa > ggsearch-P.pachyrhizi-S.cerevisiae.txt; ggsearch36 -Q -T 8 -d 1 -m10 -E 0.1 P.pachyrhizi-IsoSeq-remdup.fasta.transdecoder_dir/longest_orfs.pep mac/uniprot_sprot.fasta > ggsearch-P.pachyrhizi-uniprot.txt; ggsearch36 -Q -T 128 -d 1 -m10 -E 0.1 P.pachyrhizi-IsoSeq-remdup.fasta.transdecoder_dir/longest_orfs.pep fungidb-68.fasta > ggsearch-P.pachyrhizi-fungidb.txt &

# InterProScan
# Remove asterisks (*) at the end of amino acid sequences from the FASTA file
sed 's/\*//g' ../P.pachyrhizi/P.pachyrhizi-IsoSeq-remdup.fasta.transdecoder_dir/longest_orfs.pep > longest_orfs_cleaned.pep
# Run InterProScan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/interproscan-5.67-99.0-64-bit.tar.gz
tar -pxvzf interproscan-5.67-99.0-*-bit.tar.gz
python3 setup.py -f interproscan.properties
./interproscan.sh -i longest_orfs_cleaned.pep -f tsv --goterms
```
## annotation table
```
# parse
perl -nle 'print $1 if(/^>(\S+)/)' P.pachyrhizi-IsoSeq-remdup.fasta.transdecoder_dir/longest_orfs.pep > longest_orf-pid.txt
for file in ggsearch-*.txt; do cat "$file" | perl 15parseggsearch.pl > "${file%.txt}-parsing.txt"; done

# make annotation table
cat longest_orf-pid.txt | perl 15mkannotbl.pl ggsearch-P.pachyrhizi-human-parsing.txt | perl 15mkannotbl.pl ggsearch-P.pachyrhizi-mouse-parsing.txt | perl 15mkannotbl.pl ggsearch-P.pachyrhizi-S.cerevisiae-parsing.txt | perl 15mkannotbl.pl ggsearch-P.pachyrhizi-uniprot-parsing.txt | perl 15mkannotbl.pl ggsearch-P.pachyrhizi-fungidb-parsing.txt > ggsearch_annotbl.txt

# add InterProScan
R
library(dplyr)
library(stringr)

interpro_raw <- read.table(
  "longest_orfs_cleaned.pep.tsv",
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
# write.table(df, "annotbl.txt", sep="\t", row.names = FALSE, na = "")
```
## Quantification
```bash
# salmon (Since amino acid sequences are not accepted by Salmon, CDS were used.)
salmon index -p 64 --keepDuplicates -t P.pachyrhizi-IsoSeq-remdup.fasta.transdecoder_dir/longest_orfs.cds -i salmon_index

for i in $(seq 10130097 10130116); do
    salmon quant -i salmon_index -l A \
                 -1 "SRR${i}"_val_1.fastq \
                 -2 "SRR${i}"_val_2.fastq \
                 -p 64 \
                 -o salmon-"SRR${i}"
done
```
## Add TPM
```R
library(dplyr)
library(purrr)
library(tidyr)
sample_numbers <- 10130097:10130116
sample_names <- sprintf("SRR%d", sample_numbers)
file_paths <- sprintf("salmon-%s.sf", sample_names)
samples <- data.frame(
  sample = sample_names,
  path = file_paths
)

sample_data_list <- samples %>%
  mutate(
    data = map2(
      path, sample,
      ~ read.csv(.x, sep = "\t", check.names = FALSE) %>%
          dplyr::select(Name, TPM) %>%         
          setNames(c("Name", .y))
    )
  ) %>%
  pull(data)

combined_data <- reduce(sample_data_list, full_join, by = "Name")

tpm <- combined_data %>%
  rename_with(~ paste0("TPM_", .x), -Name) %>%
  distinct(Name, .keep_all = TRUE)

df <- df %>% left_join(tpm, by = c("Node" = "Name"))
# write.table(df, "annotbl.txt", sep="\t", na = "")
```
## DEGs
```R
#https://www.bioconductor.org/packages/release/bioc/vignettes/maSigPro/inst/doc/maSigProUsersGuide.pdf
BiocManager::install("maSigPro")
library(maSigPro)
count <- df %>% dplyr::select(Node, dplyr::matches("^TPM_")) %>% 
  dplyr::rename_with(~ sub("^TPM_", "", .x)) %>%
  column_to_rownames("Node")
edesign <- read.table("edesign.txt", header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
colnames(count)
rownames(edesign)
colnames(edesign)
rownames(count)
design <- make.design.matrix(edesign, degree = 3)
design$groups.vector
fit <- p.vector(count, design, Q = 0.05, MT.adjust = "BH", min.obs = 10)
fit$i # returns the number of significant genes
fit$alfa # gives p-value at the Q false discovery control level
fit$SELEC # is a matrix with the significant genes and their expression values
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "all")
names(sigs)
names(sigs$sig.genes)
sigs$sig.genes$g # 3078
see.genes(sigs$sig.genes, show.fit = T, dis =design$dis,
 cluster.method="hclust" ,cluster.data = 1, k = 4) 
clusters <- see.genes(sigs$sig.genes, show.fit = TRUE, dis = design$dis,
                      cluster.method = "hclust", cluster.data = 1, k = 4)
clustered_genes <- clusters$cut
pvalues_df <- sigs$sig.genes$sig.pvalues
for (i in unique(clustered_genes)) {
  genes_in_cluster <- names(clustered_genes[clustered_genes == i])  
  cluster_pvalues <- pvalues_df[rownames(pvalues_df) %in% genes_in_cluster, ]  
  
  write.csv(cluster_pvalues, file = paste0("Cluster_", i, "_genes_with_pvalues.csv"), row.names = TRUE, quote = FALSE)
}

print(table(clustered_genes))

cluster_df <- data.frame(
  Node = names(clustered_genes),
  maSigPro_Cluster = paste0("Cluster", as.integer(clustered_genes)),
  stringsAsFactors = FALSE
)

pval_df <- pvalues_df %>%
  as.data.frame() %>%
  rownames_to_column(var = "Node")

colnames(pval_df)[colnames(pval_df) != "Node"] <- paste0("maSigPro_p_", make.names(colnames(pval_df)[colnames(pval_df) != "Node"]))

add_df <- cluster_df %>%
  full_join(pval_df, by = "Node") %>%
  mutate(maSigPro_Significant = !is.na(maSigPro_Cluster))

final_df_with_DEGs <- df %>% left_join(add_df, by = "Node")
write.table(df, "annotbl-with-DEGs.txt", sep="\t", row.names = FALSE, na = "")
```
## PCA plot
```
library(tidyverse)
library(ggplot2)
library(readr)

combined_data_filtered <- combined_data %>%
  filter(if_any(-Name, ~ .x >= 1.0))
meta_data <- read_csv("meta_data.csv")
tpm_subset <- dplyr::select(combined_data_filtered, dplyr::all_of(meta_data$Run))
tpm_log <- log2(as.matrix(tpm_subset) + 1)
pca_result <- prcomp(t(tpm_log), center = TRUE, scale. = TRUE)

pca_scores <- as.data.frame(pca_result$x) %>%
  mutate(Sample = rownames(.)) %>%
  left_join(meta_data, by = c("Sample" = "Run"))

explained_var <- pca_result$sdev^2
explained_var_percent <- explained_var / sum(explained_var) * 100

colors <- c("3-dai" = "#EC407A", "7-dai" = "#FF7043", "10-dai" = "#26C6DA", "14-dai" = "#66BB6A")
pca_scores$Group <- factor(pca_scores$Group, levels = c("3-dai", "7-dai", "10-dai", "14-dai"))

p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3.8) +                                  
  theme_minimal() +
  labs(
    title = NULL,                                            
    x = sprintf("PC1 (%.2f%%)", explained_var_percent[1]),  
    y = sprintf("PC2 (%.2f%%)", explained_var_percent[2]),
    color = "Group"
  ) +
  scale_color_manual(values = colors) +
  theme(legend.position = "right")

ggsave("PCA_plot.png", plot = p, width = 11, height = 8.5, units = "in", dpi = 300, bg = "white")
```
## TopGO
```
##topGO
library(topGO)
library(tidyverse)
library(Rgraphviz)

# all GO
gene2go_long <- final_df_with_DEGs %>%
  pivot_longer(cols = all_of(go_columns), names_to = "source", values_to = "GO") %>%
  filter(!is.na(GO), GO != "") %>%
  mutate(GO_list = strsplit(GO, ";")) %>%
  unnest(GO_list) %>%
  mutate(GO_list = str_trim(GO_list)) %>%
  filter(GO_list != "") %>%
  distinct(Node, GO_list)

all_genes <- unique(gene2go_long$Node)
gene2GO_map <- split(gene2go_long$GO_list, gene2go_long$Node)
gene2GO_map <- lapply(gene2GO_map, unique)

sig_genes <- final_df_with_DEGs %>%
  dplyr::filter(maSigPro_Cluster == "Cluster4", Node %in% all_genes) %>%
  dplyr::pull(Node) %>%
  unique()

geneList <- factor(as.integer(all_genes %in% sig_genes))
names(geneList) <- all_genes

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO_map)

algos <- c("classic", "elim", "weight01")
results <- list()
for (algo in algos) {
  message(sprintf("-- %s Algorithm --", algo))
  res <- runTest(GOdata, algorithm = algo, statistic = "fisher")
  results[[algo]] <- res  
  res_named <- setNames(list(res), paste0(algo, "Fisher"))
  args <- c(
    list(GOdata),                          
    res_named,
    list(orderBy = paste0(algo, "Fisher"),
         topNodes = 200, numChar = 1000)
  )
  allRes <- do.call(topGO::GenTable, args)
  write.table(allRes,
              file = sprintf("topGO_%s_result.tsv", algo),
              sep = "\t", quote = FALSE, row.names = FALSE)
  scores <- score(res)
  png(sprintf("topGO_graph_%sFisher.png", algo), width = 1600, height = 1200, res = 300)
  showSigOfNodes(GOdata, scores, firstSigNodes = 5, useInfo = "all")
  dev.off()
}
lapply(results, function(x) head(score(x)))


# dot plot
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

res_classic  <- runTest(GOdata, algorithm = "classic",  statistic = "fisher")
res_elim     <- runTest(GOdata, algorithm = "elim",     statistic = "fisher")
res_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

crRes <- GenTable(GOdata, classicFisher = res_classic, topNodes = 200, numChar = 1000)
elRes <- GenTable(GOdata, elimFisher = res_elim, topNodes = 200, numChar = 1000)
weRes <- GenTable(GOdata, weight01Fisher = res_weight01, topNodes = 200, numChar = 1000)

topgo_dotplot <- function(GOdata, result_obj, method_label = "weight01",
                          top_n = 20, label_wrap = 40,
                          width = 8, height = 6, dpi = 300,
                          outfile = NULL) {

  tab <- GenTable(
    GOdata,
    Pvalue = result_obj,
    topNodes = max(200, top_n) ,
    numChar = 1000
  )

  df <- tab %>%
    mutate(
      Pvalue_num  = suppressWarnings(as.numeric(Pvalue)),
      Annotated   = as.numeric(Annotated),
      Significant = as.numeric(Significant),
      GeneRatio   = Significant / pmax(Annotated, 1),      
      Term_wrapped = str_wrap(Term, width = label_wrap)
    ) %>%
    filter(is.finite(Pvalue_num)) %>%
    arrange(Pvalue_num) %>%       
    slice_head(n = top_n) %>%
    mutate(Term_wrapped = factor(Term_wrapped, levels = rev(Term_wrapped)))

  if (nrow(df) == 0) {
    message("[", method_label, "] No significant terms were found")
    return(invisible(NULL))
  }

  p <- ggplot(df, aes(x = GeneRatio, y = Term_wrapped)) +
    geom_point(aes(size = Significant, color = -log10(Pvalue_num))) +
    scale_size_continuous(name = "Significant") +
    scale_color_viridis_c(name = expression(-log[10](p)), option = "D") +
    labs(
      x = "GeneRatio (Significant / Annotated)",
      y = NULL,
      title = paste0("topGO dotplot (", method_label, ")")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 7)
    )

  if (is.null(outfile)) {
    outfile <- paste0("topGO_dotplot_", method_label, ".png")
  }
  ggsave(outfile, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", outfile)

  invisible(list(table = df, plot = p, file = outfile))
}

topgo_dotplot(GOdata, res_classic,  method_label = "classic",  outfile = "topGO_dotplot_classic.png")
topgo_dotplot(GOdata, res_elim,     method_label = "elim",     outfile = "topGO_dotplot_elim.png")
topgo_dotplot(GOdata, res_weight01, method_label = "weight01", outfile = "topGO_dotplot_weight01.png")

```
## Blast2GO
```
diamond blastp -d ../nr-2024 -q P.pachyrhizi-IsoSeq-remdup.fasta.transdecoder_dir/longest_orfs.pep -o P.pachyrhizi_orf-blastp.xml -p 96 -f 5 -k 1 &
# blasted → 7,774, GO annotated → 1,895 (ここまで変更)
## topGO
data <- read.delim("blast2go-topGO.txt", stringsAsFactors = FALSE)
data$GO_IDs_clean <- gsub("(C:|P:|F:)", "", data$GO.IDs)
data_long <- data %>%
   separate_rows(GO_IDs_clean, sep = ";\\s*") %>%
   filter(GO_IDs_clean != "")
geneID2GO <- tapply(data_long$GO_IDs_clean, data_long$SeqName, unique)
sigGenes <- data %>%
  filter(Cluster_4 != "#N/A") %>%
  pull(SeqName)

allGenes <- unique(data$SeqName) 

geneList <- factor(as.integer(allGenes %in% sigGenes))
names(geneList) <- allGenes 

GOdata <- new("topGOdata",
               ontology = "BP", 
               allGenes = geneList,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO) 

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 100, numChar = 1000)
write.table(allRes, file = "topGO_blast2go_weight01.txt", sep = "\t", quote = FALSE, row.names = FALSE)

scores <- score(resultFisher) 
png("topGO_blast2go_weight01Fisher.png", width = 1600, height = 1200, res = 180)
showSigOfNodes(GOdata, scores, firstSigNodes = 5, useInfo = "all")
dev.off()

# dot plot
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

res_classic  <- runTest(GOdata, algorithm = "classic",  statistic = "fisher")
res_elim     <- runTest(GOdata, algorithm = "elim",     statistic = "fisher")
res_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

crRes <- GenTable(GOdata, classicFisher = res_classic, topNodes = 200)
elRes <- GenTable(GOdata, elimFisher = res_elim, topNodes = 200)
weRes <- GenTable(GOdata, weight01Fisher = res_weight01, topNodes = 200)

write.table(crRes, file = "topGO_elim_result_mergedGO_cluster4.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

topgo_dotplot <- function(GOdata, result_obj, method_label = "weight01",
                          top_n = 20, label_wrap = 40,
                          width = 8, height = 6, dpi = 300,
                          outfile = NULL) {

  tab <- GenTable(
    GOdata,
    Pvalue = result_obj,
    topNodes = max(200, top_n) , numChar = 1000
  )

  df <- tab %>%
    mutate(
      Pvalue_num  = suppressWarnings(as.numeric(Pvalue)),
      Annotated   = as.numeric(Annotated),
      Significant = as.numeric(Significant),
      GeneRatio   = Significant / pmax(Annotated, 1),      
      Term_wrapped = str_wrap(Term, width = label_wrap)
    ) %>%
    filter(is.finite(Pvalue_num)) %>%
    arrange(Pvalue_num) %>%       
    slice_head(n = top_n) %>%
    mutate(Term_wrapped = factor(Term_wrapped, levels = rev(Term_wrapped)))

  if (nrow(df) == 0) {
    message("[", method_label, "] No significant terms were found")
    return(invisible(NULL))
  }

  p <- ggplot(df, aes(x = GeneRatio, y = Term_wrapped)) +
    geom_point(aes(size = Significant, color = -log10(Pvalue_num))) +
    scale_size_continuous(name = "Significant") +
    scale_color_viridis_c(name = expression(-log[10](p)), option = "D") +
    labs(
      x = "GeneRatio (Significant / Annotated)",
      y = NULL,
      title = paste0("topGO dotplot (", method_label, ")")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title  = element_text(hjust = 0)
    )

  if (is.null(outfile)) {
    outfile <- paste0("topGO_dotplot_", method_label, ".png")
  }
  ggsave(outfile, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", outfile)

  invisible(list(table = df, plot = p, file = outfile))
}

topgo_dotplot(GOdata, res_classic,  method_label = "classic",  outfile = "topGO_dotplot_classic.png")
topgo_dotplot(GOdata, res_elim,     method_label = "elim",     outfile = "topGO_dotplot_elim.png")
topgo_dotplot(GOdata, res_weight01, method_label = "weight01", outfile = "topGO_dotplot_weight01.png")
```


