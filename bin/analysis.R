#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(readr)
  library(tximport); library(DESeq2); library(ggplot2)
  library(topGO); library(stringr); library(viridis)
})

opt_list <- list(
  make_option("--meta", type="character"),
  make_option("--annot", type="character"),
  make_option("--salmon_dir", type="character"),
  make_option("--out", type="character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(file.path(opt$out,"tximport"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out,"de"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out,"pca"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out,"topgo"), showWarnings=FALSE, recursive=TRUE)

meta <- read_delim(opt$meta, delim = ifelse(grepl("\\.tsv$", opt$meta),"\\t",","), show_col_types=FALSE)
# 必須列: Run, Group, LargeGroup（任意）, path ない場合は salmon_dir から推測
if (!"path" %in% colnames(meta)) {
  meta$path <- file.path(opt$salmon_dir, paste0("salmon-", meta$Run), "quant.sf")
}

files <- setNames(meta$path, meta$Run)

txi <- tximport(files, type="salmon", txOut=TRUE)
txi$length[txi$length==0] <- 1
tpm <- as.data.frame(txi$abundance) %>% tibble::rownames_to_column("Gene")
write_tsv(tpm, file.path(opt$out,"tximport/tpm.tsv"))

# DESeq2
sampleTable <- data.frame(condition = meta$Group)
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)
dds <- DESeq(dds)

# 3比較（必要ならメタから自動化可）
cmp <- list(
  c("condition","primordia","mycelia"),
  c("condition","fruiting_body","primordia"),
  c("condition","mycelia","fruiting_body")
)
nms <- c("primordia-vs-mycelia","fruiting_body-vs-primordia","mycelia-vs-fruiting_body")
for (i in seq_along(cmp)) {
  res <- results(dds, contrast = cmp[[i]])
  res <- na.omit(res)[order(res$padj),]
  write.table(res, file.path(opt$out,"de", paste0("deseq2-", nms[[i]], ".txt")),
              sep="\t", quote=FALSE)
}

# PCA（VST）
vsd <- vst(dds, blind=FALSE)
vsd_mat <- assay(vsd)
vsd_df  <- as.data.frame(t(vsd_mat))
vsd_df$Sample <- rownames(vsd_df)
pca_in <- vsd_df %>% left_join(meta, by=c("Sample"="Run"))
pc <- prcomp(pca_in %>% select(-Sample, -Group, any_of("LargeGroup")), center=TRUE, scale.=TRUE)
ev <- pc$sdev^2 / sum(pc$sdev^2) * 100
scores <- as.data.frame(pc$x) %>% mutate(Sample=pca_in$Sample, Group=pca_in$Group, LargeGroup=pca_in$LargeGroup)

p <- ggplot(scores, aes(PC1,PC2,color=LargeGroup,shape=Group)) +
  geom_point(size=3, stroke=0.6) + theme_minimal(base_size=12) +
  labs(x=sprintf("PC1 (%.2f%%)",ev[1]), y=sprintf("PC2 (%.2f%%)",ev[2])) +
  theme(legend.key.size = unit(0.6,"cm"))
ggsave(file.path(opt$out,"pca/PCA_plot.png"), p, width=11, height=8.5, dpi=300, bg="white")

# ===== topGO（ALLマージ） =====
annot <- read_tsv(opt$annot, show_col_types = FALSE)
go_cols <- c("GO_human","GO_mouse","GO_yeast","GO_uniprot","GO_from_InterPro_IDs")
go_cols <- go_cols[go_cols %in% colnames(annot)]

gene2go_long <- annot %>%
  select(Node, all_of(go_cols)) %>%
  pivot_longer(cols = all_of(go_cols), names_to = "src", values_to = "GO") %>%
  filter(!is.na(GO), GO!="") %>%
  mutate(GO_list = strsplit(GO, ";")) %>%
  unnest(GO_list) %>%
  mutate(GO_list = str_trim(GO_list)) %>%
  filter(GO_list!="") %>%
  distinct(Node, GO_list)

all_genes <- unique(gene2go_long$Node)
gene2GO_map <- split(gene2go_long$GO_list, gene2go_long$Node)
gene2GO_map <- lapply(gene2GO_map, unique)

# 有意集合の定義（例：primordia vs mycelia で |log2FC|>1 & padj<0.05）
res_pm <- read.delim(file.path(opt$out,"de","deseq2-primordia-vs-mycelia.txt"), check.names=FALSE)
res_pm$Gene <- rownames(res_pm)
sig_genes <- res_pm %>% filter(!is.na(padj), padj<0.05, abs(log2FoldChange)>1) %>% pull(Gene) %>% unique()
sig_genes <- sig_genes[sig_genes %in% all_genes]

geneList <- factor(as.integer(all_genes %in% sig_genes))
names(geneList) <- all_genes

GOdata <- new("topGOdata", ontology="BP", allGenes=geneList,
              annot=annFUN.gene2GO, gene2GO=gene2GO_map)

res_weight01 <- runTest(GOdata, algorithm="weight01", statistic="fisher")
tab <- GenTable(GOdata, Pvalue=res_weight01, topNodes=200, numChar=1000)
tab <- tab %>%
  mutate(Pvalue_num = suppressWarnings(as.numeric(Pvalue)),
         Annotated = as.numeric(Annotated),
         Significant = as.numeric(Significant),
         GeneRatio = Significant/pmax(Annotated,1)) %>%
  filter(is.finite(Pvalue_num)) %>%
  arrange(Pvalue_num) %>%
  slice_head(n=20) %>%
  mutate(Term = factor(Term, levels = rev(Term)))

p2 <- ggplot(tab, aes(x=GeneRatio, y=Term)) +
  geom_point(aes(size=Significant, color=-log10(Pvalue_num))) +
  scale_color_viridis_c(name=expression(-log[10](p))) +
  scale_size_continuous(name="Significant") +
  labs(x="GeneRatio (Significant / Annotated)", y=NULL,
       title="topGO dotplot (weight01)") +
  theme_minimal(base_size=12) + theme(axis.text.y = element_text(size=7))
ggsave(file.path(opt$out,"topgo/topGO_dotplot_weight01.png"), p2, width=9, height=6, dpi=300, bg="white")
