#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(stringr); library(readr); library(tidyr)
})

opt_list <- list(
  make_option("--annot", type="character"),
  make_option("--interpro", type="character"),
  make_option("--uniprot", type="character"),
  make_option("--out", type="character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

gg <- read.delim(opt$annot, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE, quote="")
colnames(gg) <- c("Node","human","human.1","human.2","mouse","mouse.1","mouse.2",
                  "yeast","yeast.1","yeast.2","uniprot","uniprot.1","uniprot.2",
                  "fungidb","fungidb.1","fungidb.2")

# UniProtのID|NAME整形
valid <- grepl("^sp\\|[^|]+\\|[^|]+$", gg$uniprot)
if (any(valid)) {
  sp <- do.call(rbind, strsplit(gg$uniprot[valid], "\\|"))
  gg$uniprot[valid]  <- sp[,2] # ID
  gg$uniprot.1[valid]<- sp[,3] # NAME
}
gg$fungidb.1 <- gg$fungidb.2; gg$fungidb.2 <- NULL

ipr <- read.delim(opt$interpro, header=FALSE, sep="\t", stringsAsFactors=FALSE, fill=TRUE, quote="")
colnames(ipr) <- c("Node","ID","Length","Database","Accession","Description",
                   "Start","End","Evalue","Status","Date","InterProID","InterProDescription",
                   "Extra1","Extra2")
ipr_top1 <- ipr %>%
  filter(str_detect(Evalue, "^\\d+\\.?\\d*([eE][-+]?\\d+)?$")) %>%
  mutate(Evalue = as.numeric(Evalue)) %>%
  group_by(Node) %>% slice_min(order_by=Evalue, n=1, with_ties=FALSE) %>% ungroup() %>%
  select(Node, Database, Accession, Description, InterProID, InterProDescription)

df <- gg %>% left_join(ipr_top1, by="Node")

# UniProt→GO 付与
go_uniprot <- read.delim(opt$uniprot, header=FALSE, col.names=c("uniprot_id","go_id"))
go_u <- go_uniprot %>% group_by(uniprot_id) %>% summarise(GO_uniprot=paste(unique(go_id),collapse=";"))
df <- df %>% left_join(go_u, by=c("uniprot"="uniprot_id"))

# InterPro→GO 付与（外部 interpro2go を直読み）
interpro2go_url <- "https://current.geneontology.org/ontology/external2go/interpro2go"
raw  <- readr::read_lines(interpro2go_url)
lines<- raw[!startsWith(raw, "!")]
pairs_list <- lapply(lines, function(x){
  ipr <- str_extract_all(x, "InterPro:(IPR\\d+)")[[1]]
  go  <- str_extract_all(x, "(GO:\\d+)")[[1]]
  if (length(ipr)==0 || length(go)==0) return(NULL)
  expand.grid(IPR=unique(sub("^InterPro:","",ipr)), GO=unique(go), stringsAsFactors=FALSE)
})
map_df <- bind_rows(pairs_list) %>% distinct() %>% rename(InterProID=IPR, GO_ID=GO)

split_interpro <- function(x){
  x %>% str_replace_all("InterPro:","") %>% str_split("\\s*[,;|]\\s*|\\s+")
}

df_go <- df %>%
  mutate(.rowid=row_number(), .IPR_list=split_interpro(InterProID)) %>%
  unnest_longer(.IPR_list, values_to="InterProID_expanded") %>%
  mutate(InterProID_expanded = str_trim(replace_na(InterProID_expanded,""))) %>%
  filter(InterProID_expanded!="", InterProID_expanded!="-") %>%
  mutate(InterProID_norm = str_extract(InterProID_expanded, "IPR\\d+")) %>%
  filter(!is.na(InterProID_norm)) %>%
  left_join(map_df, by=c("InterProID_norm"="InterProID")) %>%
  group_by(.rowid) %>%
  summarise(GO_from_InterPro_IDs=paste(unique(na.omit(GO_ID)),collapse=";"), .groups="drop")

df <- df %>%
  mutate(.rowid=row_number()) %>%
  left_join(df_go, by=".rowid") %>%
  select(-.rowid) %>%
  mutate(GO_from_InterPro_IDs = replace_na(GO_from_InterPro_IDs,""))

write_tsv(df, opt$out, na="")
