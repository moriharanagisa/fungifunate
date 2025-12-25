# fungifunate

A functional annotation workflow specifically designed for fungal transcriptomes.

## Overview

fungifunate integrates multiple databases and tools to provide comprehensive functional annotation for fungal RNA-seq data, enabling:
- Multi-database homology search (human, mouse, yeast, UniProt, FungiDB)
- GO term assignment from multiple sources
- Condition-specific gene annotation
- Differential expression analysis
- Functional enrichment analysis

### Prerequisites

- **R** (≥ 4.0) with packages:
  - `tximport`, `DESeq2`, `dplyr`, `tidyr`, `readr`, `ggplot2`, `topGO`, `stringr`, `biomaRt`
- **Bioinformatics tools**:
  - [SPAdes](https://github.com/ablab/spades) (for transcriptome assembly)
  - [TransDecoder](https://github.com/TransDecoder/TransDecoder) (for ORF prediction)
  - [FASTA36](https://github.com/wrpearson/fasta36) (for ggsearch)
  - [InterProScan](https://www.ebi.ac.uk/interpro/download/) (for domain annotation)
  - [Salmon](https://github.com/COMBINE-lab/salmon) (for quantification)
  - [trim_galore](https://github.com/FelixKrueger/TrimGalore) (optional, for read trimming)

### Reference Databases

Download the following reference databases:
```bash
# Ensembl protein sequences
wget http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
wget http://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz
wget http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz

# UniProt Swiss-Prot
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# FungiDB
# Visit https://fungidb.org/fungidb/ and download the latest protein sequences

# Decompress all files
gunzip *.gz
```

## Workflow

### Step 0: Quality Control and Trimming
```bash
trim_galore --paired reads_1.fastq reads_2.fastq
```

### Step 1: Transcriptome Assembly
```bash
# Assemble transcriptome using SPAdes
rnaspades.py \
  -1 reads_1_val_1.fq \
  -2 reads_2_val_2.fq \
  -o rnaspades_output
```

### Step 2: ORF Prediction
```bash
# Predict open reading frames using TransDecoder
TransDecoder.LongOrfs -t rnaspades_output/transcripts.fasta
TransDecoder.Predict -t rnaspades_output/transcripts.fasta
```

### Step 3: Homology Search and Domain Annotation

#### 3.1 Homology Search (ggsearch)
```bash
# Search against multiple databases
ggsearch36 -Q -T 96 -d 1 -m 10 -E 0.1 \
  transcripts.fasta.transdecoder_dir/longest_orfs.pep \
  Homo_sapiens.GRCh38.pep.all.fa > ggsearch-human.txt

ggsearch36 -Q -T 96 -d 1 -m 10 -E 0.1 \
  transcripts.fasta.transdecoder_dir/longest_orfs.pep \
  Mus_musculus.GRCm39.pep.all.fa > ggsearch-mouse.txt

ggsearch36 -Q -T 96 -d 1 -m 10 -E 0.1 \
  transcripts.fasta.transdecoder_dir/longest_orfs.pep \
  Saccharomyces_cerevisiae.R64-1-1.pep.all.fa > ggsearch-yeast.txt

ggsearch36 -Q -T 96 -d 1 -m 10 -E 0.1 \
  transcripts.fasta.transdecoder_dir/longest_orfs.pep \
  fungidb-68.fasta > ggsearch-fungidb.txt

ggsearch36 -Q -T 96 -d 1 -m 10 -E 0.1 \
  transcripts.fasta.transdecoder_dir/longest_orfs.pep \
  uniprot_sprot.fasta > ggsearch-uniprot.txt
```

#### 3.2 Domain Annotation (InterProScan)
```bash
# Remove stop codons (asterisks)
sed 's/\*//g' transcripts.fasta.transdecoder_dir/longest_orfs.pep > longest_orfs_cleaned.pep

# Run InterProScan
interproscan.sh -i longest_orfs_cleaned.pep -f tsv --goterms -o interpro_results.tsv
```

### Step 4: Expression Quantification
```bash
# Build Salmon index
salmon index -p 64 --keepDuplicates \
  -t transcripts.fasta.transdecoder_dir/longest_orfs.cds \
  -i salmon_index

# Quantify expression for each sample
for sample in sample1 sample2 sample3; do
    salmon quant -i salmon_index -l A \
                 -1 ${sample}_1.fastq \
                 -2 ${sample}_2.fastq \
                 -p 64 \
                 -o salmon-${sample}
done
```

### Step 5: Functional Annotation

Create annotation table integrating all databases and GO terms:
```bash
Rscript 01_create_annotation_table.R
```

**Required input files:**
- `longest_orfs.pep`
- `ggsearch-*.txt`
- `interpro_results.tsv`
- `uniprot_go_annotations.tsv`

**Output:**
- `annotation_table_with_GO.txt`

### Step 6: Differential Expression & Enrichment Analysis

Prepare sample metadata file `sample2condition.txt`:
```
sample	path	group
sample1	salmon-sample1/quant.sf	control
sample2	salmon-sample2/quant.sf	control
sample3	salmon-sample3/quant.sf	treatment
sample4	salmon-sample4/quant.sf	treatment
```

Run analysis:
```bash
Rscript 02_expression_analysis.R
```

**Required input files:**
- `annotation_table_with_GO.txt`
- `sample2condition.txt`
- Salmon `quant.sf`

**Output:**
- `annotation_table_final.txt` 
- `PCA_plot.png`
- `deseq2-*.txt`
- `topGO_*.txt` and `topGO_dotplot_*.png`

## Configuration

Edit the comparison groups in `02_expression_analysis.R`:
```r
# Analysis parameters
tpm_threshold <- 1          # Minimum mean TPM for condition-specific annotation
cv_threshold <- 1           # Maximum CV for condition-specific annotation
min_mean_tpm <- 2          # Minimum mean TPM for final selection
deseq_padj_threshold <- 1e-9  # Adjusted p-value threshold for DEGs

# Comparison groups
comparisons <- list(
  list(c("condition", "group2", "group1"), "deseq2-group2-vs-group1.txt"),
  list(c("condition", "group3", "group2"), "deseq2-group3-vs-group2.txt")
)
```

## Acknowledgments

- [DESeq2](https://bioconductor.org/packages/DESeq2/)
- [topGO](https://bioconductor.org/packages/topGO/)
- [Salmon](https://combine-lab.github.io/salmon/)
- [SPAdes](http://cab.spbu.ru/software/spades/)
- [TransDecoder](https://github.com/TransDecoder/TransDecoder)
