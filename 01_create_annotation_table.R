#!/usr/bin/env Rscript
################################################################################
# Script 1: Annotation Table Creation and GO Term Assignment
# 
# Purpose: Create comprehensive annotation table from ggsearch and InterProScan 
#          results, then add GO terms from multiple sources
#
# Required Input Files:
#   - longest_orfs.pep: FASTA file with ORF peptide sequences
#   - ggsearch-*.txt: Raw ggsearch output files for each database
#   - interpro_results.tsv: InterProScan output (TSV format)
#   - uniprot_go_annotations.tsv: UniProt GO annotations (2-column: ID, GO)
#
# Output:
#   - annotation_table_with_GO.txt: Complete annotation table with GO terms
#
# Usage:
#   Rscript 01_create_annotation_table.R
#
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(biomaRt)
  library(readr)
  library(tidyr)
})

message("================================================================================")
message("Fungal Annotation Workflow - Script 1: Annotation Table & GO Terms")
message("================================================================================\n")

################################################################################
# CONFIGURATION - Modify these paths as needed
################################################################################

# Input files
orf_peptide_file <- "longest_orfs.pep"  # FASTA file with ORF sequences
ggsearch_raw_files <- list(
  human    = "ggsearch-human.txt",
  mouse    = "ggsearch-mouse.txt",
  yeast    = "ggsearch-yeast.txt",
  uniprot  = "ggsearch-uniprot.txt",
  fungidb  = "ggsearch-fungidb.txt"
)
interpro_file <- "interpro_results.tsv"
uniprot_go_file <- "uniprot_go_annotations.tsv"

# Output file
output_file <- "annotation_table_with_GO.txt"

################################################################################
# HELPER FUNCTION: Parse ggsearch output
################################################################################

parse_ggsearch <- function(input_file) {
  # Read the ggsearch output file
  lines <- readLines(input_file)
  
  # Initialize variables
  results <- list()
  query <- ""
  dbhit <- ""
  dbdesc <- ""
  dbgsymbol <- ""
  gnw_frame <- ""
  gnw_expect <- ""
  gnw_score <- ""
  gnw_ident <- ""
  gnw_sim <- ""
  nohit <- FALSE
  
  for (line in lines) {
    # Extract query ID
    if (grepl("^>>>([^,]+), (\\d+) aa", line)) {
      query <- sub("^>>>([^,]+), .*", "\\1", line)
    }
    
    # Extract database hit
    if (grepl("^>>\\S+\\s", line)) {
      dbhit <- sub("^>>(\\S+)\\s.*", "\\1", line)
      dbdesc <- sub("^>>\\S+\\s+", "", line)
      
      # Extract gene symbol if present
      if (grepl("gene_symbol:(\\w+)", dbdesc)) {
        dbgsymbol <- sub(".*gene_symbol:(\\w+).*", "\\1", dbdesc)
      } else {
        dbgsymbol <- ""
      }
      
      # Extract description if present
      if (grepl("description:(.+)", dbdesc)) {
        dbdesc <- sub(".*description:(.+)", "\\1", dbdesc)
      }
    }
    
    # Extract alignment statistics
    if (grepl("gnw_frame:", line)) {
      gnw_frame <- sub(".*gnw_frame:\\s+(\\S+).*", "\\1", line)
    }
    if (grepl("gnw_expect:", line)) {
      gnw_expect <- sub(".*gnw_expect:\\s+(\\S+).*", "\\1", line)
    }
    if (grepl("gnw_score:", line)) {
      gnw_score <- sub(".*gnw_score:\\s+(\\S+).*", "\\1", line)
    }
    if (grepl("gnw_ident:", line)) {
      gnw_ident <- sub(".*gnw_ident:\\s+(\\S+).*", "\\1", line)
    }
    if (grepl("gnw_sim:", line)) {
      gnw_sim <- sub(".*gnw_sim:\\s+(\\S+).*", "\\1", line)
    }
    
    # Check for no hits
    if (grepl("^!! No sequences with E", line)) {
      nohit <- TRUE
    }
    
    # End of record marker
    if (grepl(">>><<<", line)) {
      if (!nohit && query != "" && dbhit != "") {
        results[[length(results) + 1]] <- data.frame(
          query = query,
          dbhit = dbhit,
          gnw_frame = gnw_frame,
          gnw_expect = gnw_expect,
          gnw_score = gnw_score,
          gnw_ident = gnw_ident,
          gnw_sim = gnw_sim,
          dbgsymbol = dbgsymbol,
          dbdesc = dbdesc,
          stringsAsFactors = FALSE
        )
      }
      # Reset variables
      nohit <- FALSE
      dbhit <- ""
      dbdesc <- ""
      dbgsymbol <- ""
      gnw_frame <- ""
      gnw_expect <- ""
      gnw_score <- ""
      gnw_ident <- ""
      gnw_sim <- ""
    }
  }
  
  # Combine all results into a data frame
  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(data.frame())
  }
}

################################################################################
# PART 1: Extract ORF IDs from FASTA file
################################################################################

message("Step 1: Extracting ORF IDs from FASTA file...")

# Check if peptide file exists
if (!file.exists(orf_peptide_file)) {
  stop("Error: ", orf_peptide_file, " not found!")
}

# Read FASTA file and extract IDs from headers
fasta_lines <- readLines(orf_peptide_file)
header_lines <- fasta_lines[grepl("^>", fasta_lines)]

# Extract IDs (everything after '>' until first whitespace)
orf_ids <- data.frame(
  Node = sub("^>(\\S+).*", "\\1", header_lines),
  stringsAsFactors = FALSE
)

message("  Extracted ", nrow(orf_ids), " ORF IDs from ", orf_peptide_file)

################################################################################
# PART 1.5: Parse ggsearch Results
################################################################################

message("\nStep 1.5: Parsing ggsearch output files...")

ggsearch_parsed <- list()
for (db_name in names(ggsearch_raw_files)) {
  raw_file <- ggsearch_raw_files[[db_name]]
  
  if (!file.exists(raw_file)) {
    warning("  Warning: ", raw_file, " not found. Skipping.")
    next
  }
  
  message("  Parsing ", db_name, " (", raw_file, ")...")
  parsed_data <- parse_ggsearch(raw_file)
  
  if (nrow(parsed_data) > 0) {
    # Select relevant columns and rename
    # Format: query | dbhit | gnw_frame | gnw_expect | gnw_score | gnw_ident | gnw_sim | dbgsymbol | dbdesc
    # We'll use: query (Node) | dbhit (ID) | dbgsymbol (Symbol) | dbdesc (Description)
    parsed_data <- parsed_data %>%
      select(query, dbhit, dbgsymbol, dbdesc) %>%
      rename(Node = query)
    
    # Save parsed file for reference (optional)
    parsed_file <- sub("\\.txt$", "-parsing.txt", raw_file)
    write.table(parsed_data, parsed_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message("    Parsed ", nrow(parsed_data), " entries, saved to ", parsed_file)
    
    ggsearch_parsed[[db_name]] <- parsed_data
  } else {
    warning("    No results parsed from ", raw_file)
  }
}

message("  Completed parsing ", length(ggsearch_parsed), " ggsearch files\n")

################################################################################
# PART 2: Create Base Annotation Table
################################################################################

message("\nStep 2: Creating base annotation table from ggsearch results...")

# Initialize annotation table
df <- orf_ids

# Function to add ggsearch results from parsed data
add_ggsearch_results <- function(df, parsed_data, prefix) {
  if (is.null(parsed_data) || nrow(parsed_data) == 0) {
    return(df)
  }
  
  # Rename columns with prefix
  colnames(parsed_data)[-1] <- paste0(prefix, c("", ".1", ".2"))
  
  df <- left_join(df, parsed_data, by = "Node")
  return(df)
}

# Add all ggsearch results
for (db_name in names(ggsearch_parsed)) {
  message("  Adding ", db_name, " annotations...")
  df <- add_ggsearch_results(df, ggsearch_parsed[[db_name]], db_name)
}

# Process UniProt IDs (split sp|P50090|KEL2_YEAST format)
if ("uniprot" %in% colnames(df)) {
  message("  Processing UniProt IDs...")
  valid_uniprot_rows <- grepl("^sp\\|[^|]+\\|[^|]+$", df$uniprot)
  if (any(valid_uniprot_rows, na.rm = TRUE)) {
    split_matrix <- do.call(rbind, strsplit(df$uniprot[valid_uniprot_rows], "\\|"))
    df$uniprot.id <- NA
    df$uniprot.name <- NA
    df$uniprot.id[valid_uniprot_rows] <- split_matrix[, 2]
    df$uniprot.name[valid_uniprot_rows] <- split_matrix[, 3]
    df$uniprot <- df$uniprot.id
    df$uniprot.1 <- df$uniprot.name
    df$uniprot.id <- NULL
    df$uniprot.name <- NULL
  }
}

# Clean up FungiDB columns if needed
if ("fungidb.2" %in% colnames(df) && "fungidb.1" %in% colnames(df)) {
  df$fungidb.1 <- df$fungidb.2
  df$fungidb.2 <- NULL
}

message("  Base annotation table created with ", nrow(df), " entries\n")

################################################################################
# PART 3: Add InterProScan Results
################################################################################

message("\nStep 3: Adding InterProScan annotations...")

if (file.exists(interpro_file)) {
  interpro_raw <- read.table(
    interpro_file, header = FALSE, sep = "\t", 
    stringsAsFactors = FALSE, fill = TRUE, quote = ""
  )
  
  colnames(interpro_raw) <- c(
    "Node", "ID", "Length", "Database", "Accession", "Description",
    "Start", "End", "Evalue", "Status", "Date", "InterProID", "InterProDescription",
    "Extra1", "Extra2"
  )
  
  # Select top hit per Node (lowest E-value)
  interpro_top1 <- interpro_raw %>%
    filter(str_detect(Evalue, "^\\d+\\.?\\d*([eE][-+]?\\d+)?$")) %>%
    mutate(Evalue = as.numeric(Evalue)) %>%
    group_by(Node) %>%
    slice_min(order_by = Evalue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(Node, Database, Accession, Description, InterProID, InterProDescription)
  
  df <- left_join(df, interpro_top1, by = "Node")
  message("  InterProScan annotations added (", nrow(interpro_top1), " entries)\n")
} else {
  warning("Warning: ", interpro_file, " not found. Skipping InterProScan.\n")
}

################################################################################
# PART 4: Add GO Terms from Multiple Sources
################################################################################

message("Step 4: Adding GO terms from multiple sources...\n")

# Function to add GO terms from Ensembl (for peptide IDs)
add_go_terms_ensembl_pep <- function(df, id_col, mart_dataset, species_colname) {
  message("  Retrieving GO terms for ", species_colname, " (Ensembl peptide)...")
  
  if (!id_col %in% colnames(df)) {
    warning("    Column '", id_col, "' not found. Skipping.")
    return(df)
  }
  
  # Clean Ensembl peptide IDs (remove version numbers)
  df[[paste0("clean_", id_col)]] <- gsub("\\..*", "", df[[id_col]])
  ids <- unique(na.omit(df[[paste0("clean_", id_col)]]))
  
  if (length(ids) == 0) {
    warning("    No valid IDs found. Skipping.")
    return(df)
  }
  
  tryCatch({
    mart <- useMart("ensembl", dataset = mart_dataset)
    
    go_raw <- getBM(
      attributes = c("ensembl_peptide_id", "go_id"),
      filters = "ensembl_peptide_id",
      values = ids,
      mart = mart
    )
    
    go_summarized <- go_raw %>%
      group_by(ensembl_peptide_id) %>%
      summarise(!!paste0("GO_", species_colname) := paste(unique(go_id), collapse = ";")) %>%
      rename(!!paste0("clean_", id_col) := ensembl_peptide_id)
    
    df <- left_join(df, go_summarized, by = paste0("clean_", id_col))
    message("    Added GO terms for ", nrow(go_summarized), " genes")
    
  }, error = function(e) {
    warning("    Error retrieving GO terms: ", e$message)
  })
  
  return(df)
}

# Function to add GO terms from Ensembl (for gene IDs)
add_go_terms_ensembl_gene <- function(df, id_col, mart_dataset, species_colname) {
  message("  Retrieving GO terms for ", species_colname, " (Ensembl gene)...")
  
  if (!id_col %in% colnames(df)) {
    warning("    Column '", id_col, "' not found. Skipping.")
    return(df)
  }
  
  ids <- unique(na.omit(df[[id_col]]))
  
  if (length(ids) == 0) {
    warning("    No valid IDs found. Skipping.")
    return(df)
  }
  
  tryCatch({
    mart <- useMart("ensembl", dataset = mart_dataset)
    
    go_raw <- getBM(
      attributes = c("ensembl_gene_id", "go_id"),
      filters = "ensembl_gene_id",
      values = ids,
      mart = mart
    )
    
    go_summarized <- go_raw %>%
      group_by(ensembl_gene_id) %>%
      summarise(!!paste0("GO_", species_colname) := paste(unique(go_id), collapse = ";")) %>%
      rename(!!id_col := ensembl_gene_id)
    
    df <- left_join(df, go_summarized, by = id_col)
    message("    Added GO terms for ", nrow(go_summarized), " genes")
    
  }, error = function(e) {
    warning("    Error retrieving GO terms: ", e$message)
  })
  
  return(df)
}

# Function to add GO terms from UniProt
add_go_terms_uniprot <- function(df, id_col, uniprot_go_file, species_colname = "uniprot") {
  message("  Adding GO terms from UniProt file...")
  
  if (!file.exists(uniprot_go_file)) {
    warning("    File '", uniprot_go_file, "' not found. Skipping.")
    return(df)
  }
  
  if (!id_col %in% colnames(df)) {
    warning("    Column '", id_col, "' not found. Skipping.")
    return(df)
  }
  
  go_df <- read.delim(uniprot_go_file, header = FALSE, 
                      col.names = c("uniprot_id", "go_id"))
  
  go_summarized <- go_df %>%
    group_by(uniprot_id) %>%
    summarise(!!paste0("GO_", species_colname) := paste(unique(go_id), collapse = ";"))
  
  df <- left_join(df, go_summarized, by = setNames("uniprot_id", id_col))
  message("    Added UniProt GO terms for ", nrow(go_summarized), " entries")
  
  return(df)
}

# Add GO terms from each source
df <- add_go_terms_ensembl_pep(df, id_col = "human", 
                               mart_dataset = "hsapiens_gene_ensembl", 
                               species_colname = "human")

df <- add_go_terms_ensembl_pep(df, id_col = "mouse", 
                               mart_dataset = "mmusculus_gene_ensembl", 
                               species_colname = "mouse")

df <- add_go_terms_ensembl_gene(df, id_col = "yeast", 
                                mart_dataset = "scerevisiae_gene_ensembl", 
                                species_colname = "yeast")

df <- add_go_terms_uniprot(df, id_col = "uniprot", 
                           uniprot_go_file = uniprot_go_file)

# Clean up temporary columns
df <- select(df, -starts_with("clean_"))

################################################################################
# PART 5: Add GO Terms from InterPro
################################################################################

message("\nStep 5: Adding GO terms from InterPro IDs...")

if ("InterProID" %in% colnames(df)) {
  tryCatch({
    # Download interpro2go mapping
    interpro2go_url <- "https://current.geneontology.org/ontology/external2go/interpro2go"
    raw <- read_lines(interpro2go_url)
    lines <- raw[!startsWith(raw, "!")]
    
    # Parse InterPro to GO mapping
    pairs_list <- lapply(lines, function(x) {
      ipr <- str_extract_all(x, "InterPro:(IPR\\d+)")[[1]]
      go  <- str_extract_all(x, "(GO:\\d+)")[[1]]
      if (length(ipr) == 0 || length(go) == 0) return(NULL)
      expand.grid(
        IPR = unique(sub("^InterPro:", "", ipr)),
        GO  = unique(go),
        stringsAsFactors = FALSE
      )
    })
    
    map_df <- bind_rows(pairs_list) %>%
      distinct() %>%
      rename(InterProID = IPR, GO_ID = GO)
    
    # Function to split InterPro IDs
    split_interpro <- function(x) {
      x %>%
        str_replace_all("InterPro:", "") %>%
        str_split("\\s*[,;|]\\s*|\\s+")
    }
    
    # Map InterPro IDs to GO terms
    df_go <- df %>%
      mutate(.rowid = row_number(),
             .IPR_list = split_interpro(InterProID)) %>%
      unnest_longer(.IPR_list, values_to = "InterProID_expanded") %>%
      mutate(
        InterProID_expanded = if_else(is.na(InterProID_expanded), "", InterProID_expanded),
        InterProID_expanded = str_trim(InterProID_expanded)
      ) %>%
      filter(InterProID_expanded != "", InterProID_expanded != "-") %>%
      mutate(InterProID_norm = str_extract(InterProID_expanded, "IPR\\d+")) %>%
      filter(!is.na(InterProID_norm)) %>%
      left_join(map_df, by = c("InterProID_norm" = "InterProID")) %>%
      group_by(.rowid) %>%
      summarize(
        GO_from_InterPro_IDs = paste(unique(na.omit(GO_ID)), collapse = ";"),
        .groups = "drop"
      )
    
    df <- df %>%
      mutate(.rowid = row_number()) %>%
      left_join(df_go, by = ".rowid") %>%
      mutate(GO_from_InterPro_IDs = replace_na(GO_from_InterPro_IDs, "")) %>%
      select(-.rowid)
    
    message("  InterPro GO terms added successfully")
    
  }, error = function(e) {
    warning("  Error adding InterPro GO terms: ", e$message)
  })
} else {
  message("  No InterProID column found. Skipping.")
}

################################################################################
# PART 6: Save Output
################################################################################

message("\nStep 6: Saving annotation table...")

write_tsv(df, output_file, na = "")

message("  Output saved to: ", output_file)
message("  Total entries: ", nrow(df))
message("  Total columns: ", ncol(df))

# Summary statistics
go_columns <- grep("^GO_", colnames(df), value = TRUE)
if (length(go_columns) > 0) {
  message("\nGO Term Statistics:")
  for (col in go_columns) {
    n_with_go <- sum(df[[col]] != "" & !is.na(df[[col]]))
    pct <- round(100 * n_with_go / nrow(df), 2)
    message("  ", col, ": ", n_with_go, " (", pct, "%)")
  }
}

message("\n================================================================================")
message("Annotation table creation completed successfully!")
message("================================================================================\n")
