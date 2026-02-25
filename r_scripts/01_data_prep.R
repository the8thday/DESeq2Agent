#!/usr/bin/env Rscript
# 01_data_prep.R - Data preparation: gene ID detection, normalization, DDS construction
# Called by: Python RScriptRunner
# Input:  config.json (path passed as first CLI argument)
# Output: data_prep_summary.json, counts_ensembl.csv, id_mapping.csv, dds.rds

suppressPackageStartupMessages({
  library(jsonlite)
  library(DESeq2)
  library(AnnotationDbi)
  library(dplyr)
  library(tibble)
})

# Null-coalescing operator (must be defined before use)
`%||%` <- function(a, b) if (!is.null(a) && !identical(a, NA)) a else b

# --- Read config ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 01_data_prep.R <config.json>")
config <- fromJSON(args[1], simplifyDataFrame = FALSE)

counts_file     <- config$counts_file
metadata_file   <- config$metadata_file
output_dir      <- config$output_dir
min_count       <- as.integer(config$min_count_threshold %||% 10L)
min_frac        <- as.numeric(config$min_samples_fraction %||% 0.5)
species         <- tolower(config$species %||% "human")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("[01_data_prep] Reading counts: ", counts_file)
counts_raw <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts_raw <- round(as.matrix(counts_raw))
storage.mode(counts_raw) <- "integer"

message("[01_data_prep] Reading metadata: ", metadata_file)
metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)

# Ensure sample order matches
common_samples <- intersect(colnames(counts_raw), rownames(metadata))
if (length(common_samples) == 0) stop("No overlapping samples between counts and metadata")
counts_raw <- counts_raw[, common_samples, drop = FALSE]
metadata   <- metadata[common_samples, , drop = FALSE]

n_samples <- ncol(counts_raw)
n_genes_raw <- nrow(counts_raw)

# --- Gene ID detection ---
detect_gene_id_type <- function(ids) {
  ids_clean <- sub("\\.\\d+$", "", ids)  # strip version suffixes
  n <- length(ids_clean)
  pct_ens    <- sum(grepl("^ENS", ids_clean)) / n   # covers all species (ENSG/ENSMUSG/ENSRNOG/ENSCAFG/...)
  pct_entrez <- sum(grepl("^\\d+$", ids_clean)) / n
  if (pct_ens > 0.5)     return("ENSEMBL")
  if (pct_entrez > 0.5)  return("ENTREZID")
  return("SYMBOL")
}

gene_ids_orig <- rownames(counts_raw)
gene_ids_clean <- sub("\\.\\d+$", "", gene_ids_orig)
id_type <- detect_gene_id_type(gene_ids_clean)
message("[01_data_prep] Detected gene ID type: ", id_type)

# --- Load annotation DB ---
if (species %in% c("human", "hs", "homo_sapiens")) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  orgdb <- org.Hs.eg.db
} else if (species %in% c("mouse", "mm", "mus_musculus")) {
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  orgdb <- org.Mm.eg.db
} else if (species %in% c("rat", "rn", "rattus_norvegicus")) {
  suppressPackageStartupMessages(library(org.Rn.eg.db))
  orgdb <- org.Rn.eg.db
} else if (species %in% c("dog", "cfa", "beagle", "canis_familiaris", "canis_lupus_familiaris")) {
  suppressPackageStartupMessages(library(org.Cf.eg.db))
  orgdb <- org.Cf.eg.db
} else {
  stop("Unsupported species: ", species,
       ". Supported: 'human', 'mouse', 'rat', 'dog'.")
}

# --- Build ID mapping ---
message("[01_data_prep] Building ID mapping...")
if (id_type == "ENSEMBL") {
  ensembl_ids <- gene_ids_clean
  symbols   <- tryCatch(
    mapIds(orgdb, keys = ensembl_ids, column = "SYMBOL",   keytype = "ENSEMBL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(ensembl_ids)), ensembl_ids)
  )
  entrez_ids <- tryCatch(
    mapIds(orgdb, keys = ensembl_ids, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(ensembl_ids)), ensembl_ids)
  )
  id_mapping <- data.frame(
    ensembl_id = ensembl_ids,
    symbol     = unname(symbols),
    entrez_id  = unname(entrez_ids),
    stringsAsFactors = FALSE
  )
  rownames(counts_raw) <- ensembl_ids
} else if (id_type == "SYMBOL") {
  ensembl_ids <- tryCatch(
    mapIds(orgdb, keys = gene_ids_clean, column = "ENSEMBL",  keytype = "SYMBOL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(gene_ids_clean)), gene_ids_clean)
  )
  entrez_ids <- tryCatch(
    mapIds(orgdb, keys = gene_ids_clean, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(gene_ids_clean)), gene_ids_clean)
  )
  pct_mapped <- mean(!is.na(ensembl_ids))
  message(sprintf("[01_data_prep] Symbol->ENSEMBL mapping: %.1f%% mapped", pct_mapped * 100))
  if (pct_mapped < 0.5) {
    warning("Less than 50% of SYMBOL IDs mapped to ENSEMBL. Proceeding with symbols as rownames.")
    ensembl_ids_final <- gene_ids_clean
  } else {
    ensembl_ids_final <- ifelse(!is.na(ensembl_ids), unname(ensembl_ids), gene_ids_clean)
  }
  id_mapping <- data.frame(
    ensembl_id = ifelse(!is.na(ensembl_ids), unname(ensembl_ids), NA),
    symbol     = gene_ids_clean,
    entrez_id  = unname(entrez_ids),
    stringsAsFactors = FALSE
  )
  rownames(counts_raw) <- ensembl_ids_final
} else {  # ENTREZID
  entrez_ids_int <- as.integer(sub("\\.0$", "", gene_ids_clean))
  ensembl_ids <- tryCatch(
    mapIds(orgdb, keys = as.character(entrez_ids_int), column = "ENSEMBL", keytype = "ENTREZID", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(entrez_ids_int)), as.character(entrez_ids_int))
  )
  symbols <- tryCatch(
    mapIds(orgdb, keys = as.character(entrez_ids_int), column = "SYMBOL",   keytype = "ENTREZID", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(entrez_ids_int)), as.character(entrez_ids_int))
  )
  ensembl_final <- ifelse(!is.na(ensembl_ids), unname(ensembl_ids), as.character(entrez_ids_int))
  id_mapping <- data.frame(
    ensembl_id = unname(ensembl_ids),
    symbol     = unname(symbols),
    entrez_id  = as.character(entrez_ids_int),
    stringsAsFactors = FALSE
  )
  rownames(counts_raw) <- ensembl_final
}

# Remove duplicate row names (keep first)
dup_rows <- duplicated(rownames(counts_raw))
if (any(dup_rows)) {
  message("[01_data_prep] Removing ", sum(dup_rows), " duplicate gene IDs")
  counts_raw <- counts_raw[!dup_rows, , drop = FALSE]
  id_mapping <- id_mapping[!dup_rows, , drop = FALSE]
}

# --- Filter low-count genes ---
min_samples_required <- ceiling(n_samples * min_frac)
keep <- rowSums(counts_raw >= min_count) >= min_samples_required
counts_filtered <- counts_raw[keep, , drop = FALSE]
id_mapping_filtered <- id_mapping[keep, , drop = FALSE]
n_genes_filtered <- nrow(counts_filtered)
message(sprintf("[01_data_prep] Kept %d/%d genes after filtering (>=%d counts in >=%d samples)",
                n_genes_filtered, n_genes_raw, min_count, min_samples_required))

# --- Construct DDS ---
# Use first contrast variable for initial DDS (or first column of metadata)
contrasts <- config$contrasts
if (!is.null(contrasts) && length(contrasts) > 0) {
  design_var <- contrasts[[1]]$variable
} else {
  design_var <- colnames(metadata)[1]
}

metadata[[design_var]] <- factor(metadata[[design_var]])
design_formula <- as.formula(paste0("~ ", design_var))

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = metadata,
  design    = design_formula
)

# --- Write outputs ---
counts_out_path  <- file.path(output_dir, "counts_ensembl.csv")
id_mapping_path  <- file.path(output_dir, "id_mapping.csv")
dds_path         <- file.path(output_dir, "dds.rds")

write.csv(as.data.frame(counts_filtered), counts_out_path)
write.csv(id_mapping_filtered, id_mapping_path, row.names = FALSE)
saveRDS(dds, dds_path)

# --- Write summary JSON ---
summary_list <- list(
  n_samples          = n_samples,
  n_genes_raw        = n_genes_raw,
  n_genes_filtered   = n_genes_filtered,
  gene_id_type       = id_type,
  species            = species,
  design_variable    = design_var,
  min_count_threshold = min_count,
  min_samples_fraction = min_frac,
  sample_ids         = colnames(counts_filtered),
  outputs            = list(
    counts_ensembl = counts_out_path,
    id_mapping     = id_mapping_path,
    dds_rds        = dds_path
  )
)
write(toJSON(summary_list, auto_unbox = TRUE, pretty = TRUE),
      file.path(output_dir, "data_prep_summary.json"))

message("[01_data_prep] Done. Outputs written to: ", output_dir)
