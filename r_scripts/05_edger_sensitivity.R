#!/usr/bin/env Rscript
# 05_edger_sensitivity.R - edgeR sensitivity analysis for cross-method validation
# Input:  config.json (same as other R scripts)
# Output: edger_{contrast}_results.csv, edger_sensitivity.json

suppressPackageStartupMessages({
  library(jsonlite)
  library(edgeR)
  library(dplyr)
})

`%||%` <- function(a, b) if (!is.null(a) && !identical(a, NA)) a else b

args <- commandArgs(trailingOnly = TRUE)
config <- fromJSON(args[1], simplifyDataFrame = FALSE)

output_dir      <- config$output_dir
contrasts_cfg   <- config$contrasts
padj_thr        <- as.numeric(config$padj_threshold %||% 0.05)
lfc_thr         <- as.numeric(config$lfc_threshold %||% 1.0)
samples_to_rm   <- as.character(unlist(config$samples_to_remove %||% list()))

# Read metadata
metadata <- read.csv(config$metadata_file, stringsAsFactors = FALSE)
rownames(metadata) <- metadata[, 1]

# Read ENSEMBL-mapped counts (produced by 01_data_prep.R)
counts_file <- file.path(output_dir, "counts_ensembl.csv")
if (!file.exists(counts_file)) {
  stop("counts_ensembl.csv not found in ", output_dir)
}
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)

# Read id_mapping for gene symbols
id_map_path <- file.path(output_dir, "id_mapping.csv")
id_map <- if (file.exists(id_map_path)) {
  read.csv(id_map_path, stringsAsFactors = FALSE)
} else NULL

get_symbol <- function(ensembl_ids) {
  if (!is.null(id_map)) {
    idx <- match(ensembl_ids, id_map$ensembl_id)
    syms <- id_map$symbol[idx]
    syms[is.na(syms)] <- ensembl_ids[is.na(syms)]
    return(syms)
  }
  return(ensembl_ids)
}

# Remove outlier samples
if (length(samples_to_rm) > 0) {
  counts <- counts[, !colnames(counts) %in% samples_to_rm, drop = FALSE]
  metadata <- metadata[!rownames(metadata) %in% samples_to_rm, , drop = FALSE]
  message("[05_edger] Removed samples: ", paste(samples_to_rm, collapse = ", "))
}

sensitivity_results <- list()

for (ct in contrasts_cfg) {
  ct_name   <- ct$name
  variable  <- ct$variable
  treatment <- ct$treatment
  control   <- ct$control
  subset_col <- ct$subset_column
  subset_val <- ct$subset_value

  message(sprintf("[05_edger] Processing contrast: %s", ct_name))

  # Subset samples if needed (same logic as 03_de_analysis.R)
  cur_meta <- metadata
  cur_counts <- counts
  if (!is.null(subset_col) && !is.null(subset_val)) {
    keep <- cur_meta[[subset_col]] == subset_val
    cur_meta <- cur_meta[keep, , drop = FALSE]
    cur_counts <- cur_counts[, rownames(cur_meta), drop = FALSE]
    message(sprintf("[05_edger] Subset to %s=%s: %d samples",
                    subset_col, subset_val, nrow(cur_meta)))
  }

  # Filter to treatment/control samples only
  keep_samples <- cur_meta[[variable]] %in% c(treatment, control)
  cur_meta <- cur_meta[keep_samples, , drop = FALSE]
  cur_counts <- cur_counts[, rownames(cur_meta), drop = FALSE]

  if (ncol(cur_counts) < 2) {
    message("[05_edger] Skipping ", ct_name, ": fewer than 2 samples")
    next
  }

  # edgeR pipeline
  group <- factor(cur_meta[[variable]], levels = c(control, treatment))
  dge <- DGEList(counts = cur_counts, group = group)
  keep_genes <- filterByExpr(dge)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  tryCatch({
    design_mat <- model.matrix(~ group)
    dge <- estimateDisp(dge, design_mat)
    fit <- glmQLFit(dge, design_mat)
    qlf <- glmQLFTest(fit, coef = 2)

    res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
    res$gene_id <- rownames(res)
    res$symbol <- get_symbol(res$gene_id)

    # Rename columns for clarity
    edger_df <- data.frame(
      gene_id = res$gene_id,
      symbol  = res$symbol,
      logFC   = res$logFC,
      logCPM  = res$logCPM,
      F       = res$F,
      PValue  = res$PValue,
      FDR     = res$FDR,
      stringsAsFactors = FALSE
    )

    # Save edgeR results
    write.csv(edger_df, file.path(output_dir, paste0("edger_", ct_name, "_results.csv")),
              row.names = FALSE)
    message(sprintf("[05_edger] Saved edger_%s_results.csv (%d genes)", ct_name, nrow(edger_df)))

    # ── Concordance with DESeq2 ────────────────────────────────────────────
    deseq2_file <- file.path(output_dir, paste0(ct_name, "_results.csv"))
    if (file.exists(deseq2_file)) {
      deseq2_df <- read.csv(deseq2_file, stringsAsFactors = FALSE)

      # Common genes
      common_genes <- intersect(edger_df$gene_id, deseq2_df$gene_id)
      n_common <- length(common_genes)

      # Significant genes
      edger_sig <- edger_df$gene_id[edger_df$FDR < padj_thr & abs(edger_df$logFC) >= lfc_thr]
      deseq2_sig <- deseq2_df$gene_id[!is.na(deseq2_df$padj) & deseq2_df$padj < padj_thr &
                                         abs(deseq2_df$log2FoldChange) >= lfc_thr]
      overlap <- intersect(edger_sig, deseq2_sig)

      # Spearman correlation of logFC on common genes
      edger_common <- edger_df[match(common_genes, edger_df$gene_id), ]
      deseq2_common <- deseq2_df[match(common_genes, deseq2_df$gene_id), ]
      spearman_cor <- cor(edger_common$logFC, deseq2_common$log2FoldChange,
                          method = "spearman", use = "complete.obs")

      # Direction agreement on common significant genes
      common_sig <- intersect(common_genes, union(edger_sig, deseq2_sig))
      if (length(common_sig) > 0) {
        e_lfc <- edger_df$logFC[match(common_sig, edger_df$gene_id)]
        d_lfc <- deseq2_df$log2FoldChange[match(common_sig, deseq2_df$gene_id)]
        direction_agree <- sum(sign(e_lfc) == sign(d_lfc), na.rm = TRUE) / length(common_sig) * 100
      } else {
        direction_agree <- NA
      }

      sensitivity_results[[ct_name]] <- list(
        contrast_name          = ct_name,
        edger_n_sig            = length(edger_sig),
        deseq2_n_sig           = length(deseq2_sig),
        overlap_n              = length(overlap),
        overlap_pct_of_deseq2  = if (length(deseq2_sig) > 0) round(length(overlap) / length(deseq2_sig) * 100, 1) else NA,
        overlap_pct_of_edger   = if (length(edger_sig) > 0) round(length(overlap) / length(edger_sig) * 100, 1) else NA,
        spearman_logfc_cor     = round(spearman_cor, 4),
        direction_agreement_pct = round(direction_agree, 1),
        n_common_genes_tested  = n_common
      )
    } else {
      message("[05_edger] DESeq2 results not found for ", ct_name, ", skipping concordance")
      sensitivity_results[[ct_name]] <- list(
        contrast_name = ct_name,
        edger_n_sig   = length(edger_sig),
        note          = "DESeq2 results not found for concordance analysis"
      )
    }
  }, error = function(e) {
    message("[05_edger] Error processing ", ct_name, ": ", e$message)
  })
}

# Write sensitivity JSON
write(toJSON(sensitivity_results, auto_unbox = TRUE, pretty = TRUE, null = "null"),
      file.path(output_dir, "edger_sensitivity.json"))

message("[05_edger] Sensitivity analysis complete.")
