#!/usr/bin/env Rscript
# 03_de_analysis.R - Differential expression with DESeq2 + lfcShrink
# Input:  config.json (includes samples_to_remove after QC decision)
# Output: per-contrast CSV files, volcano/MA plots, de_summary.json

suppressPackageStartupMessages({
  library(jsonlite)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
})

`%||%` <- function(a, b) if (!is.null(a) && !identical(a, NA)) a else b

args <- commandArgs(trailingOnly = TRUE)
config <- fromJSON(args[1], simplifyDataFrame = FALSE)

counts_file      <- file.path(config$output_dir, "counts_ensembl.csv")
metadata_file    <- config$metadata_file
output_dir       <- config$output_dir
contrasts_cfg    <- config$contrasts
samples_to_remove <- as.character(unlist(config$samples_to_remove %||% list()))
padj_thr         <- as.numeric(config$padj_threshold %||% 0.05)
lfc_thr          <- as.numeric(config$lfc_threshold %||% 1.0)

de_plots_dir <- file.path(output_dir, "de_plots")
dir.create(de_plots_dir, showWarnings = FALSE, recursive = TRUE)

message("[03_de] Reading counts and metadata...")
counts_mat <- as.matrix(read.csv(counts_file, row.names = 1, check.names = FALSE))
storage.mode(counts_mat) <- "integer"
metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)

# Remove excluded samples
if (length(samples_to_remove) > 0) {
  message("[03_de] Removing samples: ", paste(samples_to_remove, collapse = ", "))
  keep_samples <- setdiff(colnames(counts_mat), samples_to_remove)
  counts_mat <- counts_mat[, keep_samples, drop = FALSE]
  metadata   <- metadata[keep_samples, , drop = FALSE]
}

# Ensure sample alignment
common_samples <- intersect(colnames(counts_mat), rownames(metadata))
counts_mat <- counts_mat[, common_samples, drop = FALSE]
metadata   <- metadata[common_samples, , drop = FALSE]

de_results_list <- list()

for (ct in contrasts_cfg) {
  ct_name     <- ct$name
  ct_variable <- ct$variable
  ct_treat    <- ct$treatment
  ct_ctrl     <- ct$control

  message(sprintf("[03_de] Processing contrast: %s (%s vs %s)", ct_name, ct_treat, ct_ctrl))

  # Subset to relevant samples
  rel_samples <- rownames(metadata)[metadata[[ct_variable]] %in% c(ct_treat, ct_ctrl)]
  if (length(rel_samples) < 4) {
    warning(sprintf("Contrast %s: fewer than 4 samples, skipping", ct_name))
    next
  }
  counts_sub <- counts_mat[, rel_samples, drop = FALSE]
  meta_sub   <- metadata[rel_samples, , drop = FALSE]
  meta_sub[[ct_variable]] <- factor(meta_sub[[ct_variable]], levels = c(ct_ctrl, ct_treat))

  dds_ct <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = meta_sub,
    design    = as.formula(paste0("~ ", ct_variable))
  )
  dds_ct <- DESeq(dds_ct, quiet = TRUE)

  # LFC shrinkage with apeglm
  res <- tryCatch({
    coef_name <- resultsNames(dds_ct)
    coef_name <- coef_name[grepl(ct_variable, coef_name) & !coef_name %in% "Intercept"]
    if (length(coef_name) == 0) coef_name <- tail(resultsNames(dds_ct), 1)
    lfcShrink(dds_ct, coef = coef_name[1], type = "apeglm", quiet = TRUE)
  }, error = function(e) {
    message("apeglm failed, using normal shrinkage: ", e$message)
    tryCatch(
      lfcShrink(dds_ct, coef = tail(resultsNames(dds_ct), 1), type = "normal"),
      error = function(e2) results(dds_ct)
    )
  })

  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj, desc(abs(log2FoldChange)))

  # Read id_mapping for symbol annotation
  id_mapping_path <- file.path(output_dir, "id_mapping.csv")
  if (file.exists(id_mapping_path)) {
    id_map <- read.csv(id_mapping_path, stringsAsFactors = FALSE)
    res_df <- left_join(res_df, id_map, by = c("gene_id" = "ensembl_id"))
  }

  # Significant genes
  sig_up   <- res_df %>% filter(!is.na(padj), padj < padj_thr, log2FoldChange >= lfc_thr)
  sig_down <- res_df %>% filter(!is.na(padj), padj < padj_thr, log2FoldChange <= -lfc_thr)

  # Write CSVs
  write.csv(res_df,   file.path(output_dir, paste0(ct_name, "_results.csv")),  row.names = FALSE)
  write.csv(sig_up,   file.path(output_dir, paste0(ct_name, "_sig_up.csv")),   row.names = FALSE)
  write.csv(sig_down, file.path(output_dir, paste0(ct_name, "_sig_down.csv")), row.names = FALSE)

  # --- Volcano plot ---
  volcano_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      significance = case_when(
        padj < padj_thr & log2FoldChange >= lfc_thr  ~ "Up",
        padj < padj_thr & log2FoldChange <= -lfc_thr ~ "Down",
        TRUE ~ "NS"
      ),
      neg_log10_padj = -log10(pmax(padj, 1e-300)),
      label = ifelse(significance != "NS" & rank(padj) <= 20,
                     ifelse(!is.na(symbol), symbol, gene_id), "")
    )

  p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj,
                                       color = significance, label = label)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_text_repel(size = 2.5, max.overlaps = 15, show.legend = FALSE) +
    scale_color_manual(values = c("Up" = "#E74C3C", "Down" = "#3498DB", "NS" = "#95A5A6")) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", alpha = 0.5) +
    labs(title = paste("Volcano Plot:", ct_name),
         x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "") +
    theme_bw()

  ggsave(file.path(de_plots_dir, paste0(ct_name, "_volcano.png")), p_volcano,
         width = 8, height = 6, dpi = 150)

  # --- MA plot ---
  ma_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange), !is.na(baseMean)) %>%
    mutate(
      significant = padj < padj_thr & abs(log2FoldChange) >= lfc_thr,
      label = ifelse(significant & rank(padj) <= 15,
                     ifelse(!is.na(symbol), symbol, gene_id), "")
    )

  p_ma <- ggplot(ma_df, aes(x = log2(baseMean + 1), y = log2FoldChange,
                              color = significant, label = label)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_text_repel(size = 2.5, max.overlaps = 15, show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6")) +
    geom_hline(yintercept = 0, linetype = "solid", alpha = 0.8) +
    geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", alpha = 0.5) +
    labs(title = paste("MA Plot:", ct_name),
         x = "log2(Mean Expression)", y = "log2 Fold Change", color = "Significant") +
    theme_bw()

  ggsave(file.path(de_plots_dir, paste0(ct_name, "_ma.png")), p_ma,
         width = 8, height = 6, dpi = 150)

  # Store summary info
  top20_up   <- head(sig_up[, c("gene_id", intersect(c("symbol", "log2FoldChange", "padj"), colnames(sig_up)))], 20)
  top20_down <- head(sig_down[, c("gene_id", intersect(c("symbol", "log2FoldChange", "padj"), colnames(sig_down)))], 20)

  de_results_list[[ct_name]] <- list(
    contrast_name      = ct_name,
    treatment          = ct_treat,
    control            = ct_ctrl,
    n_sig_up           = nrow(sig_up),
    n_sig_down         = nrow(sig_down),
    n_total_tested     = nrow(res_df),
    top20_up_genes     = top20_up,
    top20_down_genes   = top20_down,
    padj_threshold     = padj_thr,
    lfc_threshold      = lfc_thr,
    plots = list(
      volcano = file.path(de_plots_dir, paste0(ct_name, "_volcano.png")),
      ma      = file.path(de_plots_dir, paste0(ct_name, "_ma.png"))
    )
  )
}

write(toJSON(de_results_list, auto_unbox = TRUE, pretty = TRUE, null = "null"),
      file.path(output_dir, "de_summary.json"))

message("[03_de] Done. DE results written to: ", output_dir)
