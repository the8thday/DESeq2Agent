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
  library(ggsci)
})

# ── Publication-quality theme ────────────────────────────────────────────────────
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title        = element_text(size = base_size + 1, face = "bold",
                                       hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle     = element_text(size = base_size - 1.5, hjust = 0.5,
                                       color = "grey45", margin = margin(b = 8)),
      axis.title        = element_text(size = base_size, face = "bold"),
      axis.text         = element_text(size = base_size - 1.5, color = "grey20"),
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      legend.text       = element_text(size = base_size - 2),
      legend.background = element_rect(fill = alpha("white", 0.85),
                                       color = "grey80", linewidth = 0.3),
      legend.key        = element_rect(fill = "white"),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.35),
      panel.border      = element_rect(color = "grey55", fill = NA, linewidth = 0.7),
      strip.background  = element_rect(fill = "grey95", color = "grey60"),
      strip.text        = element_text(face = "bold", size = base_size - 1)
    )
}

# DE color constants — ggsci NPG palette
COL_UP  <- pal_npg("nrc")(4)[1]   # #E64B35 coral red
COL_DN  <- pal_npg("nrc")(4)[4]   # #3C5488 navy blue
COL_NS  <- "#BDBDBD"               # light gray

`%||%` <- function(a, b) if (!is.null(a) && !identical(a, NA)) a else b

args <- commandArgs(trailingOnly = TRUE)
config <- fromJSON(args[1], simplifyDataFrame = FALSE)

counts_file       <- file.path(config$output_dir, "counts_ensembl.csv")
metadata_file     <- config$metadata_file
output_dir        <- config$output_dir
contrasts_cfg     <- config$contrasts
samples_to_remove <- as.character(unlist(config$samples_to_remove %||% list()))
padj_thr          <- as.numeric(config$padj_threshold %||% 0.05)
lfc_thr           <- as.numeric(config$lfc_threshold %||% 1.0)

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
  counts_mat   <- counts_mat[, keep_samples, drop = FALSE]
  metadata     <- metadata[keep_samples, , drop = FALSE]
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

  subset_col <- ct$subset_column %||% NULL
  subset_val <- ct$subset_value  %||% NULL

  if (!is.null(subset_col) && !is.null(subset_val) && subset_col %in% colnames(metadata)) {
    message(sprintf("[03_de]   Subsetting: %s == %s", subset_col, subset_val))
    rel_samples <- rownames(metadata)[
      metadata[[ct_variable]] %in% c(ct_treat, ct_ctrl) &
      metadata[[subset_col]]  == subset_val
    ]
  } else {
    rel_samples <- rownames(metadata)[metadata[[ct_variable]] %in% c(ct_treat, ct_ctrl)]
  }

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

  # LFC shrinkage with apeglm → normal → raw fallback
  res <- tryCatch({
    coef_name <- resultsNames(dds_ct)
    coef_name <- coef_name[grepl(ct_variable, coef_name) & coef_name != "Intercept"]
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

  # Ensure symbol column always exists
  id_mapping_path <- file.path(output_dir, "id_mapping.csv")
  if (file.exists(id_mapping_path)) {
    id_map <- read.csv(id_mapping_path, stringsAsFactors = FALSE)
    res_df <- left_join(res_df, id_map, by = c("gene_id" = "ensembl_id"))
  }
  if (!"symbol" %in% colnames(res_df)) res_df$symbol <- NA_character_

  # Significant genes
  sig_up   <- res_df %>% filter(!is.na(padj), padj < padj_thr, log2FoldChange >= lfc_thr)
  sig_down <- res_df %>% filter(!is.na(padj), padj < padj_thr, log2FoldChange <= -lfc_thr)

  # Write CSVs
  write.csv(res_df,   file.path(output_dir, paste0(ct_name, "_results.csv")),  row.names = FALSE)
  write.csv(sig_up,   file.path(output_dir, paste0(ct_name, "_sig_up.csv")),   row.names = FALSE)
  write.csv(sig_down, file.path(output_dir, paste0(ct_name, "_sig_down.csv")), row.names = FALSE)

  # ── Volcano plot ──────────────────────────────────────────────────────────────
  volcano_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      significance   = case_when(
        padj < padj_thr & log2FoldChange >= lfc_thr  ~ "Up",
        padj < padj_thr & log2FoldChange <= -lfc_thr ~ "Down",
        TRUE ~ "NS"
      ),
      neg_log10_padj = -log10(pmax(padj, 1e-300)),
      display_name   = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id),
      label          = ifelse(significance != "NS" & rank(padj) <= 20, display_name, "")
    )

  n_up   <- sum(volcano_df$significance == "Up")
  n_down <- sum(volcano_df$significance == "Down")

  ns_v  <- filter(volcano_df, significance == "NS")
  sig_v <- filter(volcano_df, significance != "NS")
  lab_v <- filter(volcano_df, label != "")

  p_volcano <- ggplot() +
    # NS points — small and muted (bottom layer)
    geom_point(data = ns_v,
               aes(x = log2FoldChange, y = neg_log10_padj),
               color = COL_NS, size = 0.7, alpha = 0.38) +
    # Significant points — larger and colored (top layer)
    geom_point(data = sig_v,
               aes(x = log2FoldChange, y = neg_log10_padj, color = significance),
               size = 1.8, alpha = 0.80) +
    # Gene labels
    geom_text_repel(data = lab_v,
                    aes(x = log2FoldChange, y = neg_log10_padj,
                        label = label, color = significance),
                    size = 2.7, max.overlaps = 20, show.legend = FALSE,
                    min.segment.length = 0, segment.size = 0.3,
                    segment.alpha = 0.55, box.padding = 0.5,
                    point.padding = 0.3, force = 1.2) +
    # Threshold reference lines
    geom_vline(xintercept = c(-lfc_thr, lfc_thr),
               linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_thr),
               linetype = "dashed", color = "grey50", linewidth = 0.5) +
    # Count annotations in corners
    annotate("text", x = Inf,  y = Inf, hjust = 1.10, vjust = 1.6,
             label = paste0("\u2191 ", n_up, " up"),
             color = COL_UP, fontface = "bold", size = 3.5) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 1.6,
             label = paste0(n_down, " down \u2193"),
             color = COL_DN, fontface = "bold", size = 3.5) +
    scale_color_manual(
      values = c("Up" = COL_UP, "Down" = COL_DN),
      name   = "",
      labels = c("Up-regulated", "Down-regulated")
    ) +
    labs(
      title    = paste("Volcano Plot:", ct_name),
      subtitle = sprintf("padj < %s  |  |log\u2082FC| \u2265 %s", padj_thr, lfc_thr),
      x = expression(log[2] ~ "Fold Change"),
      y = expression(-log[10] ~ "(adjusted p-value)")
    ) +
    theme_pub() +
    theme(legend.position = "bottom",
          legend.margin   = margin(t = -4))

  ggsave(file.path(de_plots_dir, paste0(ct_name, "_volcano.png")), p_volcano,
         width = 7, height = 6, dpi = 200)

  # ── MA plot ───────────────────────────────────────────────────────────────────
  ma_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange), !is.na(baseMean)) %>%
    mutate(
      direction  = case_when(
        padj < padj_thr & log2FoldChange >= lfc_thr  ~ "Up",
        padj < padj_thr & log2FoldChange <= -lfc_thr ~ "Down",
        TRUE ~ "NS"
      ),
      significant   = direction != "NS",
      display_name  = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id),
      label         = ifelse(significant & rank(padj) <= 15, display_name, ""),
      log2_mean     = log2(baseMean + 1)
    )

  ns_m  <- filter(ma_df, direction == "NS")
  sig_m <- filter(ma_df, direction != "NS")
  lab_m <- filter(ma_df, label != "")

  p_ma <- ggplot() +
    geom_point(data = ns_m,
               aes(x = log2_mean, y = log2FoldChange),
               color = COL_NS, size = 0.7, alpha = 0.38) +
    geom_point(data = sig_m,
               aes(x = log2_mean, y = log2FoldChange, color = direction),
               size = 1.8, alpha = 0.80) +
    geom_text_repel(data = lab_m,
                    aes(x = log2_mean, y = log2FoldChange,
                        label = label, color = direction),
                    size = 2.7, max.overlaps = 15, show.legend = FALSE,
                    min.segment.length = 0, segment.size = 0.3,
                    segment.alpha = 0.55, box.padding = 0.5, point.padding = 0.3) +
    geom_hline(yintercept = 0,
               color = "grey30", linewidth = 0.8) +
    geom_hline(yintercept = c(-lfc_thr, lfc_thr),
               linetype = "dashed", color = "grey50", linewidth = 0.5) +
    scale_color_manual(
      values = c("Up" = COL_UP, "Down" = COL_DN),
      name   = "",
      labels = c("Up-regulated", "Down-regulated")
    ) +
    labs(
      title    = paste("MA Plot:", ct_name),
      subtitle = sprintf("padj < %s  |  |log\u2082FC| \u2265 %s", padj_thr, lfc_thr),
      x = expression(log[2] ~ "(Mean Expression)"),
      y = expression(log[2] ~ "Fold Change")
    ) +
    theme_pub() +
    theme(legend.position = "bottom",
          legend.margin   = margin(t = -4))

  ggsave(file.path(de_plots_dir, paste0(ct_name, "_ma.png")), p_ma,
         width = 7, height = 6, dpi = 200)

  # Store summary info
  top20_up   <- head(sig_up[, c("gene_id",
                                 intersect(c("symbol", "log2FoldChange", "padj"),
                                           colnames(sig_up)))], 20)
  top20_down <- head(sig_down[, c("gene_id",
                                   intersect(c("symbol", "log2FoldChange", "padj"),
                                             colnames(sig_down)))], 20)

  de_results_list[[ct_name]] <- list(
    contrast_name    = ct_name,
    treatment        = ct_treat,
    control          = ct_ctrl,
    n_sig_up         = nrow(sig_up),
    n_sig_down       = nrow(sig_down),
    n_total_tested   = nrow(res_df),
    top20_up_genes   = top20_up,
    top20_down_genes = top20_down,
    padj_threshold   = padj_thr,
    lfc_threshold    = lfc_thr,
    plots = list(
      volcano = file.path(de_plots_dir, paste0(ct_name, "_volcano.png")),
      ma      = file.path(de_plots_dir, paste0(ct_name, "_ma.png"))
    )
  )
}

write(toJSON(de_results_list, auto_unbox = TRUE, pretty = TRUE, null = "null"),
      file.path(output_dir, "de_summary.json"))

message("[03_de] Done. DE results written to: ", output_dir)
