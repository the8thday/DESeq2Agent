#!/usr/bin/env Rscript
# 02_qc_analysis.R - QC analysis: PCA, correlation, dispersion, library size plots
# Input:  config.json
# Output: qc_metrics.json, qc_plots/*.png

suppressPackageStartupMessages({
  library(jsonlite)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

`%||%` <- function(a, b) if (!is.null(a) && !identical(a, NA)) a else b

args <- commandArgs(trailingOnly = TRUE)
config <- fromJSON(args[1], simplifyDataFrame = FALSE)

output_dir   <- config$output_dir
dds_path     <- file.path(output_dir, "dds.rds")
metadata_file <- config$metadata_file

qc_plots_dir <- file.path(output_dir, "qc_plots")
dir.create(qc_plots_dir, showWarnings = FALSE, recursive = TRUE)

message("[02_qc] Loading DDS from: ", dds_path)
dds <- readRDS(dds_path)
metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)
metadata <- metadata[colnames(dds), , drop = FALSE]

# --- VST/rlog for PCA & correlation ---
message("[02_qc] Computing VST...")
vsd <- tryCatch(
  vst(dds, blind = TRUE),
  error = function(e) {
    message("VST failed, falling back to rlog: ", e$message)
    rlog(dds, blind = TRUE)
  }
)
vsd_mat <- assay(vsd)

# --- Library sizes ---
lib_sizes <- colSums(counts(dds))
detected_genes <- colSums(counts(dds) > 0)

# --- PCA (top 500 variable genes) ---
message("[02_qc] Computing PCA...")
rv <- rowVars(vsd_mat)
top_genes <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca_mat <- t(vsd_mat[top_genes, ])
pca_res <- prcomp(pca_mat, scale. = FALSE)
pct_var <- round(summary(pca_res)$importance[2, ] * 100, 1)

pca_df <- as.data.frame(pca_res$x[, 1:min(4, ncol(pca_res$x))])
pca_df$sample <- rownames(pca_df)
# Add first metadata column for color
first_col <- colnames(metadata)[1]
pca_df[[first_col]] <- metadata[pca_df$sample, first_col]

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[first_col]], label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  labs(
    title = "PCA Plot (Top 500 Variable Genes)",
    x = paste0("PC1: ", pct_var[1], "% variance"),
    y = paste0("PC2: ", pct_var[2], "% variance")
  ) +
  theme_bw() +
  theme(legend.title = element_text(size = 10))

ggsave(file.path(qc_plots_dir, "pca.png"), p_pca, width = 8, height = 6, dpi = 150)

# --- Correlation heatmap ---
message("[02_qc] Computing correlation heatmap...")
cor_mat <- cor(vsd_mat, method = "pearson")
annot_df <- metadata[, 1, drop = FALSE]
colnames(annot_df) <- colnames(metadata)[1]

tryCatch({
  pheatmap(
    cor_mat,
    annotation_col = annot_df,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    breaks = seq(min(cor_mat), 1, length.out = 101),
    main = "Sample Correlation Heatmap",
    filename = file.path(qc_plots_dir, "heatmap.png"),
    width = 8, height = 7
  )
}, error = function(e) message("Heatmap warning: ", e$message))

# --- MDS ---
message("[02_qc] Computing MDS...")
dist_mat <- dist(t(vsd_mat))
mds_res  <- cmdscale(dist_mat, k = 2)
mds_df   <- data.frame(
  MDS1 = mds_res[, 1], MDS2 = mds_res[, 2],
  sample = rownames(mds_res)
)
mds_df[[first_col]] <- metadata[mds_df$sample, first_col]

p_mds <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = .data[[first_col]], label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = "MDS Plot") +
  theme_bw()

ggsave(file.path(qc_plots_dir, "mds.png"), p_mds, width = 8, height = 6, dpi = 150)

# --- Dispersion plot ---
message("[02_qc] Running DESeq for dispersion plot...")
dds_for_disp <- tryCatch(
  DESeq(dds, quiet = TRUE),
  error = function(e) {
    message("DESeq() for dispersion failed: ", e$message)
    NULL
  }
)
if (!is.null(dds_for_disp)) {
  png(file.path(qc_plots_dir, "dispersion.png"), width = 800, height = 600)
  plotDispEsts(dds_for_disp, main = "Dispersion Estimates")
  dev.off()
}

# --- Library sizes plot ---
lib_df <- data.frame(
  sample = names(lib_sizes),
  library_size = lib_sizes / 1e6
)
lib_df[[first_col]] <- metadata[lib_df$sample, first_col]
lib_df$sample <- factor(lib_df$sample, levels = lib_df$sample[order(lib_df$library_size)])

p_lib <- ggplot(lib_df, aes(x = sample, y = library_size, fill = .data[[first_col]])) +
  geom_col() +
  coord_flip() +
  labs(title = "Library Sizes", x = "Sample", y = "Million Reads") +
  theme_bw()

ggsave(file.path(qc_plots_dir, "library_sizes.png"), p_lib, width = 8, height = 6, dpi = 150)

# --- Detected genes plot ---
det_df <- data.frame(
  sample = names(detected_genes),
  n_genes = detected_genes
)
det_df[[first_col]] <- metadata[det_df$sample, first_col]
det_df$sample <- factor(det_df$sample, levels = det_df$sample[order(det_df$n_genes)])

p_det <- ggplot(det_df, aes(x = sample, y = n_genes, fill = .data[[first_col]])) +
  geom_col() +
  coord_flip() +
  labs(title = "Detected Genes per Sample", x = "Sample", y = "Number of Genes Detected") +
  theme_bw()

ggsave(file.path(qc_plots_dir, "detected_genes.png"), p_det, width = 8, height = 6, dpi = 150)

# --- IQR-based outlier detection ---
detect_outliers_iqr <- function(metric_vec, names_vec, metric_name, threshold = 1.5) {
  q1 <- quantile(metric_vec, 0.25)
  q3 <- quantile(metric_vec, 0.75)
  iqr <- q3 - q1
  lower <- q1 - threshold * iqr
  upper <- q3 + threshold * iqr
  outlier_idx <- which(metric_vec < lower | metric_vec > upper)
  lapply(outlier_idx, function(i) list(
    sample    = names_vec[i],
    metric    = metric_name,
    value     = metric_vec[i],
    threshold = list(lower = lower, upper = upper)
  ))
}

outlier_flags <- c(
  detect_outliers_iqr(lib_sizes, names(lib_sizes), "library_size"),
  detect_outliers_iqr(detected_genes, names(detected_genes), "detected_genes"),
  detect_outliers_iqr(cor_mat[lower.tri(cor_mat)],
                      rep("pairwise", sum(lower.tri(cor_mat))),
                      "sample_correlation")
)
# PCA outliers: samples far from centroid
pc1_vals <- pca_res$x[, 1]
pc2_vals <- pca_res$x[, 2]
pca_dist <- sqrt((pc1_vals - mean(pc1_vals))^2 + (pc2_vals - mean(pc2_vals))^2)
outlier_flags <- c(outlier_flags,
                   detect_outliers_iqr(pca_dist, names(pca_dist), "pca_distance"))

# Deduplicate by sample name
outlier_samples_flagged <- unique(sapply(outlier_flags, function(x) x$sample))
outlier_samples_flagged <- outlier_samples_flagged[outlier_samples_flagged != "pairwise"]

# --- Write qc_metrics.json ---
qc_metrics <- list(
  n_samples           = ncol(dds),
  n_genes_filtered    = nrow(dds),
  library_sizes       = as.list(lib_sizes),
  detected_genes      = as.list(detected_genes),
  pca_variance_pct    = as.list(pct_var[1:min(4, length(pct_var))]),
  sample_correlation  = lapply(rownames(cor_mat), function(s) {
    list(sample = s, correlations = as.list(setNames(cor_mat[s, ], colnames(cor_mat))))
  }),
  outlier_flags       = outlier_flags,
  outlier_samples_flagged = outlier_samples_flagged,
  plots               = list(
    pca            = file.path(qc_plots_dir, "pca.png"),
    heatmap        = file.path(qc_plots_dir, "heatmap.png"),
    mds            = file.path(qc_plots_dir, "mds.png"),
    dispersion     = file.path(qc_plots_dir, "dispersion.png"),
    library_sizes  = file.path(qc_plots_dir, "library_sizes.png"),
    detected_genes = file.path(qc_plots_dir, "detected_genes.png")
  )
)

write(toJSON(qc_metrics, auto_unbox = TRUE, pretty = TRUE, null = "null"),
      file.path(output_dir, "qc_metrics.json"))

message("[02_qc] Done. QC metrics written.")
