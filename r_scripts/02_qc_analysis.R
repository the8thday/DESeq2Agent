#!/usr/bin/env Rscript
# 02_qc_analysis.R - QC analysis: PCA, correlation, dispersion, library size plots
# Input:  config.json
# Output: qc_metrics.json, qc_plots/*.pdf

suppressPackageStartupMessages({
  library(jsonlite)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(ggsci)
  library(ggdendro)
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
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.key        = element_rect(fill = "white"),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.35),
      panel.border      = element_rect(color = "grey55", fill = NA, linewidth = 0.7),
      strip.background  = element_rect(fill = "grey95", color = "grey60"),
      strip.text        = element_text(face = "bold", size = base_size - 1)
    )
}

# ggsci NPG palette (up to 10 groups)
GROUP_COLORS <- pal_npg("nrc")(10)

# Save ggplot as PDF (vector format, publication-quality)
save_gg <- function(path, plot_obj, ...) {
  ggsave(path, plot_obj, ...)
}

`%||%` <- function(a, b) if (!is.null(a) && !identical(a, NA)) a else b

args   <- commandArgs(trailingOnly = TRUE)
config <- fromJSON(args[1], simplifyDataFrame = FALSE)

output_dir    <- config$output_dir
dds_path      <- file.path(output_dir, "dds.rds")
metadata_file <- config$metadata_file

qc_plots_dir <- file.path(output_dir, "qc_plots")
dir.create(qc_plots_dir, showWarnings = FALSE, recursive = TRUE)

message("[02_qc] Loading DDS from: ", dds_path)
dds      <- readRDS(dds_path)
metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)
metadata <- metadata[colnames(dds), , drop = FALSE]

# Auto-detect informative grouping columns
n_samp <- ncol(dds)
color_cols <- Filter(function(col) {
  n_uniq <- length(unique(metadata[[col]]))
  n_uniq >= 2 && n_uniq <= 10 && n_uniq < n_samp
}, colnames(metadata))
if (length(color_cols) == 0) color_cols <- colnames(metadata)[1]  # fallback
for (col in color_cols) metadata[[col]] <- as.character(metadata[[col]])
primary_col <- color_cols[1]   # used for heatmap, library sizes, detected genes, hclust

# ── VST normalization ────────────────────────────────────────────────────────────
message("[02_qc] Computing VST...")
vsd <- tryCatch(
  vst(dds, blind = TRUE),
  error = function(e) {
    message("VST failed, falling back to rlog: ", e$message)
    rlog(dds, blind = TRUE)
  }
)
vsd_mat <- assay(vsd)

# ── Library sizes & detected genes ──────────────────────────────────────────────
lib_sizes      <- colSums(counts(dds))
detected_genes <- colSums(counts(dds) > 0)

# ── PCA (top 2000 variable genes) ────────────────────────────────────────────────
message("[02_qc] Computing PCA...")
rv      <- rowVars(vsd_mat)
top_idx <- order(rv, decreasing = TRUE)[seq_len(min(2000, length(rv)))]
pca_res <- prcomp(t(vsd_mat[top_idx, ]), scale. = FALSE)
pct_var <- round(summary(pca_res)$importance[2, ] * 100, 1)

pca_df <- as.data.frame(pca_res$x[, 1:min(4, ncol(pca_res$x))])
pca_df$sample <- rownames(pca_df)

for (col in color_cols) {
  pca_df[[col]] <- metadata[pca_df$sample, col]
  min_grp_size  <- min(table(pca_df[[col]]))
  draw_ellipse  <- min_grp_size >= 3

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                                color = .data[[col]],
                                label = sample)) +
    { if (draw_ellipse) stat_ellipse(aes(group = .data[[col]]),
                 level = 0.9, type = "t",
                 linetype = "dashed", linewidth = 0.6, show.legend = FALSE) } +
    geom_point(size = 4, alpha = 0.88) +
    geom_text_repel(size = 3.2, max.overlaps = 20, show.legend = FALSE,
                    min.segment.length = 0, segment.size = 0.3,
                    segment.alpha = 0.6, box.padding = 0.45, point.padding = 0.3) +
    scale_color_manual(values = GROUP_COLORS, name = col) +
    labs(
      title    = "PCA — Top 2000 Variable Genes",
      subtitle = sprintf("PC1 explains %s%% | PC2 explains %s%% of variance",
                         pct_var[1], pct_var[2]),
      x = paste0("PC1  (", pct_var[1], "%)"),
      y = paste0("PC2  (", pct_var[2], "%)")
    ) +
    theme_pub()

  save_gg(file.path(qc_plots_dir, paste0("pca_", col, ".pdf")), p_pca,
          width = 7, height = 5.5)
}

# ── Sample correlation heatmap ───────────────────────────────────────────────────
message("[02_qc] Computing correlation heatmap...")
cor_mat  <- cor(vsd_mat, method = "pearson")
annot_df <- metadata[, primary_col, drop = FALSE]

n_grp     <- length(unique(annot_df[[primary_col]]))
grp_cols  <- setNames(GROUP_COLORS[seq_len(n_grp)], sort(unique(annot_df[[primary_col]])))
annot_colors <- setNames(list(grp_cols), primary_col)

cor_min <- max(0.5, floor(min(cor_mat) * 100) / 100 - 0.01)
tryCatch({
  pheatmap(
    cor_mat,
    annotation_col    = annot_df,
    annotation_colors = annot_colors,
    color             = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
    breaks            = seq(cor_min, 1.0, length.out = 101),
    main              = "Sample-to-Sample Correlation (Pearson, VST)",
    fontsize          = 11,
    fontsize_row      = 10,
    fontsize_col      = 10,
    border_color      = "white",
    treeheight_row    = 35,
    treeheight_col    = 35,
    filename          = file.path(qc_plots_dir, "heatmap.pdf"),
    width = 8, height = 7
  )
}, error = function(e) message("Heatmap warning: ", e$message))

# ── MDS ─────────────────────────────────────────────────────────────────────────
message("[02_qc] Computing MDS...")
dist_mat <- dist(t(vsd_mat))
mds_res  <- cmdscale(dist_mat, k = 2)
mds_df   <- data.frame(MDS1 = mds_res[, 1], MDS2 = mds_res[, 2],
                        sample = rownames(mds_res))

for (col in color_cols) {
  mds_df[[col]] <- metadata[mds_df$sample, col]
  min_grp_size  <- min(table(mds_df[[col]]))
  draw_ellipse  <- min_grp_size >= 3

  p_mds <- ggplot(mds_df, aes(x = MDS1, y = MDS2,
                                color = .data[[col]],
                                label = sample)) +
    { if (draw_ellipse) stat_ellipse(aes(group = .data[[col]]),
                 level = 0.9, type = "t",
                 linetype = "dashed", linewidth = 0.6, show.legend = FALSE) } +
    geom_point(size = 4, alpha = 0.88) +
    geom_text_repel(size = 3.2, max.overlaps = 20, show.legend = FALSE,
                    min.segment.length = 0, segment.size = 0.3,
                    segment.alpha = 0.6, box.padding = 0.45, point.padding = 0.3) +
    scale_color_manual(values = GROUP_COLORS, name = col) +
    labs(
      title    = "Multidimensional Scaling (MDS)",
      subtitle = "Based on VST-normalized Euclidean distances",
      x = "MDS Dimension 1", y = "MDS Dimension 2"
    ) +
    theme_pub()

  save_gg(file.path(qc_plots_dir, paste0("mds_", col, ".pdf")), p_mds,
          width = 7, height = 5.5)
}

# ── Hierarchical clustering ──────────────────────────────────────────────────────
message("[02_qc] Computing hierarchical clustering...")
hc        <- hclust(dist_mat, method = "complete")
dend_data <- dendro_data(hc, type = "rectangle")
label_df  <- label(dend_data)
label_df[[primary_col]] <- metadata[as.character(label_df$label), primary_col]

p_hclust <- ggplot() +
  geom_segment(data = segment(dend_data),
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey35", linewidth = 0.55) +
  geom_text(data = label_df,
            aes(x = x, y = 0, label = label,
                color = .data[[primary_col]]),
            hjust = 1, angle = 90, size = 3.0, show.legend = TRUE) +
  scale_color_manual(
    values = setNames(GROUP_COLORS[seq_len(n_grp)],
                      sort(unique(label_df[[primary_col]]))),
    name = primary_col
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.35, 0.05))) +
  labs(
    title    = "Sample Hierarchical Clustering",
    subtitle = "Complete linkage \u00b7 Euclidean distance on VST counts",
    x = NULL, y = "Distance"
  ) +
  theme_pub() +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_gg(file.path(qc_plots_dir, "hclust.pdf"), p_hclust,
        width = max(6, ncol(dds) * 0.45), height = 6)

# ── Dispersion (custom ggplot) ───────────────────────────────────────────────────
message("[02_qc] Running DESeq for dispersion plot...")
dds_for_disp <- tryCatch(
  DESeq(dds, quiet = TRUE),
  error = function(e) { message("DESeq() for dispersion failed: ", e$message); NULL }
)

if (!is.null(dds_for_disp)) {
  disp_df <- tryCatch({
    df <- data.frame(
      mean_norm  = rowMeans(counts(dds_for_disp, normalized = TRUE)),
      disp_est   = mcols(dds_for_disp)$dispGeneEst,
      disp_fit   = mcols(dds_for_disp)$dispFit,
      disp_final = dispersions(dds_for_disp),
      outlier    = mcols(dds_for_disp)$dispOutlier %||% FALSE
    )
    filter(df, !is.na(disp_est), mean_norm > 0, disp_est > 0)
  }, error = function(e) NULL)

  if (!is.null(disp_df) && nrow(disp_df) > 0) {
    disp_fit_df <- filter(disp_df, !is.na(disp_fit)) %>% arrange(mean_norm)

    p_disp <- ggplot(disp_df, aes(x = mean_norm)) +
      geom_point(data = filter(disp_df, !outlier),
                 aes(y = disp_est),
                 color = "#999999", size = 0.55, alpha = 0.45) +
      geom_point(data = filter(disp_df, outlier),
                 aes(y = disp_est),
                 color = "#E64B35", size = 1.1, alpha = 0.8, shape = 17) +
      geom_point(aes(y = disp_final),
                 color = "#3C5488", size = 0.7, alpha = 0.55) +
      geom_line(data = disp_fit_df,
                aes(y = disp_fit),
                color = "#E64B35", linewidth = 1.1) +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title    = "Dispersion Estimates",
        subtitle = "Grey \u25cf gene-wise  |  Red line: fitted trend  |  Blue \u25cf final (MAP)  |  Red \u25b2 outliers",
        x = "Mean of Normalized Counts",
        y = "Dispersion"
      ) +
      theme_pub()

    save_gg(file.path(qc_plots_dir, "dispersion.pdf"), p_disp,
            width = 8, height = 5.5)
  } else {
    pdf(file.path(qc_plots_dir, "dispersion.pdf"), width = 8, height = 6)
    plotDispEsts(dds_for_disp, main = "Dispersion Estimates")
    dev.off()
  }
}

# ── Library sizes ────────────────────────────────────────────────────────────────
lib_df <- data.frame(sample = names(lib_sizes), library_size = lib_sizes / 1e6)
lib_df[[primary_col]] <- metadata[lib_df$sample, primary_col]
lib_df$sample <- factor(lib_df$sample, levels = lib_df$sample[order(lib_df$library_size)])

p_lib <- ggplot(lib_df, aes(x = sample, y = library_size,
                              fill = .data[[primary_col]])) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f M", library_size)),
            hjust = -0.1, size = 3.2, color = "grey25") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = GROUP_COLORS, name = primary_col) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title    = "Library Sizes",
    subtitle = "Total mapped reads per sample",
    x = NULL, y = "Million Reads"
  ) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

save_gg(file.path(qc_plots_dir, "library_sizes.pdf"), p_lib,
        width = 7, height = 5)

# ── Detected genes ───────────────────────────────────────────────────────────────
det_df <- data.frame(sample = names(detected_genes), n_genes = detected_genes)
det_df[[primary_col]] <- metadata[det_df$sample, primary_col]
det_df$sample <- factor(det_df$sample, levels = det_df$sample[order(det_df$n_genes)])

p_det <- ggplot(det_df, aes(x = sample, y = n_genes,
                              fill = .data[[primary_col]])) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = format(n_genes, big.mark = ",")),
            hjust = -0.1, size = 3.2, color = "grey25") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = GROUP_COLORS, name = primary_col) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title    = "Detected Genes per Sample",
    subtitle = "Genes with at least 1 raw count",
    x = NULL, y = "Number of Genes Detected"
  ) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

save_gg(file.path(qc_plots_dir, "detected_genes.pdf"), p_det,
        width = 7, height = 5)

# ── IQR-based outlier detection ──────────────────────────────────────────────────
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
  detect_outliers_iqr(lib_sizes,      names(lib_sizes),      "library_size"),
  detect_outliers_iqr(detected_genes, names(detected_genes), "detected_genes"),
  detect_outliers_iqr(cor_mat[lower.tri(cor_mat)],
                      rep("pairwise", sum(lower.tri(cor_mat))),
                      "sample_correlation")
)
pc1_vals <- pca_res$x[, 1]
pc2_vals <- pca_res$x[, 2]
pca_dist <- sqrt((pc1_vals - mean(pc1_vals))^2 + (pc2_vals - mean(pc2_vals))^2)
outlier_flags <- c(outlier_flags,
                   detect_outliers_iqr(pca_dist, names(pca_dist), "pca_distance"))

outlier_samples_flagged <- unique(sapply(outlier_flags, function(x) x$sample))
outlier_samples_flagged <- outlier_samples_flagged[outlier_samples_flagged != "pairwise"]

# ── Write qc_metrics.json ────────────────────────────────────────────────────────
qc_metrics <- list(
  n_samples               = ncol(dds),
  n_genes_filtered        = nrow(dds),
  library_sizes           = as.list(lib_sizes),
  detected_genes          = as.list(detected_genes),
  pca_variance_pct        = as.list(pct_var[1:min(4, length(pct_var))]),
  sample_correlation      = lapply(rownames(cor_mat), function(s) {
    list(sample = s, correlations = as.list(setNames(cor_mat[s, ], colnames(cor_mat))))
  }),
  outlier_flags           = outlier_flags,
  outlier_samples_flagged = outlier_samples_flagged,
  plots = c(
    setNames(
      lapply(color_cols, function(col) file.path(qc_plots_dir, paste0("pca_", col, ".pdf"))),
      paste0("pca_", color_cols)
    ),
    setNames(
      lapply(color_cols, function(col) file.path(qc_plots_dir, paste0("mds_", col, ".pdf"))),
      paste0("mds_", color_cols)
    ),
    list(
      heatmap        = file.path(qc_plots_dir, "heatmap.pdf"),
      hclust         = file.path(qc_plots_dir, "hclust.pdf"),
      dispersion     = file.path(qc_plots_dir, "dispersion.pdf"),
      library_sizes  = file.path(qc_plots_dir, "library_sizes.pdf"),
      detected_genes = file.path(qc_plots_dir, "detected_genes.pdf")
    )
  )
)

write(toJSON(qc_metrics, auto_unbox = TRUE, pretty = TRUE, null = "null"),
      file.path(output_dir, "qc_metrics.json"))

message("[02_qc] Done. QC metrics written.")
