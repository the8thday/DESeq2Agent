#!/usr/bin/env Rscript
# 04_enrichment.R - GSEA and ORA enrichment analysis via clusterProfiler
# Input:  config.json
# Output: enrichment_results_{contrast}.json, enrichment/plots/{contrast}_*.png

suppressPackageStartupMessages({
  library(jsonlite)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
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
      axis.text.y       = element_text(size = base_size - 2),
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

`%||%` <- function(a, b) if (!is.null(a) && !identical(a, NA)) a else b

args <- commandArgs(trailingOnly = TRUE)
config <- fromJSON(args[1], simplifyDataFrame = FALSE)

output_dir      <- config$output_dir
contrasts_cfg   <- config$contrasts
species         <- tolower(config$species %||% "human")
padj_thr        <- as.numeric(config$padj_threshold %||% 0.05)
lfc_thr         <- as.numeric(config$lfc_threshold %||% 1.0)
id_mapping_path <- file.path(output_dir, "id_mapping.csv")

if (species %in% c("human", "hs", "homo_sapiens")) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  orgdb    <- org.Hs.eg.db
  kegg_org <- "hsa"
} else if (species %in% c("mouse", "mm", "mus_musculus")) {
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  orgdb    <- org.Mm.eg.db
  kegg_org <- "mmu"
} else if (species %in% c("rat", "rn", "rattus_norvegicus")) {
  suppressPackageStartupMessages(library(org.Rn.eg.db))
  orgdb    <- org.Rn.eg.db
  kegg_org <- "rno"
} else if (species %in% c("dog", "cfa", "beagle", "canis_familiaris", "canis_lupus_familiaris")) {
  suppressPackageStartupMessages(library(org.Cf.eg.db))
  orgdb    <- org.Cf.eg.db
  kegg_org <- "cfa"
} else {
  stop("Unsupported species: ", species,
       ". Supported: 'human', 'mouse', 'rat', 'dog'.")
}

enrich_plots_dir <- file.path(output_dir, "enrichment", "plots")
dir.create(enrich_plots_dir, showWarnings = FALSE, recursive = TRUE)

id_map <- if (file.exists(id_mapping_path)) {
  read.csv(id_mapping_path, stringsAsFactors = FALSE)
} else NULL

# Convert ENSEMBL ids to ENTREZ — returns NA for unmapped IDs (same length as input)
ensembl_to_entrez <- function(ensembl_ids) {
  if (!is.null(id_map)) {
    idx    <- match(ensembl_ids, id_map$ensembl_id)
    entrez <- id_map$entrez_id[idx]
  } else {
    entrez <- tryCatch(
      mapIds(orgdb, keys = ensembl_ids, column = "ENTREZID",
             keytype = "ENSEMBL", multiVals = "first"),
      error = function(e) rep(NA_character_, length(ensembl_ids))
    )
  }
  as.character(entrez)  # NAs preserved; callers filter with !is.na()
}

# Return only non-NA, non-empty entrez IDs
entrez_vec <- function(ensembl_ids) {
  e <- ensembl_to_entrez(ensembl_ids)
  e[!is.na(e) & e != "NA" & nchar(e) > 0]
}

safe_gsea_go <- function(gene_list) {
  tryCatch(
    gseGO(geneList      = gene_list,
          OrgDb         = orgdb,
          ont           = "ALL",
          keyType       = "ENTREZID",
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.1,
          verbose       = FALSE),
    error = function(e) { message("gseGO failed: ", e$message); NULL }
  )
}

safe_gsea_kegg <- function(gene_list) {
  tryCatch(
    gseKEGG(geneList      = gene_list,
            organism      = kegg_org,
            pAdjustMethod = "BH",
            pvalueCutoff  = 0.1,
            verbose       = FALSE),
    error = function(e) { message("gseKEGG failed: ", e$message); NULL }
  )
}

safe_enrich_go <- function(gene_ids, universe_ids) {
  if (length(gene_ids) < 10) return(NULL)
  tryCatch(
    enrichGO(gene          = gene_ids,
             universe      = universe_ids,
             OrgDb         = orgdb,
             ont           = "ALL",
             keyType       = "ENTREZID",
             pAdjustMethod = "BH",
             pvalueCutoff  = padj_thr,
             qvalueCutoff  = 0.2,
             readable      = TRUE),
    error = function(e) { message("enrichGO failed: ", e$message); NULL }
  )
}

safe_enrich_kegg <- function(gene_ids, universe_ids) {
  if (length(gene_ids) < 10) return(NULL)
  tryCatch(
    enrichKEGG(gene          = gene_ids,
               universe      = universe_ids,
               organism      = kegg_org,
               pAdjustMethod = "BH",
               pvalueCutoff  = padj_thr),
    error = function(e) { message("enrichKEGG failed: ", e$message); NULL }
  )
}

result_to_df <- function(enrich_obj, top_n = 20) {
  if (is.null(enrich_obj)) return(NULL)
  df <- tryCatch(as.data.frame(enrich_obj), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  head(df[order(df$p.adjust), ], top_n)
}

save_plot_safe <- function(plot_obj, path, w = 12, h = 8) {
  tryCatch(
    ggsave(path, plot_obj, width = w, height = h, dpi = 200),
    error = function(e) message("Plot save failed for ", path, ": ", e$message)
  )
}

for (ct in contrasts_cfg) {
  ct_name <- ct$name
  message(sprintf("[04_enrichment] Processing contrast: %s", ct_name))

  results_file <- file.path(output_dir, paste0(ct_name, "_results.csv"))
  if (!file.exists(results_file)) {
    message("Results file not found, skipping: ", results_file)
    next
  }

  res_df <- read.csv(results_file, stringsAsFactors = FALSE)
  if (!"gene_id" %in% colnames(res_df)) res_df$gene_id <- rownames(res_df)

  # Build ranked gene list (all genes with log2FC)
  res_ranked <- res_df %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    arrange(desc(log2FoldChange))

  universe_entrez <- entrez_vec(res_df$gene_id)

  # Build ranked list: map each gene to entrez (NA-preserving), then filter
  ranked_entrez <- ensembl_to_entrez(res_ranked$gene_id)
  entrez_map <- data.frame(
    gene_id = res_ranked$gene_id,
    entrez  = ranked_entrez,
    lfc     = res_ranked$log2FoldChange,
    stringsAsFactors = FALSE
  )
  entrez_map <- entrez_map[!is.na(entrez_map$entrez) &
                             entrez_map$entrez != "NA" &
                             nchar(entrez_map$entrez) > 0 &
                             !duplicated(entrez_map$entrez), ]

  gene_list <- setNames(entrez_map$lfc, entrez_map$entrez)
  gene_list <- sort(gene_list, decreasing = TRUE)

  # Significant genes for ORA
  sig_genes    <- res_df[!is.na(res_df$padj) & res_df$padj < padj_thr &
                           abs(res_df$log2FoldChange) >= lfc_thr, ]
  sig_entrez   <- entrez_vec(sig_genes$gene_id)
  ora_skipped  <- length(sig_entrez) < 10
  if (ora_skipped) message(sprintf("[04_enrichment] <10 sig genes (%d), skipping ORA", length(sig_entrez)))

  # ── GSEA ─────────────────────────────────────────────────────────────────────
  message("[04_enrichment] Running GSEA GO...")
  gsea_go   <- if (length(gene_list) > 10) safe_gsea_go(gene_list) else NULL
  message("[04_enrichment] Running GSEA KEGG...")
  gsea_kegg <- if (length(gene_list) > 10) safe_gsea_kegg(gene_list) else NULL

  # ── ORA ──────────────────────────────────────────────────────────────────────
  message("[04_enrichment] Running ORA GO...")
  ora_go   <- if (!ora_skipped) safe_enrich_go(sig_entrez, universe_entrez) else NULL
  message("[04_enrichment] Running ORA KEGG...")
  ora_kegg <- if (!ora_skipped) safe_enrich_kegg(sig_entrez, universe_entrez) else NULL

  # ── Plots ─────────────────────────────────────────────────────────────────────
  ct_plot_dir <- file.path(enrich_plots_dir, ct_name)
  dir.create(ct_plot_dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    tryCatch({
      p <- dotplot(gsea_go, showCategory = 15, split = ".sign") +
        facet_grid(. ~ .sign) +
        labs(
          title    = paste("GSEA \u2014 Gene Ontology:", ct_name),
          subtitle = "Normalized Enrichment Score (NES), BH-adjusted p < 0.1"
        ) +
        theme_pub() +
        theme(axis.text.y = element_text(size = 8.5),
              legend.position = "right")
      save_plot_safe(p, file.path(ct_plot_dir, "gsea_go_dotplot.png"), w = 14, h = 8)
    }, error = function(e) message("GSEA GO dotplot failed: ", e$message))
  }

  if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
    tryCatch({
      p <- dotplot(gsea_kegg, showCategory = 15, split = ".sign") +
        facet_grid(. ~ .sign) +
        labs(
          title    = paste("GSEA \u2014 KEGG Pathways:", ct_name),
          subtitle = "Normalized Enrichment Score (NES), BH-adjusted p < 0.1"
        ) +
        theme_pub() +
        theme(axis.text.y = element_text(size = 8.5),
              legend.position = "right")
      save_plot_safe(p, file.path(ct_plot_dir, "gsea_kegg_dotplot.png"), w = 14, h = 8)
    }, error = function(e) message("GSEA KEGG dotplot failed: ", e$message))
  }

  if (!is.null(ora_go) && nrow(as.data.frame(ora_go)) > 0) {
    tryCatch({
      p <- dotplot(ora_go, showCategory = 15) +
        labs(
          title    = paste("ORA \u2014 Gene Ontology:", ct_name),
          subtitle = sprintf("Significant genes: %d  |  BH-adjusted p < %s", length(sig_entrez), padj_thr)
        ) +
        theme_pub() +
        theme(axis.text.y = element_text(size = 9),
              legend.position = "right")
      save_plot_safe(p, file.path(ct_plot_dir, "ora_go_dotplot.png"), w = 10, h = 7)
    }, error = function(e) message("ORA GO dotplot failed: ", e$message))
  }

  if (!is.null(ora_kegg) && nrow(as.data.frame(ora_kegg)) > 0) {
    tryCatch({
      p <- dotplot(ora_kegg, showCategory = 15) +
        labs(
          title    = paste("ORA \u2014 KEGG Pathways:", ct_name),
          subtitle = sprintf("Significant genes: %d  |  BH-adjusted p < %s", length(sig_entrez), padj_thr)
        ) +
        theme_pub() +
        theme(axis.text.y = element_text(size = 9),
              legend.position = "right")
      save_plot_safe(p, file.path(ct_plot_dir, "ora_kegg_dotplot.png"), w = 10, h = 7)
    }, error = function(e) message("ORA KEGG dotplot failed: ", e$message))
  }

  # ── Write JSON ────────────────────────────────────────────────────────────────
  enrich_result <- list(
    contrast_name = ct_name,
    ora_skipped   = ora_skipped,
    gsea_go       = result_to_df(gsea_go),
    gsea_kegg     = result_to_df(gsea_kegg),
    ora_go        = result_to_df(ora_go),
    ora_kegg      = result_to_df(ora_kegg),
    plots         = list(
      gsea_go_dotplot   = file.path(ct_plot_dir, "gsea_go_dotplot.png"),
      gsea_kegg_dotplot = file.path(ct_plot_dir, "gsea_kegg_dotplot.png"),
      ora_go_dotplot    = file.path(ct_plot_dir, "ora_go_dotplot.png"),
      ora_kegg_dotplot  = file.path(ct_plot_dir, "ora_kegg_dotplot.png")
    )
  )

  write(toJSON(enrich_result, auto_unbox = TRUE, pretty = TRUE, null = "null"),
        file.path(output_dir, paste0("enrichment_results_", ct_name, ".json")))

  message(sprintf("[04_enrichment] Contrast %s done.", ct_name))
}

message("[04_enrichment] All contrasts processed.")
