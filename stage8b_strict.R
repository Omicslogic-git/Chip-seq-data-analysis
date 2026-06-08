#!/usr/bin/env Rscript
# =============================================================================
# GSE68952 ChIP-seq Pipeline — STAGE 8b: Strict Pathway Analysis
# =============================================================================
# Re-runs GO + KEGG enrichment using STRICT gene lists (|log2FC| >= 1.0)
# instead of the loose cutoff (|log2FC| >= 0.585) used in Stage 8a.
#
# RATIONALE:
#   Stage 8a used 1,230 genes (loose FC) → GO dominated by morphogenesis terms
#   Stage 8b uses 232 lost + 71 gained genes (strict FC) → cleaner neuronal
#   enrichment expected: dopamine signalling, synaptic transmission, etc.
#
# STRICT GENE COUNTS (from Stage 8a):
#   Lost  (|FC|>=1.0, FDR<0.05): 237 peaks → 232 genes
#   Gained(|FC|>=1.0, FDR<0.05):  72 peaks →  71 genes
#
# OUTPUT: 08_annotation/strict/
#   GSE68952_strict_GO_loss.csv
#   GSE68952_strict_GO_gain.csv
#   GSE68952_strict_KEGG_loss.csv
#   GSE68952_strict_KEGG_gain.csv
#   plots/  — dotplots, barplots, cnetplots
#
# HOW TO RUN:
#   conda activate diffbind_env
#   Rscript stage8b_strict.R 2>&1 | tee run_stage8b.log
# =============================================================================

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggplot2)
  library(enrichplot)
})

cat("=============================================================\n")
cat("  GSE68952 Stage 8b: Strict Pathway Analysis (|FC| >= 1.0)\n")
cat("  Date:", format(Sys.time()), "\n")
cat("=============================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================
ANNO_DIR   <- "08_annotation"
OUT_DIR    <- file.path(ANNO_DIR, "strict")
PLOT_DIR   <- file.path(OUT_DIR, "plots")
FC_STRICT  <- 1.0
FDR_CUTOFF <- 0.05

dir.create(OUT_DIR,  showWarnings=FALSE, recursive=TRUE)
dir.create(PLOT_DIR, showWarnings=FALSE, recursive=TRUE)

# =============================================================================
# STEP 1 — Load strict annotated gene lists from Stage 8a
# =============================================================================
cat("── Step 1: Load strict annotated gene lists ──\n")

loss_file <- file.path(ANNO_DIR, "GSE68952_loss_annotated_strict.csv")
gain_file <- file.path(ANNO_DIR, "GSE68952_gain_annotated_strict.csv")
all_file  <- file.path(ANNO_DIR, "GSE68952_peaks_annotated.csv")

for (f in c(loss_file, gain_file, all_file)) {
  if (!file.exists(f)) stop("Missing file: ", f, "\nRun stage8_annotation.R first.")
}

anno_loss <- read.csv(loss_file, stringsAsFactors=FALSE)
anno_gain <- read.csv(gain_file, stringsAsFactors=FALSE)
anno_all  <- read.csv(all_file,  stringsAsFactors=FALSE)

cat(sprintf("[OK] Loss peaks (strict): %d\n", nrow(anno_loss)))
cat(sprintf("[OK] Gain peaks (strict): %d\n", nrow(anno_gain)))
cat(sprintf("[OK] All peaks (universe): %d\n", nrow(anno_all)))

# Extract gene IDs — filter to peaks within 10 kb of TSS
get_genes <- function(df, label, tss_kb=10) {
  near <- df[!is.na(df$distanceToTSS) & abs(df$distanceToTSS) <= tss_kb*1000, ]
  symbols   <- unique(na.omit(near$SYMBOL))
  entrezids <- unique(na.omit(near$geneId))
  cat(sprintf("  %s: %d peaks within %d kb TSS → %d genes\n",
              label, nrow(near), tss_kb, length(symbols)))
  list(symbols=symbols, entrezids=entrezids, peaks=near)
}

loss_genes <- get_genes(anno_loss, "Lost (strict)")
gain_genes <- get_genes(anno_gain, "Gained (strict)")
universe   <- unique(na.omit(anno_all$geneId))

cat(sprintf("\n  Gene universe: %d genes\n\n", length(universe)))

# Print gene lists
cat("  Lost genes (strict, first 30):\n  ")
cat(paste(head(loss_genes$symbols, 30), collapse=", "), "\n\n")
cat("  Gained genes (strict, first 20):\n  ")
cat(paste(head(gain_genes$symbols, 20), collapse=", "), "\n\n")

# Validate positive controls
cat("── Positive control check (strict list) ──\n")
hd_targets <- c("DRD2","BDNF","PENK","HTT","ADORA2A","PDE10A",
                 "CAMK2A","CALB1","DARPP32","PPP1R1B","GRM4","GRM2",
                 "CHRNA4","CHRNA5","GABRQ","NPY","AVP","OXTR",
                 "NFATC1","NFATC2","NFATC4")
for (g in hd_targets) {
  in_loss <- g %in% loss_genes$symbols
  in_gain <- g %in% gain_genes$symbols
  if (in_loss) cat(sprintf("  [LOSS]  %-12s reduced H3K4me3 in HD\n", g))
  if (in_gain) cat(sprintf("  [GAIN]  %-12s gained H3K4me3 in HD\n",  g))
}
cat("\n")

# =============================================================================
# STEP 2 — GO Enrichment (strict gene lists)
# =============================================================================
cat("── Step 2: GO Enrichment (strict) ──\n")

run_go_strict <- function(entrez_ids, universe, label, ont="BP") {
  cat(sprintf("  GO-%s for %s (%d genes)...\n", ont, label, length(entrez_ids)))
  if (length(entrez_ids) < 5) {
    cat(sprintf("  [SKIP] Too few genes\n"))
    return(NULL)
  }
  ego <- enrichGO(
    gene          = entrez_ids,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE,
    minGSSize     = 5,
    maxGSSize     = 500
  )
  if (is.null(ego) || nrow(ego) == 0) {
    cat(sprintf("  [NOTE] No significant terms\n"))
    return(NULL)
  }
  # Simplify to remove highly redundant GO terms
  ego_simp <- tryCatch(
    simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min),
    error = function(e) ego
  )
  cat(sprintf("  [OK] %d terms (%d after simplify)\n", nrow(ego), nrow(ego_simp)))
  ego_simp
}

go_loss_bp <- run_go_strict(loss_genes$entrezids, universe, "Lost (strict)", "BP")
go_gain_bp <- run_go_strict(gain_genes$entrezids, universe, "Gained (strict)", "BP")
go_loss_mf <- run_go_strict(loss_genes$entrezids, universe, "Lost MF (strict)", "MF")
cat("\n")

# Write results
if (!is.null(go_loss_bp)) {
  write.csv(as.data.frame(go_loss_bp),
            file.path(OUT_DIR, "GSE68952_strict_GO_loss_BP.csv"), row.names=FALSE)
  cat("[OK] GO-BP loss written\n")
}
if (!is.null(go_gain_bp)) {
  write.csv(as.data.frame(go_gain_bp),
            file.path(OUT_DIR, "GSE68952_strict_GO_gain_BP.csv"), row.names=FALSE)
  cat("[OK] GO-BP gain written\n")
}
if (!is.null(go_loss_mf)) {
  write.csv(as.data.frame(go_loss_mf),
            file.path(OUT_DIR, "GSE68952_strict_GO_loss_MF.csv"), row.names=FALSE)
  cat("[OK] GO-MF loss written\n")
}
cat("\n")

# =============================================================================
# STEP 3 — KEGG Enrichment (strict)
# =============================================================================
cat("── Step 3: KEGG Enrichment (strict) ──\n")

run_kegg_strict <- function(entrez_ids, universe, label) {
  cat(sprintf("  KEGG for %s (%d genes)...\n", label, length(entrez_ids)))
  if (length(entrez_ids) < 5) return(NULL)
  ekegg <- enrichKEGG(
    gene          = entrez_ids,
    organism      = "hsa",
    universe      = universe,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 5
  )
  if (is.null(ekegg) || nrow(ekegg) == 0) {
    cat(sprintf("  [NOTE] No significant pathways\n"))
    return(NULL)
  }
  ekegg <- setReadable(ekegg, OrgDb=org.Hs.eg.db, keyType="ENTREZID")
  cat(sprintf("  [OK] %d pathways\n", nrow(ekegg)))
  ekegg
}

kegg_loss <- run_kegg_strict(loss_genes$entrezids, universe, "Lost (strict)")
kegg_gain <- run_kegg_strict(gain_genes$entrezids, universe, "Gained (strict)")

if (!is.null(kegg_loss)) {
  write.csv(as.data.frame(kegg_loss),
            file.path(OUT_DIR, "GSE68952_strict_KEGG_loss.csv"), row.names=FALSE)
  cat("[OK] KEGG loss written\n")
}
if (!is.null(kegg_gain)) {
  write.csv(as.data.frame(kegg_gain),
            file.path(OUT_DIR, "GSE68952_strict_KEGG_gain.csv"), row.names=FALSE)
  cat("[OK] KEGG gain written\n")
}
cat("\n")

# =============================================================================
# STEP 4 — Plots
# =============================================================================
cat("── Step 4: Plots ──\n")

# GO dotplot — loss (simplified)
if (!is.null(go_loss_bp) && nrow(go_loss_bp) > 0) {
  p <- dotplot(go_loss_bp, showCategory=20, font.size=9) +
    ggtitle("GO-BP (Simplified) — Lost H3K4me3 in HD (|FC|>=1.0)")
  ggsave(file.path(PLOT_DIR, "01_GO_dotplot_loss_strict.png"),
         p, width=10, height=8, dpi=150)
  cat("[OK] GO dotplot loss (strict)\n")

  # Barplot
  p2 <- barplot(go_loss_bp, showCategory=20, font.size=9) +
    ggtitle("GO-BP — Lost H3K4me3 in HD (|FC|>=1.0)")
  ggsave(file.path(PLOT_DIR, "02_GO_barplot_loss_strict.png"),
         p2, width=10, height=8, dpi=150)
  cat("[OK] GO barplot loss (strict)\n")

  # cnetplot — gene-concept network (shows which genes drive each GO term)
  if (nrow(go_loss_bp) >= 3) {
    p3 <- tryCatch({
      cnetplot(go_loss_bp, showCategory=10, foldChange=NULL,
               node_label="category") +
        ggtitle("Gene-Concept Network — Lost H3K4me3 (|FC|>=1.0)")
    }, error = function(e) NULL)
    if (!is.null(p3)) {
      ggsave(file.path(PLOT_DIR, "03_GO_cnetplot_loss_strict.png"),
             p3, width=12, height=10, dpi=150)
      cat("[OK] GO cnetplot loss (strict)\n")
    }
  }
}

# GO dotplot — gain
if (!is.null(go_gain_bp) && nrow(go_gain_bp) > 0) {
  p4 <- dotplot(go_gain_bp, showCategory=20, font.size=9) +
    ggtitle("GO-BP (Simplified) — Gained H3K4me3 in HD (|FC|>=1.0)")
  ggsave(file.path(PLOT_DIR, "04_GO_dotplot_gain_strict.png"),
         p4, width=10, height=8, dpi=150)
  cat("[OK] GO dotplot gain (strict)\n")
}

# KEGG dotplot — loss
if (!is.null(kegg_loss) && nrow(kegg_loss) > 0) {
  p5 <- dotplot(kegg_loss, showCategory=20, font.size=9) +
    ggtitle("KEGG Pathways — Lost H3K4me3 in HD (|FC|>=1.0)")
  ggsave(file.path(PLOT_DIR, "05_KEGG_dotplot_loss_strict.png"),
         p5, width=10, height=7, dpi=150)
  cat("[OK] KEGG dotplot loss (strict)\n")

  # KEGG cnetplot
  p6 <- tryCatch({
    cnetplot(kegg_loss, showCategory=8, node_label="category") +
      ggtitle("KEGG Gene-Concept Network — Lost H3K4me3 (|FC|>=1.0)")
  }, error = function(e) NULL)
  if (!is.null(p6)) {
    ggsave(file.path(PLOT_DIR, "06_KEGG_cnetplot_loss_strict.png"),
           p6, width=12, height=10, dpi=150)
    cat("[OK] KEGG cnetplot loss (strict)\n")
  }
}

# KEGG dotplot — gain
if (!is.null(kegg_gain) && nrow(kegg_gain) > 0) {
  p7 <- dotplot(kegg_gain, showCategory=20, font.size=9) +
    ggtitle("KEGG Pathways — Gained H3K4me3 in HD (|FC|>=1.0)")
  ggsave(file.path(PLOT_DIR, "07_KEGG_dotplot_gain_strict.png"),
         p7, width=10, height=7, dpi=150)
  cat("[OK] KEGG dotplot gain (strict)\n")
}

# GO Molecular Function dotplot (loss)
if (!is.null(go_loss_mf) && nrow(go_loss_mf) > 0) {
  p8 <- dotplot(go_loss_mf, showCategory=15, font.size=9) +
    ggtitle("GO Molecular Function — Lost H3K4me3 in HD (|FC|>=1.0)")
  ggsave(file.path(PLOT_DIR, "08_GO_MF_dotplot_loss_strict.png"),
         p8, width=10, height=7, dpi=150)
  cat("[OK] GO-MF dotplot loss (strict)\n")
}
cat("\n")

# =============================================================================
# STEP 5 — Print top results to terminal
# =============================================================================
cat("── Step 5: Top pathway results (strict) ──\n\n")

if (!is.null(go_loss_bp) && nrow(go_loss_bp) > 0) {
  go_df <- as.data.frame(go_loss_bp)
  cat("  Top 15 GO-BP terms — Lost H3K4me3 (strict |FC|>=1.0):\n")
  cat(sprintf("  %-52s %8s %6s\n", "Description", "p.adj", "Count"))
  cat("  ", strrep("-", 70), "\n", sep="")
  top15 <- head(go_df[order(go_df$p.adjust), ], 15)
  for (i in seq_len(nrow(top15))) {
    cat(sprintf("  %-52s %8.2e %6s\n",
                substr(top15$Description[i], 1, 52),
                top15$p.adjust[i],
                top15$Count[i]))
  }
  cat("\n")
}

if (!is.null(kegg_loss) && nrow(kegg_loss) > 0) {
  kegg_df <- as.data.frame(kegg_loss)
  cat("  Top 10 KEGG Pathways — Lost H3K4me3 (strict |FC|>=1.0):\n")
  cat(sprintf("  %-52s %8s %6s\n", "Pathway", "p.adj", "Count"))
  cat("  ", strrep("-", 70), "\n", sep="")
  top10 <- head(kegg_df[order(kegg_df$p.adjust), ], 10)
  for (i in seq_len(nrow(top10))) {
    cat(sprintf("  %-52s %8.2e %6s\n",
                substr(top10$Description[i], 1, 52),
                top10$p.adjust[i],
                top10$Count[i]))
  }
  cat("\n")
}

# =============================================================================
# STEP 6 — Comparison: loose vs strict
# =============================================================================
cat("── Step 6: Loose vs Strict comparison ──\n\n")

loose_go_file <- file.path(ANNO_DIR, "GSE68952_GO_loss.csv")
if (file.exists(loose_go_file)) {
  loose_go <- read.csv(loose_go_file)
  strict_go <- if (!is.null(go_loss_bp)) as.data.frame(go_loss_bp) else data.frame()

  cat(sprintf("  GO-BP terms (loose |FC|>=0.585) : %d\n", nrow(loose_go)))
  cat(sprintf("  GO-BP terms (strict |FC|>=1.0)  : %d\n", nrow(strict_go)))

  if (nrow(strict_go) > 0) {
    # Terms shared between both lists
    shared <- intersect(loose_go$Description, strict_go$Description)
    cat(sprintf("  Shared terms                    : %d\n", length(shared)))
    cat(sprintf("  Strict-only terms               : %d (higher confidence)\n\n",
                nrow(strict_go) - length(shared)))

    cat("  Top 5 strict-only GO terms (not in loose list):\n")
    strict_only <- strict_go[!strict_go$Description %in% loose_go$Description, ]
    strict_only <- strict_only[order(strict_only$p.adjust), ]
    for (i in seq_len(min(5, nrow(strict_only)))) {
      cat(sprintf("    %s (p.adj=%.2e)\n",
                  strict_only$Description[i], strict_only$p.adjust[i]))
    }
  }
}
cat("\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("=============================================================\n")
cat("  STAGE 8b SUMMARY\n")
cat("=============================================================\n\n")
cat(sprintf("  Strict gene lists used:\n"))
cat(sprintf("    Lost  (|FC|>=1.0): %d genes\n", length(loss_genes$symbols)))
cat(sprintf("    Gained(|FC|>=1.0): %d genes\n", length(gain_genes$symbols)))
cat(sprintf("  GO-BP terms  (loss) : %s\n",
    if (!is.null(go_loss_bp)) nrow(go_loss_bp) else "0 (no enrichment)"))
cat(sprintf("  GO-BP terms  (gain) : %s\n",
    if (!is.null(go_gain_bp)) nrow(go_gain_bp) else "0 (no enrichment)"))
cat(sprintf("  KEGG pathways(loss) : %s\n",
    if (!is.null(kegg_loss)) nrow(kegg_loss) else "0 (no enrichment)"))
cat(sprintf("  KEGG pathways(gain) : %s\n",
    if (!is.null(kegg_gain)) nrow(kegg_gain) else "0 (no enrichment)"))
cat(sprintf("\n  Output: %s/\n", OUT_DIR))
cat(sprintf("  Plots : %s/\n\n", PLOT_DIR))
cat("  Next step: bash stage9_homer.sh\n")
cat(sprintf("\n  %s: Stage 8b finished.\n", format(Sys.time())))
