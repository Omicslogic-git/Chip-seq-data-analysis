#!/usr/bin/env Rscript
# =============================================================================
# GSE68952 ChIP-seq Pipeline — STAGE 8: Annotation & Pathway Analysis
# =============================================================================
# Study    : H3K4me3 ChIP-seq · Huntington Disease vs Control · hg38
#
# STAGE 7 RESULTS (baked in):
#   Total differential (FDR<0.05)     : 8,203
#   Lost in HD                        : 4,063 (49.5%)
#   Gained in HD                      : 4,140 (50.5%)
#   PCA outliers                      : C4, HD2, HD4 (low-read samples)
#   BDNF                              : present, not differential (FDR=0.749)
#   DRD2/PENK                         : not in consensus — check full peak set
#
# ANALYSIS PLAN:
#   Primary set   : FDR < 0.05 + |log2FC| >= 1.0  (~stricter, cleaner list)
#   Secondary set : FDR < 0.05 only               (all 8,203 for broad pathways)
#   Annotation    : ChIPseeker → nearest TSS, genomic feature distribution
#   Pathways      : clusterProfiler GO + KEGG on lost peaks (HD biology focus)
#   Motifs        : Prepare BED files for Homer (run in Stage 9)
#
# OUTPUT:
#   08_annotation/
#     GSE68952_peaks_annotated.csv          — all differential peaks annotated
#     GSE68952_loss_annotated.csv           — lost peaks with gene names
#     GSE68952_gain_annotated.csv           — gained peaks with gene names
#     GSE68952_loss_genes.txt               — gene list for lost peaks
#     GSE68952_gain_genes.txt               — gene list for gained peaks
#     GSE68952_GO_loss.csv                  — GO enrichment, lost peaks
#     GSE68952_GO_gain.csv                  — GO enrichment, gained peaks
#     GSE68952_KEGG_loss.csv                — KEGG enrichment, lost peaks
#     plots/                               — all annotation + pathway plots
#
# REQUIREMENTS (diffbind_env or chipseq_hg38):
#   ChIPseeker   TxDb.Hsapiens.UCSC.hg38.knownGene   org.Hs.eg.db
#   clusterProfiler   ggplot2   pheatmap   ReactomePA (optional)
#
# HOW TO RUN:
#   conda activate diffbind_env
#   Rscript stage8_annotation.R 2>&1 | tee run_stage8.log
# =============================================================================

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
})

cat("=============================================================\n")
cat("  GSE68952 Stage 8: Annotation & Pathway Analysis\n")
cat("  Date:", format(Sys.time()), "\n")
cat("=============================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================
DIFFBIND_DIR   <- "07_diffbind"
OUT_DIR        <- "08_annotation"
PLOT_DIR       <- file.path(OUT_DIR, "plots")
HOMER_DIR      <- "09_homer"         # BED files prepared here for Stage 9

# Thresholds
FDR_CUTOFF     <- 0.05
FC_STRICT      <- 1.0    # |log2FC| >= 1.0 for primary gene list
FC_LOOSE       <- 0.585  # |log2FC| >= 0.585 (log2(1.5)) for secondary

# Annotation database
TXdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

dir.create(OUT_DIR,   showWarnings=FALSE, recursive=TRUE)
dir.create(PLOT_DIR,  showWarnings=FALSE, recursive=TRUE)
dir.create(HOMER_DIR, showWarnings=FALSE, recursive=TRUE)

# =============================================================================
# STEP 1 — Load DiffBind results
# =============================================================================
cat("── Step 1: Load DiffBind results ──\n")

results_file <- file.path(DIFFBIND_DIR, "GSE68952_diffbind_all_peaks.csv")
if (!file.exists(results_file)) {
  stop("DiffBind results not found: ", results_file,
       "\nRun stage7_diffbind.R first.")
}

res <- read.csv(results_file, stringsAsFactors=FALSE)
cat(sprintf("[OK] Loaded %d peaks from %s\n", nrow(res), results_file))

# Standardise column names (DiffBind output varies by version)
# Expected: seqnames/chr, start, end, Fold, FDR, p.value or pvalue
colnames(res) <- tolower(colnames(res))
if ("seqnames" %in% colnames(res)) res$chr <- res$seqnames
if (!"chr" %in% colnames(res))     res$chr <- res[,1]
if (!"fold" %in% colnames(res) && "log2foldchange" %in% colnames(res))
  res$fold <- res$log2foldchange

cat(sprintf("  Columns: %s\n", paste(colnames(res)[1:10], collapse=", ")))
cat(sprintf("  Total peaks: %d\n", nrow(res)))
cat(sprintf("  FDR < 0.05: %d\n", sum(res$fdr < FDR_CUTOFF, na.rm=TRUE)))
cat(sprintf("  FDR < 0.05 + |FC| >= 1.0: %d\n",
    sum(res$fdr < FDR_CUTOFF & abs(res$fold) >= FC_STRICT, na.rm=TRUE)))
cat("\n")

# =============================================================================
# STEP 2 — Convert to GRanges and annotate with ChIPseeker
# annotatePeak() assigns each peak to:
#   - Nearest TSS (and distance)
#   - Genomic feature: Promoter, 5'UTR, Exon, Intron, 3'UTR, Downstream, Intergenic
# tssRegion=c(-3000, 3000): define promoter as TSS ± 3 kb
# =============================================================================
cat("── Step 2: Annotate peaks with ChIPseeker ──\n")
cat("  Using TxDb.Hsapiens.UCSC.hg38.knownGene\n")
cat("  Promoter defined as TSS ± 3 kb\n\n")

# Convert all differential peaks to GRanges
res_sig  <- res[res$fdr < FDR_CUTOFF, ]
res_gr   <- GRanges(
  seqnames = res_sig$chr,
  ranges   = IRanges(start=res_sig$start, end=res_sig$end),
  Fold     = res_sig$fold,
  FDR      = res_sig$fdr
)

# Annotate
peak_anno <- annotatePeak(
  res_gr,
  tssRegion        = c(-3000, 3000),
  TxDb             = TXdb,
  annoDb           = "org.Hs.eg.db",
  addFlankGeneInfo = TRUE,
  verbose          = FALSE
)

anno_df <- as.data.frame(peak_anno)
cat(sprintf("[OK] Annotated %d peaks\n", nrow(anno_df)))

# Separate into loss and gain sets (using strict FC cutoff)
anno_loss <- anno_df[anno_df$Fold < -FC_STRICT, ]
anno_gain <- anno_df[anno_df$Fold >  FC_STRICT, ]
anno_loss_loose <- anno_df[anno_df$Fold < -FC_LOOSE, ]
anno_gain_loose <- anno_df[anno_df$Fold >  FC_LOOSE, ]

cat(sprintf("  Lost (|FC|>=1.0):   %d peaks\n", nrow(anno_loss)))
cat(sprintf("  Gained (|FC|>=1.0): %d peaks\n", nrow(anno_gain)))
cat(sprintf("  Lost (|FC|>=0.585): %d peaks\n", nrow(anno_loss_loose)))
cat(sprintf("  Gained (|FC|>=0.585): %d peaks\n\n", nrow(anno_gain_loose)))

# Write annotated CSVs
write.csv(anno_df,        file.path(OUT_DIR, "GSE68952_peaks_annotated.csv"),      row.names=FALSE)
write.csv(anno_loss,      file.path(OUT_DIR, "GSE68952_loss_annotated_strict.csv"),row.names=FALSE)
write.csv(anno_gain,      file.path(OUT_DIR, "GSE68952_gain_annotated_strict.csv"),row.names=FALSE)
write.csv(anno_loss_loose,file.path(OUT_DIR, "GSE68952_loss_annotated.csv"),       row.names=FALSE)
write.csv(anno_gain_loose,file.path(OUT_DIR, "GSE68952_gain_annotated.csv"),       row.names=FALSE)
cat("[OK] Annotated CSVs written\n\n")

# =============================================================================
# STEP 3 — ChIPseeker annotation plots
# =============================================================================
cat("── Step 3: ChIPseeker annotation plots ──\n")

# Annotate separate GRanges for loss and gain (loose cutoff for plots)
loss_gr <- GRanges(
  seqnames = anno_loss_loose$seqnames,
  ranges   = IRanges(start=anno_loss_loose$start, end=anno_loss_loose$end)
)
gain_gr <- GRanges(
  seqnames = anno_gain_loose$seqnames,
  ranges   = IRanges(start=anno_gain_loose$start, end=anno_gain_loose$end)
)

loss_anno <- annotatePeak(loss_gr, tssRegion=c(-3000,3000),
                          TxDb=TXdb, verbose=FALSE)
gain_anno <- annotatePeak(gain_gr, tssRegion=c(-3000,3000),
                          TxDb=TXdb, verbose=FALSE)

peak_list <- list("Lost in HD"   = loss_anno,
                  "Gained in HD" = gain_anno)

# 1. Feature distribution bar chart
png(file.path(PLOT_DIR, "01_feature_distribution.png"),
    width=1000, height=600, res=120)
plotAnnoBar(peak_list,
            title="H3K4me3 Genomic Feature Distribution — HD vs Control")
dev.off()
cat("[OK] Feature distribution: 01_feature_distribution.png\n")

# 2. TSS distance pie charts
png(file.path(PLOT_DIR, "02_TSS_distance_loss.png"),
    width=800, height=750, res=120)
par(oma=c(0,0,2,0))  # outer margin at top for title
plotAnnoPie(loss_anno)
mtext("Lost in HD — Genomic Features", outer=TRUE, cex=1.2, font=2)
dev.off()

png(file.path(PLOT_DIR, "03_TSS_distance_gain.png"),
    width=800, height=750, res=120)
par(oma=c(0,0,2,0))
plotAnnoPie(gain_anno)
mtext("Gained in HD — Genomic Features", outer=TRUE, cex=1.2, font=2)
dev.off()
cat("[OK] TSS distance pie charts written\n")

# 3. TSS distance distribution
png(file.path(PLOT_DIR, "04_TSS_distance_histogram.png"),
    width=900, height=500, res=120)
plotDistToTSS(peak_list,
              title="Distance to Nearest TSS — H3K4me3 Differential Peaks")
dev.off()
cat("[OK] TSS distance histogram: 04_TSS_distance_histogram.png\n\n")

# =============================================================================
# STEP 4 — Extract gene lists
# Use SYMBOL for readable names, ENTREZID for pathway enrichment
# distanceToTSS filter: keep peaks within 10 kb of TSS for promoter-focused
# analysis (H3K4me3 is a promoter mark — distal peaks are less reliable)
# =============================================================================
cat("── Step 4: Extract gene lists ──\n")

extract_genes <- function(anno_df, label, tss_cutoff=10000) {
  # Filter by TSS distance
  near_tss <- anno_df[abs(anno_df$distanceToTSS) <= tss_cutoff, ]
  cat(sprintf("  %s: %d peaks within %d kb of TSS\n",
              label, nrow(near_tss), tss_cutoff/1000))

  # Get unique gene symbols and Entrez IDs
  symbols  <- unique(na.omit(near_tss$SYMBOL))
  entrezids <- unique(na.omit(near_tss$geneId))

  list(peaks=near_tss, symbols=symbols, entrezids=entrezids)
}

loss_genes <- extract_genes(anno_loss_loose, "Lost in HD")
gain_genes <- extract_genes(anno_gain_loose, "Gained in HD")

# Write gene lists
writeLines(loss_genes$symbols,
           file.path(OUT_DIR, "GSE68952_loss_genes.txt"))
writeLines(gain_genes$symbols,
           file.path(OUT_DIR, "GSE68952_gain_genes.txt"))

cat(sprintf("\n  Unique genes — Lost:   %d\n", length(loss_genes$symbols)))
cat(sprintf("  Unique genes — Gained: %d\n\n",  length(gain_genes$symbols)))

# Check published positive controls
positive_controls <- c("BDNF","DRD2","PENK","HTT","DARPP32","PPP1R1B",
                        "ADORA2A","PDE10A","CNR1","DRDP32")
cat("  Positive control genes (published HD H3K4me3 targets):\n")
for (g in positive_controls) {
  in_loss  <- g %in% loss_genes$symbols
  in_gain  <- g %in% gain_genes$symbols
  in_any   <- g %in% c(loss_genes$symbols, gain_genes$symbols)
  if (in_loss)  cat(sprintf("    [LOSS]  %s — reduced H3K4me3 in HD (expected)\n", g))
  if (in_gain)  cat(sprintf("    [GAIN]  %s — gained H3K4me3 in HD\n", g))
  if (!in_any)  cat(sprintf("    [--]    %s — not in differential list\n", g))
}
cat("\n")

# =============================================================================
# STEP 5 — GO Enrichment Analysis
# Biological Process ontology — focus on lost peaks (H3K4me3 loss = gene
# silencing in HD = neuronal/synaptic pathway downregulation expected)
# universe = all genes in consensus peak set (not just differential)
# =============================================================================
cat("── Step 5: GO Enrichment Analysis ──\n")

# Background gene universe — all annotated peaks
universe_ids <- unique(na.omit(as.data.frame(peak_anno)$geneId))
cat(sprintf("  Gene universe: %d genes\n", length(universe_ids)))

run_go <- function(entrez_ids, label, ont="BP") {
  cat(sprintf("  Running GO-%s for %s (%d genes)...\n",
              ont, label, length(entrez_ids)))
  if (length(entrez_ids) < 10) {
    cat(sprintf("  [SKIP] Too few genes for GO analysis\n"))
    return(NULL)
  }
  ego <- enrichGO(
    gene          = entrez_ids,
    universe      = universe_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(ego) == 0) {
    cat(sprintf("  [NOTE] No significant GO terms for %s\n", label))
    return(NULL)
  }
  cat(sprintf("  [OK] %d significant GO terms\n", nrow(ego)))
  ego
}

go_loss <- run_go(loss_genes$entrezids, "Lost in HD",   ont="BP")
go_gain <- run_go(gain_genes$entrezids, "Gained in HD", ont="BP")

# Write GO results
if (!is.null(go_loss)) {
  write.csv(as.data.frame(go_loss),
            file.path(OUT_DIR, "GSE68952_GO_loss.csv"), row.names=FALSE)
  cat("[OK] GO loss results written\n")
}
if (!is.null(go_gain)) {
  write.csv(as.data.frame(go_gain),
            file.path(OUT_DIR, "GSE68952_GO_gain.csv"), row.names=FALSE)
  cat("[OK] GO gain results written\n")
}
cat("\n")

# =============================================================================
# STEP 6 — KEGG Pathway Analysis
# =============================================================================
cat("── Step 6: KEGG Pathway Analysis ──\n")

run_kegg <- function(entrez_ids, label) {
  cat(sprintf("  Running KEGG for %s (%d genes)...\n",
              label, length(entrez_ids)))
  if (length(entrez_ids) < 10) return(NULL)

  ekegg <- enrichKEGG(
    gene          = entrez_ids,
    organism      = "hsa",
    universe      = universe_ids,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  if (is.null(ekegg) || nrow(ekegg) == 0) {
    cat(sprintf("  [NOTE] No significant KEGG pathways for %s\n", label))
    return(NULL)
  }
  # Convert Entrez IDs to gene symbols for readability
  ekegg <- setReadable(ekegg, OrgDb=org.Hs.eg.db, keyType="ENTREZID")
  cat(sprintf("  [OK] %d significant KEGG pathways\n", nrow(ekegg)))
  ekegg
}

kegg_loss <- run_kegg(loss_genes$entrezids, "Lost in HD")
kegg_gain <- run_kegg(gain_genes$entrezids, "Gained in HD")

if (!is.null(kegg_loss)) {
  write.csv(as.data.frame(kegg_loss),
            file.path(OUT_DIR, "GSE68952_KEGG_loss.csv"), row.names=FALSE)
  cat("[OK] KEGG loss results written\n")
}
if (!is.null(kegg_gain)) {
  write.csv(as.data.frame(kegg_gain),
            file.path(OUT_DIR, "GSE68952_KEGG_gain.csv"), row.names=FALSE)
  cat("[OK] KEGG gain results written\n")
}
cat("\n")

# =============================================================================
# STEP 7 — Pathway visualisation plots
# =============================================================================
cat("── Step 7: Pathway plots ──\n")

# GO dotplot — Lost peaks
if (!is.null(go_loss) && nrow(go_loss) > 0) {
  p1 <- dotplot(go_loss, showCategory=20, font.size=9) +
    ggtitle("GO Biological Process — Lost H3K4me3 in HD") +
    theme(plot.title=element_text(size=11))
  ggsave(file.path(PLOT_DIR, "05_GO_dotplot_loss.png"),
         p1, width=10, height=8, dpi=150)
  cat("[OK] GO dotplot (loss): 05_GO_dotplot_loss.png\n")

  # GO barplot
  p2 <- barplot(go_loss, showCategory=20, font.size=9) +
    ggtitle("GO Biological Process — Lost H3K4me3 in HD")
  ggsave(file.path(PLOT_DIR, "06_GO_barplot_loss.png"),
         p2, width=10, height=8, dpi=150)
  cat("[OK] GO barplot (loss): 06_GO_barplot_loss.png\n")
}

# GO dotplot — Gained peaks
if (!is.null(go_gain) && nrow(go_gain) > 0) {
  p3 <- dotplot(go_gain, showCategory=20, font.size=9) +
    ggtitle("GO Biological Process — Gained H3K4me3 in HD") +
    theme(plot.title=element_text(size=11))
  ggsave(file.path(PLOT_DIR, "07_GO_dotplot_gain.png"),
         p3, width=10, height=8, dpi=150)
  cat("[OK] GO dotplot (gain): 07_GO_dotplot_gain.png\n")
}

# KEGG dotplot — Lost peaks
if (!is.null(kegg_loss) && nrow(kegg_loss) > 0) {
  p4 <- dotplot(kegg_loss, showCategory=20, font.size=9) +
    ggtitle("KEGG Pathways — Lost H3K4me3 in HD")
  ggsave(file.path(PLOT_DIR, "08_KEGG_dotplot_loss.png"),
         p4, width=10, height=7, dpi=150)
  cat("[OK] KEGG dotplot (loss): 08_KEGG_dotplot_loss.png\n")
}

# KEGG dotplot — Gained peaks
if (!is.null(kegg_gain) && nrow(kegg_gain) > 0) {
  p5 <- dotplot(kegg_gain, showCategory=20, font.size=9) +
    ggtitle("KEGG Pathways — Gained H3K4me3 in HD")
  ggsave(file.path(PLOT_DIR, "09_KEGG_dotplot_gain.png"),
         p5, width=10, height=7, dpi=150)
  cat("[OK] KEGG dotplot (gain): 09_KEGG_dotplot_gain.png\n")
}

# Upset plot — gene overlap across top GO terms
# Requires ggupset package
if (!requireNamespace("ggupset", quietly=TRUE)) {
  message("Installing ggupset...")
  install.packages("ggupset", repos="https://cran.r-project.org", quiet=TRUE)
}

if (!is.null(go_loss) && nrow(go_loss) >= 2 &&
    requireNamespace("ggupset", quietly=TRUE)) {
  p_upset <- tryCatch(
    enrichplot::upsetplot(go_loss, n=10),
    error = function(e) NULL
  )
  if (!is.null(p_upset)) {
    ggsave(file.path(PLOT_DIR, "10_upsetplot_loss.png"),
           p_upset, width=12, height=7, dpi=120)
    cat("[OK] Upset plot: 10_upsetplot_loss.png\n\n")
  } else {
    cat("[SKIP] Upset plot could not be generated\n\n")
  }
} else {
  cat("[SKIP] Upset plot skipped (ggupset not available or < 2 GO terms)\n\n")
}

# =============================================================================
# STEP 8 — Prepare BED files for Homer motif analysis (Stage 9)
# Homer needs: chr start end name score strand (6-column BED)
# Separate files for: all diff peaks, loss only, gain only
# Also prepare a background set (non-differential peaks from consensus)
# =============================================================================
cat("── Step 8: Prepare Homer BED files ──\n")

make_homer_bed <- function(df, filename, label) {
  # Robust column access — ChIPseeker output column names vary by version
  get_col <- function(df, candidates) {
    for (nm in candidates) {
      if (nm %in% colnames(df) && length(df[[nm]]) == nrow(df))
        return(df[[nm]])
    }
    stop(paste("Could not find column from:", paste(candidates, collapse=",")))
  }

  chr_col   <- get_col(df, c("seqnames","chr","chrom","Chr"))
  start_col <- get_col(df, c("start","Start","chromStart"))
  end_col   <- get_col(df, c("end","End","chromEnd"))
  fdr_col   <- get_col(df, c("FDR","fdr","padj","p.adjust","qvalue"))

  bed <- data.frame(
    chr    = as.character(chr_col),
    start  = as.integer(start_col) - 1L,  # Homer uses 0-based coordinates
    end    = as.integer(end_col),
    name   = paste0(label, "_", seq_len(nrow(df))),
    score  = round(-log10(pmax(as.numeric(fdr_col), 1e-300)), 2),
    strand = "."
  )
  write.table(bed, filename, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  cat(sprintf("  [OK] %s: %d peaks -> %s\n", label, nrow(bed), filename))
}

# All significant differential peaks
make_homer_bed(anno_df,
               file.path(HOMER_DIR, "all_differential.bed"),
               "diff")

# Lost peaks (strict FC)
if (nrow(anno_loss) > 0)
  make_homer_bed(anno_loss,
                 file.path(HOMER_DIR, "loss_strict.bed"),
                 "loss")

# Gained peaks (strict FC)
if (nrow(anno_gain) > 0)
  make_homer_bed(anno_gain,
                 file.path(HOMER_DIR, "gain_strict.bed"),
                 "gain")

# Lost peaks near TSS (promoter-focused — best for H3K4me3 motif analysis)
loss_tss <- anno_loss[abs(anno_loss$distanceToTSS) <= 3000, ]
if (nrow(loss_tss) > 0)
  make_homer_bed(loss_tss,
                 file.path(HOMER_DIR, "loss_promoter.bed"),
                 "loss_tss")

gain_tss <- anno_gain[abs(anno_gain$distanceToTSS) <= 3000, ]
if (nrow(gain_tss) > 0)
  make_homer_bed(gain_tss,
                 file.path(HOMER_DIR, "gain_promoter.bed"),
                 "gain_tss")

# Background: non-differential consensus peaks (FDR > 0.5)
bg <- res[res$fdr > 0.5, ]
if (nrow(bg) > 5000) bg <- bg[sample(nrow(bg), 5000), ]  # cap at 5000
# Robust background BED construction
bg_chr   <- if ("chr"      %in% colnames(bg)) bg$chr      else
            if ("seqnames" %in% colnames(bg)) as.character(bg$seqnames) else bg[,1]
bg_start <- if ("start"    %in% colnames(bg)) bg$start    else bg[,2]
bg_end   <- if ("end"      %in% colnames(bg)) bg$end      else bg[,3]
bg_bed <- data.frame(
  chr=as.character(bg_chr), start=as.integer(bg_start)-1L,
  end=as.integer(bg_end),
  name=paste0("bg_", seq_len(nrow(bg))), score=0, strand="."
)
write.table(bg_bed,
            file.path(HOMER_DIR, "background.bed"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat(sprintf("  [OK] background: %d peaks → %s/background.bed\n\n",
            nrow(bg_bed), HOMER_DIR))

# =============================================================================
# STEP 9 — Print top results summary
# =============================================================================
cat("── Step 9: Top results summary ──\n\n")

# Top 20 GO terms (lost peaks)
if (!is.null(go_loss) && nrow(go_loss) > 0) {
  go_df <- as.data.frame(go_loss)
  cat("  Top 15 GO Biological Process terms — Lost H3K4me3 in HD:\n")
  cat(sprintf("  %-50s %8s %8s\n", "Description", "p.adj", "Count"))
  cat(paste0("  ", strrep("-", 70)), "\n")
  top15 <- head(go_df[order(go_df$p.adjust), ], 15)
  for (i in seq_len(nrow(top15))) {
    cat(sprintf("  %-50s %8.2e %8s\n",
                substr(top15$Description[i], 1, 50),
                top15$p.adjust[i],
                top15$Count[i]))
  }
  cat("\n")
}

# Top 10 KEGG pathways (lost peaks)
if (!is.null(kegg_loss) && nrow(kegg_loss) > 0) {
  kegg_df <- as.data.frame(kegg_loss)
  cat("  Top 10 KEGG Pathways — Lost H3K4me3 in HD:\n")
  cat(sprintf("  %-50s %8s %8s\n", "Pathway", "p.adj", "Count"))
  cat(paste0("  ", strrep("-", 70)), "\n")
  top10 <- head(kegg_df[order(kegg_df$p.adjust), ], 10)
  for (i in seq_len(nrow(top10))) {
    cat(sprintf("  %-50s %8.2e %8s\n",
                substr(top10$Description[i], 1, 50),
                top10$p.adjust[i],
                top10$Count[i]))
  }
  cat("\n")
}

# =============================================================================
# STEP 10 — Final summary
# =============================================================================
cat("=============================================================\n")
cat("  STAGE 8 SUMMARY\n")
cat("=============================================================\n\n")
cat(sprintf("  Peaks annotated          : %d\n", nrow(anno_df)))
cat(sprintf("  Lost peaks (|FC|>=1.0)   : %d → %d genes\n",
            nrow(anno_loss), length(unique(na.omit(anno_loss$SYMBOL)))))
cat(sprintf("  Gained peaks (|FC|>=1.0) : %d → %d genes\n",
            nrow(anno_gain), length(unique(na.omit(anno_gain$SYMBOL)))))
cat(sprintf("  GO terms (loss)          : %s\n",
            if (!is.null(go_loss)) nrow(go_loss) else "0"))
cat(sprintf("  KEGG pathways (loss)     : %s\n",
            if (!is.null(kegg_loss)) nrow(kegg_loss) else "0"))
cat("\n  Output files:\n")
cat(sprintf("    %s/GSE68952_peaks_annotated.csv\n", OUT_DIR))
cat(sprintf("    %s/GSE68952_loss_genes.txt\n",       OUT_DIR))
cat(sprintf("    %s/GSE68952_gain_genes.txt\n",       OUT_DIR))
cat(sprintf("    %s/GSE68952_GO_loss.csv\n",          OUT_DIR))
cat(sprintf("    %s/GSE68952_KEGG_loss.csv\n",        OUT_DIR))
cat(sprintf("    %s/plots/  (%d plots)\n",            OUT_DIR, 10))
cat(sprintf("    %s/        (Homer BED files)\n\n",   HOMER_DIR))
cat("  Next step: bash stage9_homer.sh\n")
cat(sprintf("\n  %s: Stage 8 finished.\n", format(Sys.time())))
