# =============================================================================
# GSE68952 ChIP-seq Pipeline — STAGE 7: Differential Binding (DiffBind/DESeq2)
# =============================================================================
# Study    : H3K4me3 ChIP-seq · Huntington Disease vs Control · hg38
#
# STAGE 6 RESULTS (baked in):
#   Per-sample peaks    : 22,772 – 27,621 (all 12 samples)
#   Control avg peaks   : 25,090
#   HD avg peaks        : 23,563  (6.5% fewer — global H3K4me3 loss visible)
#   Consensus peaks     : 25,217  (present in >= 2 samples, blacklist-filtered)
#   Paper expectation   : ~28,608 total, ~2,830 differential (FDR < 0.05)
#                         ~55% down in HD
#
# ANALYSIS DESIGN:
#   Comparison  : HD (n=6) vs Control (n=6)
#   Peak set    : 06_peaks/consensus/GSE68952_consensus_peaks.bed
#   BAMs        : 04_aligned/SRR*.markdup.bam (duplicates marked, not removed)
#   Input       : shared SRR2028284 for all 12 ChIP samples
#   Model       : DESeq2 (recommended for ChIP-seq differential binding)
#   Normalise   : RLE (DESeq2 default — handles library size differences)
#   FDR cutoff  : 0.05
#   FC cutoff   : 1.5 (log2FC > 0.585 or < -0.585)
#
# OUTPUT FILES:
#   07_diffbind/GSE68952_diffbind_results.csv     — all differential peaks
#   07_diffbind/GSE68952_diffbind_gain.bed        — gained in HD (up)
#   07_diffbind/GSE68952_diffbind_loss.bed        — lost in HD (down)
#   07_diffbind/GSE68952_diffbind_all_counts.csv  — raw counts all 25,217 peaks
#   07_diffbind/plots/                            — QC and result plots
#
# REQUIREMENTS (chipseq_hg38 conda environment):
#   R >= 4.3   DiffBind >= 3.12   DESeq2 >= 1.42   ggplot2   pheatmap
#
# HOW TO RUN:
#   conda activate chipseq_hg38
#   Rscript stage7_diffbind.R 2>&1 | tee run_stage7.log
# =============================================================================

suppressPackageStartupMessages({
  library(DiffBind)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

cat("=============================================================\n")
cat("  GSE68952 Stage 7: Differential Binding Analysis\n")
cat("  Date:", format(Sys.time()), "\n")
cat("=============================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================
BAM_DIR        <- "04_aligned"
PEAK_DIR       <- "06_peaks"
CONSENSUS_BED  <- file.path(PEAK_DIR, "consensus", "GSE68952_consensus_peaks.bed")
OUT_DIR        <- "07_diffbind"
PLOT_DIR       <- file.path(OUT_DIR, "plots")

FDR_CUTOFF     <- 0.05
FC_CUTOFF      <- log2(1.5)   # log2(1.5) = 0.585

dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STEP 1 — Build DiffBind sample sheet
# DiffBind requires: SampleID, Tissue, Factor, Condition, Replicate,
#                   bamReads, ControlID, bamControl, Peaks, PeakCaller
# Shared Input: same BAM listed for all 12 ChIP samples
# =============================================================================
cat("── Step 1: Build sample sheet ──\n")

INPUT_BAM <- file.path(BAM_DIR, "SRR2028284.markdup.bam")

samples <- data.frame(
  SampleID  = c("C1","C2","C3","C4","C5","C6",
                "HD1","HD2","HD3","HD4","HD5","HD6"),
  Tissue    = "PFC_neuron",
  Factor    = "H3K4me3",
  Condition = c(rep("Control",6), rep("HD",6)),
  Replicate = c(1:6, 1:6),
  bamReads  = file.path(BAM_DIR, c(
    "SRR2028278.markdup.bam","SRR2028279.markdup.bam",
    "SRR2028280.markdup.bam","SRR2028281.markdup.bam",
    "SRR2028282.markdup.bam","SRR2028283.markdup.bam",
    "SRR2028285.markdup.bam","SRR2028286.markdup.bam",
    "SRR2028287.markdup.bam","SRR2028288.markdup.bam",
    "SRR2028289.markdup.bam","SRR2028290.markdup.bam"
  )),
  ControlID  = "Input",
  bamControl = INPUT_BAM,
  Peaks      = c(
    file.path(PEAK_DIR,"SRR2028278_C1","C1_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028279_C2","C2_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028280_C3","C3_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028281_C4","C4_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028282_C5","C5_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028283_C6","C6_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028285_HD1","HD1_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028286_HD2","HD2_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028287_HD3","HD3_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028288_HD4","HD4_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028289_HD5","HD5_peaks.broadPeak"),
    file.path(PEAK_DIR,"SRR2028290_HD6","HD6_peaks.broadPeak")
  ),
  PeakCaller = "macs",
  stringsAsFactors = FALSE
)

# Verify all files exist before proceeding
cat("Verifying input files...\n")
missing <- c()
for (f in c(samples$bamReads, INPUT_BAM, samples$Peaks)) {
  if (!file.exists(f)) missing <- c(missing, f)
}
if (length(missing) > 0) {
  cat("MISSING FILES:\n")
  cat(paste(" ", missing, collapse="\n"), "\n")
  stop("Fix missing files before proceeding.")
} else {
  cat("[OK] All BAMs and peak files present\n")
}

# Write sample sheet to disk for reference
write.csv(samples, file.path(OUT_DIR, "GSE68952_diffbind_samplesheet.csv"),
          row.names = FALSE, quote = FALSE)
cat("[OK] Sample sheet written:", file.path(OUT_DIR, "GSE68952_diffbind_samplesheet.csv"), "\n\n")

# =============================================================================
# STEP 2 — Load data into DiffBind
# dba() reads the sample sheet and peak files
# peakFormat="MACS"  : broadPeak format from MACS3
# consensus          : use our pre-built consensus peak set
# =============================================================================
cat("── Step 2: Load into DiffBind ──\n")

# Load with per-sample peaks first (for QC plots)
dba_obj <- dba(
  sampleSheet = samples,
  peakFormat  = "MACS",
  scoreCol    = 5,          # use score column from broadPeak
  bRemoveM    = TRUE,       # remove chrM peaks
  bRemoveRandom = TRUE      # remove random/unplaced contigs
)

cat("\nDiffBind object summary:\n")
print(dba_obj)

# =============================================================================
# STEP 3 — DiffBind QC plots (before counting)
# Correlation heatmap and PCA using peak overlap scores
# This confirms sample grouping before we invest in read counting
# =============================================================================
cat("\n── Step 3: Pre-count QC plots ──\n")

# Correlation heatmap (peak overlap)
png(file.path(PLOT_DIR, "01_correlation_heatmap_peaks.png"),
    width=800, height=700, res=100)
dba.plotHeatmap(dba_obj,
                main="GSE68952 — Peak Overlap Correlation",
                ColAttributes=DBA_CONDITION)
dev.off()
cat("[OK] Correlation heatmap (peaks): 01_correlation_heatmap_peaks.png\n")

# PCA (peak overlap scores)
png(file.path(PLOT_DIR, "02_PCA_peaks.png"), width=700, height=600, res=100)
dba.plotPCA(dba_obj,
            attributes = DBA_CONDITION,
            label      = DBA_ID)
dev.off()
cat("[OK] PCA (peaks): 02_PCA_peaks.png\n\n")

# =============================================================================
# STEP 4 — Count reads in consensus peaks
# dba.count() re-counts reads in ALL consensus regions for ALL samples.
# summits=FALSE     : use full peak regions (not summit ± fixed window)
#                     appropriate for broad H3K4me3 peaks
# score=DBA_SCORE_READS : raw read counts (DESeq2 handles normalisation)
# removeDup=TRUE    : exclude marked duplicate reads from counting
#                     (duplicates were marked by Picard in Stage 4)
# =============================================================================
cat("── Step 4: Count reads in consensus peaks ──\n")
cat("Using consensus peak set:", CONSENSUS_BED, "\n")
cat("This step takes ~20-40 minutes for 13 BAMs...\n\n")

# Count reads across all samples using consensus peaks
dba_counted <- dba.count(
  dba_obj,
  peaks      = CONSENSUS_BED,
  score      = DBA_SCORE_READS,
  summits    = FALSE,
  bRemoveDuplicates = TRUE,
  bParallel  = FALSE        # set TRUE if BiocParallel is configured
)

cat("\nPost-count summary:\n")
print(dba_counted)

# Save counted object (checkpoint — counting is the slowest step)
saveRDS(dba_counted, file.path(OUT_DIR, "GSE68952_dba_counted.rds"))
cat("[OK] Counted DBA object saved:", file.path(OUT_DIR, "GSE68952_dba_counted.rds"), "\n\n")

# Write raw counts matrix
counts_mat <- dba.peakset(dba_counted, bRetrieve=TRUE)
counts_df  <- as.data.frame(counts_mat)
write.csv(counts_df,
          file.path(OUT_DIR, "GSE68952_diffbind_all_counts.csv"),
          row.names = FALSE)
cat("[OK] Raw counts written:", file.path(OUT_DIR, "GSE68952_diffbind_all_counts.csv"), "\n\n")

# =============================================================================
# STEP 5 — Normalise and set up contrast
# Method: DESeq2 (RLE normalisation — handles library size differences)
# Contrast: HD vs Control
# =============================================================================
cat("── Step 5: Normalise + set contrast ──\n")

dba_norm <- dba.normalize(
  dba_counted,
  method    = DBA_DESEQ2,
  normalize = DBA_NORM_RLE,    # RLE = median ratio normalisation
  library   = DBA_LIBSIZE_FULL # use full library size for scaling
)
cat("[OK] Normalisation complete (DESeq2/RLE)\n")

# Set contrast: HD vs Control
dba_contrast <- dba.contrast(
  dba_norm,
  minMembers = 2,
  categories = DBA_CONDITION,
  reorderMeta = list(Condition="Control")  # Control = reference level
)
cat("[OK] Contrast set: HD vs Control\n\n")

# =============================================================================
# STEP 6 — Differential binding analysis
# =============================================================================
cat("── Step 6: DESeq2 differential binding analysis ──\n")

dba_analyzed <- dba.analyze(
  dba_contrast,
  method = DBA_DESEQ2,
  bParallel = FALSE
)

cat("\nAnalysis summary:\n")
print(dba_analyzed)

# =============================================================================
# STEP 7 — Extract results
# =============================================================================
cat("\n── Step 7: Extract results ──\n")

# All significant peaks (FDR < 0.05)
db_results <- dba.report(
  dba_analyzed,
  method    = DBA_DESEQ2,
  th        = FDR_CUTOFF,
  bCounts   = TRUE,
  fold      = 0            # return all FDR < 0.05 regardless of FC
)

n_total <- length(db_results)
n_gain  <- sum(db_results$Fold > 0)   # gained in HD (up)
n_loss  <- sum(db_results$Fold < 0)   # lost in HD (down)
pct_loss <- round(n_loss / n_total * 100, 1)

cat(sprintf("\nDifferential peaks (FDR < %.2f):\n", FDR_CUTOFF))
cat(sprintf("  Total differential : %d\n", n_total))
cat(sprintf("  Gained in HD (up)  : %d  (%.1f%%)\n", n_gain,  n_gain/n_total*100))
cat(sprintf("  Lost in HD (down)  : %d  (%.1f%%)\n", n_loss,  pct_loss))
cat(sprintf("  Paper expects      : ~2,830 differential, ~55%% down\n\n"))

# Validate key loci
cat("── Validation: checking published positive controls ──\n")
# Use GRanges-native overlap for validation — more robust than data.frame indexing
suppressPackageStartupMessages(library(GenomicRanges))

loci_gr <- GRanges(seqnames = c("chr11","chr11","chr8","chr4"),
                   ranges   = IRanges(
                     start = c(27700000, 113410000, 57094000, 3074000),
                     end   = c(27730000, 113440000, 57110000, 3085000)),
                   gene    = c("BDNF","DRD2","PENK","HTT"))

# Get all peaks with stats (FDR threshold = 1 = return everything)
all_gr <- dba.report(dba_analyzed, method=DBA_DESEQ2, th=1, bCounts=FALSE)

for (i in seq_along(loci_gr)) {
  gene   <- loci_gr$gene[i]
  region <- loci_gr[i]
  hits_all <- subsetByOverlaps(all_gr, region)
  hits_sig <- subsetByOverlaps(db_results, region)
  if (length(hits_sig) > 0) {
    direction <- ifelse(hits_sig$Fold[1] < 0, "DOWN in HD", "UP in HD")
    cat(sprintf("  [OK]   %s: %d differential peak(s) — %s (Fold=%.2f, FDR=%.4f)\n",
                gene, length(hits_sig), direction,
                hits_sig$Fold[1], hits_sig$FDR[1]))
  } else if (length(hits_all) > 0) {
    cat(sprintf("  [NOTE] %s: peak present but not differential (FDR=%.3f)\n",
                gene, hits_all$FDR[1]))
  } else {
    cat(sprintf("  [WARN] %s: no peak in consensus set at this locus\n", gene))
  }
}
cat("\n")

# =============================================================================
# STEP 8 — Write result files
# =============================================================================
cat("── Step 8: Write output files ──\n")

# Full results CSV
results_full <- as.data.frame(dba.report(dba_analyzed,
                                          method=DBA_DESEQ2,
                                          th=1,
                                          bCounts=TRUE,
                                          fold=0))
write.csv(results_full,
          file.path(OUT_DIR, "GSE68952_diffbind_all_peaks.csv"),
          row.names = FALSE)
cat("[OK] All peaks (with stats):", file.path(OUT_DIR, "GSE68952_diffbind_all_peaks.csv"), "\n")

# Significant results CSV
results_sig_df <- as.data.frame(db_results)
write.csv(results_sig_df,
          file.path(OUT_DIR, "GSE68952_diffbind_results.csv"),
          row.names = FALSE)
cat("[OK] Significant results:", file.path(OUT_DIR, "GSE68952_diffbind_results.csv"), "\n")

# Gained peaks BED (up in HD)
gained <- db_results[db_results$Fold > FC_CUTOFF]
if (length(gained) > 0) {
  gained_df <- as.data.frame(gained)
  gained_bed <- data.frame(
    chr   = as.character(seqnames(gained)),
    start = start(gained),
    end   = end(gained),
    name  = paste0("gained_", seq_len(length(gained))),
    score = round(-log10(gained$FDR), 2),
    strand = "."
  )
  write.table(gained_bed,
              file.path(OUT_DIR, "GSE68952_diffbind_gain.bed"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  cat(sprintf("[OK] Gained peaks BED (%d peaks): GSE68952_diffbind_gain.bed\n", nrow(gained_bed)))
}

# Lost peaks BED (down in HD)
lost <- db_results[db_results$Fold < -FC_CUTOFF]
if (length(lost) > 0) {
  lost_bed <- data.frame(
    chr   = as.character(seqnames(lost)),
    start = start(lost),
    end   = end(lost),
    name  = paste0("lost_", seq_len(length(lost))),
    score = round(-log10(lost$FDR), 2),
    strand = "."
  )
  write.table(lost_bed,
              file.path(OUT_DIR, "GSE68952_diffbind_loss.bed"),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  cat(sprintf("[OK] Lost peaks BED (%d peaks): GSE68952_diffbind_loss.bed\n", nrow(lost_bed)))
}
cat("\n")

# =============================================================================
# STEP 9 — Plots
# =============================================================================
cat("── Step 9: Generate plots ──\n")

# 1. MA plot
png(file.path(PLOT_DIR, "03_MA_plot.png"), width=800, height=600, res=100)
dba.plotMA(dba_analyzed, method=DBA_DESEQ2)
dev.off()
cat("[OK] MA plot: 03_MA_plot.png\n")

# 2. Volcano plot (ggplot2)
volcano_df <- results_full
volcano_df$chr <- as.character(volcano_df[,1])
volcano_df$sig <- "NS"
volcano_df$sig[volcano_df$FDR < FDR_CUTOFF & volcano_df$Fold >  FC_CUTOFF] <- "Gained"
volcano_df$sig[volcano_df$FDR < FDR_CUTOFF & volcano_df$Fold < -FC_CUTOFF] <- "Lost"
volcano_df$neg_log10_FDR <- pmin(-log10(volcano_df$FDR + 1e-300), 50)

p_volcano <- ggplot(volcano_df, aes(x=Fold, y=neg_log10_FDR, colour=sig)) +
  geom_point(alpha=0.4, size=0.8) +
  scale_colour_manual(values=c("Gained"="#E41A1C","Lost"="#377EB8","NS"="grey60")) +
  geom_vline(xintercept=c(-FC_CUTOFF, FC_CUTOFF), linetype="dashed", colour="black") +
  geom_hline(yintercept=-log10(FDR_CUTOFF), linetype="dashed", colour="black") +
  labs(title="GSE68952 Volcano Plot — HD vs Control (H3K4me3)",
       x="log2 Fold Change (HD/Control)",
       y="-log10(FDR)",
       colour="Direction") +
  theme_bw(base_size=13) +
  annotate("text", x=max(volcano_df$Fold)*0.8, y=45,
           label=paste0("Gained: ", n_gain), colour="#E41A1C", size=4) +
  annotate("text", x=min(volcano_df$Fold)*0.8, y=45,
           label=paste0("Lost: ", n_loss), colour="#377EB8", size=4)

ggsave(file.path(PLOT_DIR, "04_volcano_plot.png"), p_volcano,
       width=8, height=6, dpi=150)
cat("[OK] Volcano plot: 04_volcano_plot.png\n")

# 3. Correlation heatmap (read counts)
png(file.path(PLOT_DIR, "05_correlation_heatmap_counts.png"),
    width=800, height=700, res=100)
dba.plotHeatmap(dba_analyzed,
                method    = DBA_DESEQ2,
                contrast  = 1,
                correlations = TRUE)
dev.off()
cat("[OK] Correlation heatmap (counts): 05_correlation_heatmap_counts.png\n")

# 4. PCA (read counts, post-normalisation)
png(file.path(PLOT_DIR, "06_PCA_counts.png"), width=700, height=600, res=100)
dba.plotPCA(dba_analyzed,
            method     = DBA_DESEQ2,
            contrast   = 1,
            attributes = DBA_CONDITION,
            label      = DBA_ID)
dev.off()
cat("[OK] PCA (counts): 06_PCA_counts.png\n")

# 5. Box plot of normalised counts
png(file.path(PLOT_DIR, "07_boxplot_counts.png"), width=900, height=600, res=100)
dba.plotBox(dba_analyzed,
            method  = DBA_DESEQ2)
dev.off()
cat("[OK] Box plot: 07_boxplot_counts.png\n")

# 6. Heatmap of top 100 differential peaks
if (n_total >= 10) {
  top_n <- min(100, n_total)
  top_peaks <- db_results[order(db_results$FDR)[1:top_n]]

  png(file.path(PLOT_DIR, "08_heatmap_top_differential.png"),
      width=900, height=1000, res=100)
  dba.plotHeatmap(dba_analyzed,
                  method   = DBA_DESEQ2,
                  contrast = 1,
                  correlations = FALSE,
                  sites    = top_peaks)
  dev.off()
  cat(sprintf("[OK] Top %d peaks heatmap: 08_heatmap_top_differential.png\n", top_n))
}
cat("\n")

# =============================================================================
# STEP 10 — Final summary
# =============================================================================
cat("=============================================================\n")
cat("  STAGE 7 SUMMARY\n")
cat("=============================================================\n\n")
cat(sprintf("  Consensus peaks analysed : %d\n",   nrow(results_full)))
cat(sprintf("  Differential (FDR<0.05) : %d\n",    n_total))
cat(sprintf("  Gained in HD (up)       : %d (%.1f%%)\n", n_gain, n_gain/n_total*100))
cat(sprintf("  Lost in HD (down)       : %d (%.1f%%)\n", n_loss, pct_loss))
cat(sprintf("  Paper expected          : ~2,830 differential, ~55%% down\n\n"))

cat("  Output files:\n")
cat(sprintf("    %s/GSE68952_diffbind_results.csv    — significant peaks\n", OUT_DIR))
cat(sprintf("    %s/GSE68952_diffbind_all_peaks.csv  — all peaks with stats\n", OUT_DIR))
cat(sprintf("    %s/GSE68952_diffbind_gain.bed       — gained in HD\n", OUT_DIR))
cat(sprintf("    %s/GSE68952_diffbind_loss.bed       — lost in HD\n", OUT_DIR))
cat(sprintf("    %s/GSE68952_dba_counted.rds         — saved DBA object\n", OUT_DIR))
cat(sprintf("    %s/plots/                           — all QC + result plots\n\n", OUT_DIR))

cat("  Next step: Rscript stage8_annotation.R\n")
cat(sprintf("\n  %s: Stage 7 finished.\n", format(Sys.time())))
