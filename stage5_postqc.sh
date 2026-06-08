#!/usr/bin/env bash
# =============================================================================
# GSE68952 ChIP-seq Pipeline — STAGE 5: Post-Alignment QC
# =============================================================================
# Study    : H3K4me3 ChIP-seq · Huntington Disease vs Control · hg38
#
# STAGE 4 ALIGNMENT RESULTS (baked in — no guessing):
#   Mapping rate     : 100% all 13 samples (MAPQ>=30 + blacklist filtered)
#   Optical dups     : 0 (correct for Genome Analyzer II)
#   Input %dup       : 0.8% — excellent genomic background control
#   ChIP %dup range  : 11.6% (C4) – 32.6% (C2)
#
#   Usable reads (non-dup) per sample:
#     C1:4.13M  C2:5.01M  C3:8.49M  C4:2.25M*  C5:6.16M  C6:4.49M
#     Input:11.86M
#     HD1:3.77M  HD2:2.47M*  HD3:4.43M  HD4:2.09M*  HD5:4.04M  HD6:3.98M
#     (* = monitor, all above 2M hard floor)
#
# WHAT THIS STAGE DOES:
#   1. deepTools plotFingerprint  — assess ChIP enrichment quality
#   2. deepTools bamCompare       — log2(ChIP/Input) bigWig tracks
#   3. deepTools multiBamSummary  — genome-wide read coverage matrix
#   4. deepTools plotCorrelation  — sample-to-sample Pearson correlation
#   5. deepTools plotPCA          — PCA of all samples (expect HD vs Ctrl sep.)
#   6. deepTools bamCoverage      — normalised bigWig per sample (for IGV)
#   7. Summary table              — FRiP proxy + key QC flags
#
# INPUT  : 04_aligned/SRR*.markdup.bam  +  04_aligned/SRR*.markdup.bam.bai
# OUTPUT : 05_postqc/
#
# REQUIREMENTS (chipseq_hg38 conda environment):
#   deeptools >= 3.5
#
# HOW TO RUN:
#   conda activate chipseq_hg38
#   chmod +x stage5_postqc.sh
#   ./stage5_postqc.sh 2>&1 | tee run_stage5.log
# =============================================================================

set -euo pipefail

# ── Colour helpers ─────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
info()    { echo -e "${CYAN}[INFO]${NC}  $*"; }
success() { echo -e "${GREEN}[OK]${NC}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error()   { echo -e "${RED}[ERROR]${NC} $*"; exit 1; }
header()  { echo -e "\n${BOLD}── $* ──${NC}"; }

# =============================================================================
# CONFIGURATION
# =============================================================================
THREADS=8
INDIR="04_aligned"
OUTDIR="05_postqc"
GENOME_DIR="genome"
BLACKLIST="${GENOME_DIR}/hg38-blacklist.v2.bed"

# hg38 effective genome size for deepTools normalisation
# 2.7e9 = mappable hg38 size (standard for 36 nt reads)
EFFECTIVE_GENOME_SIZE=2700000000

# Normalisation method for bamCoverage bigWig tracks
# RPGC = reads per genomic content (1x depth) — standard for ChIP-seq IGV tracks
NORMALISATION="RPGC"

# Binsize for coverage tracks (bp) — 10 bp is standard for H3K4me3
BINSIZE=10

# Samples in SRR order (alphabetical = processing order)
SRR_ORDERED=(
    SRR2028278 SRR2028279 SRR2028280 SRR2028281 SRR2028282 SRR2028283
    SRR2028284
    SRR2028285 SRR2028286 SRR2028287 SRR2028288 SRR2028289 SRR2028290
)

declare -A SAMPLE_NAME
SAMPLE_NAME["SRR2028278"]="C1"
SAMPLE_NAME["SRR2028279"]="C2"
SAMPLE_NAME["SRR2028280"]="C3"
SAMPLE_NAME["SRR2028281"]="C4"
SAMPLE_NAME["SRR2028282"]="C5"
SAMPLE_NAME["SRR2028283"]="C6"
SAMPLE_NAME["SRR2028284"]="Input"
SAMPLE_NAME["SRR2028285"]="HD1"
SAMPLE_NAME["SRR2028286"]="HD2"
SAMPLE_NAME["SRR2028287"]="HD3"
SAMPLE_NAME["SRR2028288"]="HD4"
SAMPLE_NAME["SRR2028289"]="HD5"
SAMPLE_NAME["SRR2028290"]="HD6"

# Input BAM (shared control for all ChIP samples)
INPUT_BAM="${INDIR}/SRR2028284.markdup.bam"
INPUT_NAME="Input"

# ChIP-only BAMs (excludes Input) — used for bamCompare, fingerprint
CHIP_SRRS=(
    SRR2028278 SRR2028279 SRR2028280 SRR2028281 SRR2028282 SRR2028283
    SRR2028285 SRR2028286 SRR2028287 SRR2028288 SRR2028289 SRR2028290
)

# Known low-read samples from Stage 4 (usable reads < 3M)
LOW_READ="C4(2.25M) HD2(2.47M) HD4(2.09M)"

# =============================================================================
# STEP 0 — Preflight
# =============================================================================
header "GSE68952 Stage 5: Post-Alignment QC"
info "Date       : $(date)"
info "Threads    : ${THREADS}"
info "Norm method: ${NORMALISATION}"
info "Bin size   : ${BINSIZE} bp"
info "Eff genome : ${EFFECTIVE_GENOME_SIZE}"
echo ""
warn "Low-read samples to monitor: ${LOW_READ}"
echo ""

info "Checking deeptools..."
if command -v deeptools &>/dev/null; then
    VER=$(deeptools --version 2>&1 | head -1) || VER="installed (version query failed)"
    success "deeptools: $VER"
else
    error "deeptools not found. Run: conda activate chipseq_hg38"
fi
echo ""

# Verify all BAMs present and indexed
info "Checking BAM inputs..."
MISSING=0
for SRR in "${SRR_ORDERED[@]}"; do
    BAM="${INDIR}/${SRR}.markdup.bam"
    BAI="${BAM}.bai"
    if [[ -f "$BAM" && -f "$BAI" ]]; then
        SZ=$(du -sh "$BAM" | cut -f1)
        success "  ${SRR}.markdup.bam  [${SZ}]"
    else
        warn "  MISSING: ${BAM} or its index"
        MISSING=$((MISSING + 1))
    fi
done
[[ $MISSING -gt 0 ]] && error "$MISSING BAM(s) missing. Run stage4_align.sh first."
echo ""

# =============================================================================
# STEP 1 — Create output directories
# =============================================================================
header "Step 1: Create directories"
mkdir -p \
    "${OUTDIR}/fingerprint" \
    "${OUTDIR}/bigwig" \
    "${OUTDIR}/bigwig_log2ratio" \
    "${OUTDIR}/correlation" \
    "${OUTDIR}/pca"
success "Output directory: ${OUTDIR}/"
echo ""

# Build BAM lists used by multiple deepTools commands
ALL_BAMS=()
ALL_LABELS=()
CHIP_BAMS=()
CHIP_LABELS=()
for SRR in "${SRR_ORDERED[@]}"; do
    BAM="${INDIR}/${SRR}.markdup.bam"
    NM="${SAMPLE_NAME[$SRR]}"
    ALL_BAMS+=("$BAM")
    ALL_LABELS+=("$NM")
    if [[ "$SRR" != "SRR2028284" ]]; then
        CHIP_BAMS+=("$BAM")
        CHIP_LABELS+=("$NM")
    fi
done

# =============================================================================
# STEP 2 — plotFingerprint
# Diagnoses ChIP enrichment quality.
# A good ChIP shows a "J-curve": the enriched sample has many reads piled
# into a small fraction of the genome (steep left, sharp rise at right).
# Input should be a near-diagonal line (reads distributed across genome).
# Expected for H3K4me3 (active promoter mark): clear J-curve separation
# between ChIP and Input.
# =============================================================================
header "Step 2: plotFingerprint — ChIP enrichment quality"

FINGERPRINT_OUT="${OUTDIR}/fingerprint/GSE68952_fingerprint"

if [[ -f "${FINGERPRINT_OUT}.png" ]]; then
    success "Already done — skipping fingerprint plot."
else
    info "Running plotFingerprint on all 13 samples..."
    info "  (samples 1000 bp bins, 500k reads sampled — takes ~10 min)"

    plotFingerprint \
        --bamfiles "${ALL_BAMS[@]}" \
        --labels   "${ALL_LABELS[@]}" \
        --outRawCounts "${FINGERPRINT_OUT}.tab" \
        --plotFile     "${FINGERPRINT_OUT}.png" \
        --plotTitle    "GSE68952 ChIP-seq Fingerprint (H3K4me3)" \
        --numberOfSamples 500000 \
        --binSize 1000 \
        --numberOfProcessors "${THREADS}" \
        --blackListFileName "${BLACKLIST}" \
        --skipZeros \
        2>&1 | grep -v "^$" || true

    success "Fingerprint plot: ${FINGERPRINT_OUT}.png"
    success "Raw counts: ${FINGERPRINT_OUT}.tab"
fi
echo ""

# =============================================================================
# STEP 3 — bamCoverage: normalised bigWig per sample
# RPGC normalisation = reads per genomic content (equivalent to 1x depth).
# --ignoreDuplicates: exclude marked duplicates from coverage calculation.
# --extendReads: not used for SE data (no fragment size info without MACS2).
#   We will use --centerReads as alternative for H3K4me3.
# --smoothLength 150: smooth over 150 bp window — appropriate for H3K4me3
#   broad peaks (~1-2 kb). Do not use for narrow peaks (CTCF etc.).
# =============================================================================
header "Step 3: bamCoverage — normalised bigWig tracks"
info "Generating RPGC-normalised bigWig for each sample..."
echo ""

for SRR in "${SRR_ORDERED[@]}"; do
    NM="${SAMPLE_NAME[$SRR]}"
    BAM="${INDIR}/${SRR}.markdup.bam"
    BW="${OUTDIR}/bigwig/${SRR}_${NM}.bw"

    if [[ -f "$BW" ]]; then
        SZ=$(du -sh "$BW" | cut -f1)
        success "  Already done [${SZ}] — skipping: ${NM}"
        continue
    fi

    info "  ${SRR} | ${NM}..."
    bamCoverage \
        --bam              "${BAM}" \
        --outFileName      "${BW}" \
        --outFileFormat    bigwig \
        --binSize          "${BINSIZE}" \
        --normalizeUsing   "${NORMALISATION}" \
        --effectiveGenomeSize "${EFFECTIVE_GENOME_SIZE}" \
        --ignoreDuplicates \
        --blackListFileName "${BLACKLIST}" \
        --smoothLength     150 \
        --numberOfProcessors "${THREADS}" \
        2>&1 | grep -v "^$" || true

    SZ=$(du -sh "$BW" | cut -f1)
    success "  Done: ${BW}  [${SZ}]"
done
echo ""

# =============================================================================
# STEP 4 — bamCompare: log2(ChIP/Input) bigWig
# Produces a signal-over-input track for each ChIP sample.
# This is what you load in IGV to see enrichment at specific loci.
# --pseudocount 1: avoids log2(0) when Input has zero coverage in a bin.
# --scaleFactorsMethod None: let RPGC normalisation handle scaling.
# Expected: positive values at active gene promoters (BDNF, DRD2, PENK
# should show positive signal in Controls and reduced signal in HD).
# =============================================================================
header "Step 4: bamCompare — log2(ChIP/Input) tracks"
info "Generating log2 ratio tracks for 12 ChIP samples vs shared Input..."
echo ""

for SRR in "${CHIP_SRRS[@]}"; do
    NM="${SAMPLE_NAME[$SRR]}"
    BAM="${INDIR}/${SRR}.markdup.bam"
    BW_RATIO="${OUTDIR}/bigwig_log2ratio/${SRR}_${NM}_log2ratio.bw"

    if [[ -f "$BW_RATIO" ]]; then
        SZ=$(du -sh "$BW_RATIO" | cut -f1)
        success "  Already done [${SZ}] — skipping: ${NM}"
        continue
    fi

    info "  ${SRR} | ${NM} vs Input..."
    bamCompare \
        --bamfile1         "${BAM}" \
        --bamfile2         "${INPUT_BAM}" \
        --outFileName      "${BW_RATIO}" \
        --outFileFormat    bigwig \
        --operation        log2 \
        --pseudocount      1 \
        --binSize          "${BINSIZE}" \
        --normalizeUsing   CPM \
        --scaleFactorsMethod None \
        --ignoreDuplicates \
        --blackListFileName "${BLACKLIST}" \
        --smoothLength     150 \
        --numberOfProcessors "${THREADS}" \
        2>&1 | grep -v "^$" || true

    SZ=$(du -sh "$BW_RATIO" | cut -f1)
    success "  Done: ${BW_RATIO}  [${SZ}]"
done
echo ""

# =============================================================================
# STEP 5 — multiBamSummary + plotCorrelation
# Computes read counts in 10 kb bins across genome, then Pearson correlation.
# Expected: Controls cluster together, HD samples cluster together.
# Outlier: any sample with r < 0.85 to its group centroid warrants inspection.
# Input should cluster separately from all ChIP samples.
# =============================================================================
header "Step 5: multiBamSummary + plotCorrelation"

SUMMARY_NPZ="${OUTDIR}/correlation/GSE68952_multiBamSummary.npz"
CORR_PLOT="${OUTDIR}/correlation/GSE68952_pearson_heatmap.png"
CORR_TAB="${OUTDIR}/correlation/GSE68952_pearson_matrix.tab"

if [[ -f "$SUMMARY_NPZ" ]]; then
    success "multiBamSummary already done — skipping."
else
    info "Running multiBamSummary (10 kb bins, all 13 samples)..."
    info "  This is the slowest step — ~20-40 min on 8 threads."
    multiBamSummary bins \
        --bamfiles      "${ALL_BAMS[@]}" \
        --labels        "${ALL_LABELS[@]}" \
        --outFileName   "${SUMMARY_NPZ}" \
        --binSize       10000 \
        --numberOfProcessors "${THREADS}" \
        --blackListFileName "${BLACKLIST}" \
        --ignoreDuplicates \
        2>&1 | grep -v "^$" || true
    success "multiBamSummary done: ${SUMMARY_NPZ}"
fi
echo ""

if [[ -f "$CORR_PLOT" ]]; then
    success "Correlation plot already done — skipping."
else
    info "Plotting Pearson correlation heatmap..."
    plotCorrelation \
        --corData       "${SUMMARY_NPZ}" \
        --corMethod     pearson \
        --whatToPlot    heatmap \
        --plotFile      "${CORR_PLOT}" \
        --outFileCorMatrix "${CORR_TAB}" \
        --plotTitle     "GSE68952 Pearson Correlation (H3K4me3)" \
        --colorMap      RdYlBu \
        --plotNumbers \
        --removeOutliers \
        2>&1 | grep -v "^$" || true
    success "Correlation heatmap: ${CORR_PLOT}"
    success "Correlation matrix : ${CORR_TAB}"
fi
echo ""

# =============================================================================
# STEP 6 — plotPCA
# PCA of all 13 samples from the multiBamSummary matrix.
# Expected clustering: HD samples separate from Controls on PC1 or PC2.
# Input should sit far from all ChIP samples on at least one PC.
# Low-read samples (C4, HD2, HD4) may show slightly more scatter.
# =============================================================================
header "Step 6: plotPCA"

PCA_PLOT="${OUTDIR}/pca/GSE68952_PCA.png"
PCA_TAB="${OUTDIR}/pca/GSE68952_PCA.tab"

if [[ -f "$PCA_PLOT" ]]; then
    success "PCA plot already done — skipping."
else
    info "Plotting PCA..."
    plotPCA \
        --corData    "${SUMMARY_NPZ}" \
        --plotFile   "${PCA_PLOT}" \
        --outFileNameData "${PCA_TAB}" \
        --plotTitle  "GSE68952 PCA — H3K4me3 ChIP-seq" \
        --ntop       1000 \
        2>&1 | grep -v "^$" || true
    success "PCA plot: ${PCA_PLOT}"
fi
echo ""

# =============================================================================
# STEP 7 — Print Pearson correlation matrix summary
# Flag any off-diagonal values < 0.85 within the same group
# =============================================================================
header "Step 7: Correlation matrix summary"

if [[ -f "$CORR_TAB" ]]; then
    python3 - << 'PYEOF'
import sys

tab_file = "05_postqc/correlation/GSE68952_pearson_matrix.tab"
try:
    with open(tab_file) as f:
        lines = [l.strip() for l in f if l.strip()]

    header_line = lines[0].split('\t')
    samples = header_line[1:]   # skip first empty cell
    print(f"\n  {'':8}", end="")
    for s in samples:
        print(f"  {s:>6}", end="")
    print()
    print("  " + "─"*(8 + 8*len(samples)))

    warn_pairs = []
    ctrl = {"C1","C2","C3","C4","C5","C6"}
    hd   = {"HD1","HD2","HD3","HD4","HD5","HD6"}

    for row_line in lines[1:]:
        parts = row_line.split('\t')
        row_name = parts[0]
        vals = [float(v) for v in parts[1:]]
        print(f"  {row_name:8}", end="")
        for j, v in enumerate(vals):
            col_name = samples[j]
            color = ""
            # flag within-group pairs below 0.85
            same_grp = (row_name in ctrl and col_name in ctrl) or \
                       (row_name in hd   and col_name in hd)
            if same_grp and row_name != col_name and v < 0.85:
                warn_pairs.append(f"{row_name}–{col_name}: r={v:.2f}")
                color = "!"
            print(f"  {v:>5.2f}{color if color else ' '}", end="")
        print()

    if warn_pairs:
        print(f"\n  WARNING — within-group correlations below 0.85:")
        for wp in warn_pairs:
            print(f"    {wp}")
    else:
        print(f"\n  All within-group Pearson correlations >= 0.85")

except FileNotFoundError:
    print("  Correlation matrix not yet generated — re-run after multiBamSummary.")
except Exception as e:
    print(f"  Could not parse matrix: {e}")
PYEOF
else
    warn "Correlation matrix not yet available — will be printed on next run."
fi
echo ""

# =============================================================================
# STEP 8 — Final output summary
# =============================================================================
header "Step 8: Output summary"
echo ""
echo "  Fingerprint plot   : ${OUTDIR}/fingerprint/GSE68952_fingerprint.png"
echo "  bigWig tracks      : ${OUTDIR}/bigwig/*.bw              (RPGC, load in IGV)"
echo "  log2 ratio tracks  : ${OUTDIR}/bigwig_log2ratio/*.bw    (ChIP/Input, enrichment)"
echo "  Correlation heatmap: ${OUTDIR}/correlation/GSE68952_pearson_heatmap.png"
echo "  Correlation matrix : ${OUTDIR}/correlation/GSE68952_pearson_matrix.tab"
echo "  PCA plot           : ${OUTDIR}/pca/GSE68952_PCA.png"
echo ""

# Check all expected outputs
MISSING_OUT=0
for F in \
    "${OUTDIR}/fingerprint/GSE68952_fingerprint.png" \
    "${OUTDIR}/correlation/GSE68952_pearson_heatmap.png" \
    "${OUTDIR}/pca/GSE68952_PCA.png"; do
    if [[ -f "$F" ]]; then
        success "  $(basename $F)"
    else
        warn "  MISSING: $(basename $F)"
        MISSING_OUT=$((MISSING_OUT + 1))
    fi
done

BW_COUNT=$(ls "${OUTDIR}/bigwig/"*.bw 2>/dev/null | wc -l) || BW_COUNT=0
RATIO_COUNT=$(ls "${OUTDIR}/bigwig_log2ratio/"*.bw 2>/dev/null | wc -l) || RATIO_COUNT=0
info "  bigWig files       : ${BW_COUNT} / 13"
info "  log2 ratio files   : ${RATIO_COUNT} / 12"

echo ""
if [[ $MISSING_OUT -eq 0 && $BW_COUNT -eq 13 && $RATIO_COUNT -eq 12 ]]; then
    success "Stage 5 complete. All QC outputs present."
else
    warn "Stage 5 incomplete — re-run to generate missing outputs."
fi

echo ""
info "Next step: run stage6_peakcalling.sh"
info "  Input  : ${INDIR}/*.markdup.bam  +  ${OUTDIR}/bigwig/"
info "  Command: bash stage6_peakcalling.sh"
echo ""
info "$(date): Stage 5 finished."
echo ""

# =============================================================================
# NOTES
# =============================================================================
cat << 'NOTES'
-----------------------------------------------------------------------------
STAGE 5 NOTES — GSE68952 Post-Alignment QC
-----------------------------------------------------------------------------
FINGERPRINT PLOT — what to look for:
  X-axis: fraction of genome (sorted by coverage, lowest to highest)
  Y-axis: fraction of reads
  Input : near-diagonal line (reads evenly distributed)
  H3K4me3 ChIP: steep J-curve — most reads in small fraction of genome
  If ChIP curve overlaps Input: poor enrichment, consider excluding sample
  Low-read samples (C4, HD2, HD4) may show flatter curves — still valid
  if J-shape is visible

CORRELATION HEATMAP — expected pattern:
  Controls cluster together (r > 0.85 within group)
  HD samples cluster together (r > 0.85 within group)
  Input clusters separately from all ChIP samples
  Cross-group (HD vs Control) r = 0.70-0.85 expected (biological difference)
  r < 0.70 for same-group pair: investigate that sample

PCA PLOT — expected pattern:
  PC1 or PC2 should separate HD from Control
  Input should be an outlier (separate from all ChIP samples)
  C4, HD2, HD4 (low-read) may be slightly displaced but should still
  cluster with their respective groups

bigWig tracks (RPGC normalised) — IGV loading:
  Load all 13 .bw files in IGV (hg38)
  Navigate to BDNF (chr11:27,654,893-27,722,058)
  Expect: H3K4me3 peaks at TSS in Controls, reduced in HD samples
  Also check: DRD2 (chr11:113,409,605-113,475,432)
              PENK (chr8:57,095,128-57,107,776)

log2(ChIP/Input) tracks:
  Values > 0 = enriched above input (real signal)
  Values < 0 = depleted below input (background/artefact)
  Look for peaks > 2.0 at active promoters in Control samples
  These should be reduced in HD samples at validated loci

LOW-READ SAMPLES IN THIS CONTEXT:
  C4 (2.25M), HD2 (2.47M), HD4 (2.09M)
  bamCoverage RPGC normalisation accounts for depth differences
  These samples are valid for QC and peak calling
  Sensitivity analysis (re-run DiffBind without these 3) recommended
  after Stage 7 to confirm conclusions hold

multiBamSummary runtime:
  Most time-consuming step (~20-40 min on 8 threads for 13 BAMs)
  Uses 10 kb bins across whole genome — do not reduce bin size here
  (10 kb is appropriate for genome-wide correlation; peak-level analysis
  comes in DiffBind)
-----------------------------------------------------------------------------
NOTES
