#!/usr/bin/env bash
# =============================================================================
# GSE68952 ChIP-seq Pipeline — STAGE 6: Peak Calling (MACS3)
# =============================================================================
# Study    : H3K4me3 ChIP-seq · Huntington Disease vs Control · hg38
#
# STAGE 5 QC SUMMARY (informs peak calling decisions):
#   Fingerprint    : Strong J-curves all 12 ChIP samples — real enrichment
#   Correlation    : All ChIP-ChIP pairs r >= 0.89 — no outliers
#   PCA            : PC1=ChIP/Input, PC2 shows HD/Control biological signal
#   C3 depth note  : Highest read depth (8.49M) — may yield more peaks;
#                    DiffBind DESeq2 normalisation handles this
#   All 12 ChIP + 1 Input proceed — no exclusions
#
# PEAK CALLING STRATEGY:
#   Tool     : MACS3 callpeak
#   Mode     : --broad (H3K4me3 = broad promoter mark, ~1-3 kb domains)
#   Control  : SRR2028284 (shared Input DNA) for all 12 ChIP samples
#   Genome   : hs (hg38, effective size 2.7e9)
#   q-value  : 0.05 (strong enrichment confirmed by fingerprint)
#   Dups     : --keep-dup auto (MACS3 models local duplication itself)
#
# OUTPUTS PER SAMPLE (in 06_peaks/{SRR}/):
#   {name}_peaks.broadPeak   — peak calls (BED-like, 6+3 format)
#   {name}_peaks.gappedPeak  — gapped peak format (broad + narrow combined)
#   {name}_treat_pileup.bdg  — treatment pileup bedGraph
#   {name}_control_lambda.bdg— local lambda bedGraph
#   {name}_model.r           — R script for fragment size model plot
#
# CONSENSUS PEAK SET (in 06_peaks/consensus/):
#   Step 1: Merge all 12 per-sample broadPeak files
#   Step 2: Remove blacklist regions
#   Step 3: Keep only peaks present in >= 2 samples (reproducibility filter)
#   Step 4: Write final consensus BED for DiffBind
#
# INPUT  : 04_aligned/SRR*.markdup.bam
# OUTPUT : 06_peaks/
#
# REQUIREMENTS (chipseq_hg38 conda environment):
#   macs3 >= 2.2.9   bedtools >= 2.31   python >= 3.11
#
# HOW TO RUN:
#   conda activate chipseq_hg38
#   chmod +x stage6_peakcalling.sh
#   ./stage6_peakcalling.sh 2>&1 | tee run_stage6.log
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
INDIR="04_aligned"
OUTDIR="06_peaks"
GENOME_DIR="genome"
BLACKLIST="${GENOME_DIR}/hg38-blacklist.v2.bed"
SAMPLESHEET="samplesheet.csv"

# MACS3 parameters
MACS3_GENOME="hs"           # hg38 effective genome size shorthand
MACS3_QVALUE="0.05"         # FDR threshold for peak calling
MACS3_FORMAT="BAM"          # input format
MACS3_KEEPDUP="auto"        # let MACS3 model local duplication

# Consensus peak filter: minimum number of samples a region must appear in
# 2/12 = ~17% — standard reproducibility threshold for 6+6 design
CONSENSUS_MIN_OVERLAP=2

# Shared Input BAM
INPUT_BAM="${INDIR}/SRR2028284.markdup.bam"
INPUT_NAME="Input"

# ChIP samples in order (no Input)
declare -A CHIP_META
#            SRR           Name  Group
CHIP_META["SRR2028278"]="C1:Control"
CHIP_META["SRR2028279"]="C2:Control"
CHIP_META["SRR2028280"]="C3:Control"
CHIP_META["SRR2028281"]="C4:Control"
CHIP_META["SRR2028282"]="C5:Control"
CHIP_META["SRR2028283"]="C6:Control"
CHIP_META["SRR2028285"]="HD1:HD"
CHIP_META["SRR2028286"]="HD2:HD"
CHIP_META["SRR2028287"]="HD3:HD"
CHIP_META["SRR2028288"]="HD4:HD"
CHIP_META["SRR2028289"]="HD5:HD"
CHIP_META["SRR2028290"]="HD6:HD"

CHIP_ORDERED=(
    SRR2028278 SRR2028279 SRR2028280 SRR2028281 SRR2028282 SRR2028283
    SRR2028285 SRR2028286 SRR2028287 SRR2028288 SRR2028289 SRR2028290
)

# =============================================================================
# STEP 0 — Preflight
# =============================================================================
header "GSE68952 Stage 6: Peak Calling (MACS3)"
info "Date      : $(date)"
info "Genome    : ${MACS3_GENOME} (hg38)"
info "q-value   : ${MACS3_QVALUE}"
info "Mode      : --broad (H3K4me3)"
info "Control   : shared Input (${INPUT_BAM})"
info "Consensus : present in >= ${CONSENSUS_MIN_OVERLAP} samples"
echo ""

# Tool check — safe version query
info "Checking required tools..."
for tool in macs3 bedtools python3; do
    if command -v "$tool" &>/dev/null; then
        VER=$("$tool" --version 2>&1 | head -1) || VER="installed"
        success "  $tool : $VER"
    else
        error "$tool not found. Run: conda activate chipseq_hg38"
    fi
done
echo ""

# Check Input BAM
info "Checking Input BAM..."
if [[ -f "$INPUT_BAM" && -f "${INPUT_BAM}.bai" ]]; then
    SZ=$(du -sh "$INPUT_BAM" | cut -f1)
    success "  ${INPUT_BAM}  [${SZ}]"
else
    error "Input BAM missing or not indexed: ${INPUT_BAM}"
fi
echo ""

# Check ChIP BAMs
info "Checking ChIP BAMs..."
MISSING=0
for SRR in "${CHIP_ORDERED[@]}"; do
    BAM="${INDIR}/${SRR}.markdup.bam"
    META="${CHIP_META[$SRR]}"
    NM=$(echo "$META" | cut -d: -f1)
    if [[ -f "$BAM" && -f "${BAM}.bai" ]]; then
        SZ=$(du -sh "$BAM" | cut -f1)
        success "  ${NM} (${SRR})  [${SZ}]"
    else
        warn "  MISSING: $BAM"
        MISSING=$((MISSING + 1))
    fi
done
[[ $MISSING -gt 0 ]] && error "$MISSING BAM(s) missing. Run stage4_align.sh first."

# Check blacklist
if [[ ! -f "$BLACKLIST" ]]; then
    error "Blacklist not found: ${BLACKLIST}"
fi
echo ""

# =============================================================================
# STEP 1 — Create output directories
# =============================================================================
header "Step 1: Create directories"
mkdir -p "${OUTDIR}/consensus" "${OUTDIR}/logs"
for SRR in "${CHIP_ORDERED[@]}"; do
    META="${CHIP_META[$SRR]}"
    NM=$(echo "$META" | cut -d: -f1)
    mkdir -p "${OUTDIR}/${SRR}_${NM}"
done
success "Output directory: ${OUTDIR}/"
echo ""

# =============================================================================
# STEP 2 — MACS3 peak calling per sample
# --broad         : call broad peaks (appropriate for H3K4me3)
# --broad-cutoff  : q-value cutoff for broad region merging (same as -q)
# --nomodel       : skip fragment size modelling for SE ChIP-seq
#                   (model requires paired-end or specific fragment sizes)
#                   ALTERNATIVE: omit --nomodel and let MACS3 model it;
#                   both are valid but --nomodel is faster and more robust
#                   for 36 nt SE data
# --extsize 200   : extend reads to 200 bp (estimated fragment size for
#                   H3K4me3 from nucleosome-sized fragments ~147 bp + linker)
# --keep-dup auto : MACS3 calculates expected duplication from Poisson model
#                   and removes excess. Safe for ChIP-seq with marked dups.
# =============================================================================
header "Step 2: MACS3 peak calling (per sample)"
echo ""

CALL_COUNT=0
SKIP_COUNT=0

for SRR in "${CHIP_ORDERED[@]}"; do
    META="${CHIP_META[$SRR]}"
    NM=$(echo  "$META" | cut -d: -f1)
    GRP=$(echo "$META" | cut -d: -f2)
    BAM="${INDIR}/${SRR}.markdup.bam"
    SAMPLE_OUTDIR="${OUTDIR}/${SRR}_${NM}"
    PEAK_FILE="${SAMPLE_OUTDIR}/${NM}_peaks.broadPeak"
    MACS3_LOG="${OUTDIR}/logs/${SRR}_${NM}_macs3.log"

    echo "────────────────────────────────────────────────────────────"
    info "Processing: ${SRR}  |  ${NM}  |  ${GRP}"

    # Skip if peaks already called
    if [[ -f "$PEAK_FILE" ]]; then
        NPEAKS=$(wc -l < "$PEAK_FILE") || NPEAKS=0
        success "Already done (${NPEAKS} peaks) — skipping: ${NM}"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    info "  Running MACS3 callpeak..."
    macs3 callpeak \
        --treatment  "${BAM}" \
        --control    "${INPUT_BAM}" \
        --name       "${NM}" \
        --outdir     "${SAMPLE_OUTDIR}" \
        --format     "${MACS3_FORMAT}" \
        --gsize      "${MACS3_GENOME}" \
        --qvalue     "${MACS3_QVALUE}" \
        --broad \
        --broad-cutoff "${MACS3_QVALUE}" \
        --nomodel \
        --extsize    200 \
        --keep-dup   "${MACS3_KEEPDUP}" \
        2> "${MACS3_LOG}" || true

    if [[ -f "$PEAK_FILE" ]]; then
        NPEAKS=$(wc -l < "$PEAK_FILE") || NPEAKS=0
        success "  Peaks called: ${NPEAKS}  →  ${PEAK_FILE}"
    else
        warn "  broadPeak file not created — check log: ${MACS3_LOG}"
    fi

    CALL_COUNT=$((CALL_COUNT + 1))
    echo ""
done

# =============================================================================
# STEP 3 — Per-sample peak count summary
# =============================================================================
header "Step 3: Per-sample peak count summary"
echo ""
printf "%-14s %-8s %-10s %-12s\n" "SRR" "Sample" "Group" "Peaks"
printf '%0.s-' {1..50}; echo ""

CTRL_PEAKS=()
HD_PEAKS=()
ALL_PEAK_FILES=()

for SRR in "${CHIP_ORDERED[@]}"; do
    META="${CHIP_META[$SRR]}"
    NM=$(echo  "$META" | cut -d: -f1)
    GRP=$(echo "$META" | cut -d: -f2)
    PEAK_FILE="${OUTDIR}/${SRR}_${NM}/${NM}_peaks.broadPeak"

    if [[ -f "$PEAK_FILE" ]]; then
        NPEAKS=$(wc -l < "$PEAK_FILE") || NPEAKS=0
        ALL_PEAK_FILES+=("$PEAK_FILE")
        printf "%-14s %-8s %-10s %-12s\n" "$SRR" "$NM" "$GRP" "$NPEAKS"
        if [[ "$GRP" == "Control" ]]; then
            CTRL_PEAKS+=("$NPEAKS")
        else
            HD_PEAKS+=("$NPEAKS")
        fi
    else
        printf "%-14s %-8s %-10s %-12s\n" "$SRR" "$NM" "$GRP" "MISSING"
    fi
done

# Print group averages
if [[ ${#CTRL_PEAKS[@]} -gt 0 ]]; then
    CTRL_AVG=$(python3 -c "
vals=[${CTRL_PEAKS[*]:-0}]
print(f'{sum(vals)/len(vals):.0f}')
" 2>/dev/null) || CTRL_AVG="?"
    HD_AVG=$(python3 -c "
vals=[${HD_PEAKS[*]:-0}]
print(f'{sum(vals)/len(vals):.0f}')
" 2>/dev/null) || HD_AVG="?"
    echo ""
    info "  Control average peaks : ${CTRL_AVG}"
    info "  HD average peaks      : ${HD_AVG}"
    info "  Expected from paper   : ~28,608 total (combined)"
fi
echo ""

# =============================================================================
# STEP 4 — Build consensus peak set
# Strategy:
#   1. Concatenate all 12 per-sample broadPeak BED files
#   2. Sort and merge overlapping/adjacent peaks (bedtools merge)
#   3. Remove ENCODE blacklist regions (bedtools intersect -v)
#   4. Count how many original per-sample peaks overlap each merged region
#   5. Keep only regions covered by >= CONSENSUS_MIN_OVERLAP samples
#   6. Write final consensus BED (used as DiffBind peak set in Stage 7)
# =============================================================================
header "Step 4: Build consensus peak set"

MERGED_ALL="${OUTDIR}/consensus/all_samples_merged.bed"
MERGED_NOBLACKLIST="${OUTDIR}/consensus/all_samples_merged_noBL.bed"
CONSENSUS_BED="${OUTDIR}/consensus/GSE68952_consensus_peaks.bed"
CONSENSUS_ANNOTATED="${OUTDIR}/consensus/GSE68952_consensus_peaks_annotated.bed"

# 4a. Concatenate, sort, merge all peak files
info "  Concatenating and merging all per-sample peaks..."
if [[ ${#ALL_PEAK_FILES[@]} -eq 0 ]]; then
    # Re-collect in case STEP 3 was skipped
    for SRR in "${CHIP_ORDERED[@]}"; do
        META="${CHIP_META[$SRR]}"
        NM=$(echo "$META" | cut -d: -f1)
        PF="${OUTDIR}/${SRR}_${NM}/${NM}_peaks.broadPeak"
        [[ -f "$PF" ]] && ALL_PEAK_FILES+=("$PF")
    done
fi

# broadPeak cols: chr start end name score strand signalValue pValue qValue
# Take only chr, start, end (cols 1-3) for merging
cat "${ALL_PEAK_FILES[@]}" \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' \
    | grep -v "^#" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin \
    > "${MERGED_ALL}"

N_MERGED=$(wc -l < "$MERGED_ALL") || N_MERGED=0
success "  Merged regions: ${N_MERGED}"

# 4b. Remove blacklist regions
info "  Removing blacklist regions..."
bedtools intersect \
    -a "${MERGED_ALL}" \
    -b "${BLACKLIST}" \
    -v \
    > "${MERGED_NOBLACKLIST}"

N_NOBLACKLIST=$(wc -l < "$MERGED_NOBLACKLIST") || N_NOBLACKLIST=0
N_REMOVED=$((N_MERGED - N_NOBLACKLIST))
success "  After blacklist removal: ${N_NOBLACKLIST} (removed ${N_REMOVED})"

# 4c. Count per-sample overlaps and filter by minimum overlap count
info "  Filtering: keeping regions in >= ${CONSENSUS_MIN_OVERLAP} samples..."

python3 - << PYEOF || true
import subprocess, os

merged_file = "${MERGED_NOBLACKLIST}"
peak_files  = [$(printf '"%s",' "${ALL_PEAK_FILES[@]}"| sed 's/,$//')]
min_overlap = ${CONSENSUS_MIN_OVERLAP}
out_bed     = "${CONSENSUS_BED}"
out_ann     = "${CONSENSUS_ANNOTATED}"

sample_names = []
for pf in peak_files:
    # Extract sample name from path: 06_peaks/SRR_NAME/NAME_peaks.broadPeak
    nm = os.path.basename(os.path.dirname(pf)).split('_',1)[1]
    sample_names.append(nm)

# Read merged regions
regions = []
with open(merged_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            regions.append((parts[0], int(parts[1]), int(parts[2])))

# For each merged region, count how many sample peak files overlap it
# Use bedtools intersect -c for efficiency
import tempfile, subprocess

# Write merged BED to temp file (already done, use merged_file)
# Run bedtools intersect -c for each sample peak file
counts = [0] * len(regions)

for pf in peak_files:
    result = subprocess.run(
        ['bedtools', 'intersect', '-a', merged_file, '-b', pf, '-c'],
        capture_output=True, text=True
    )
    for i, line in enumerate(result.stdout.strip().split('\n')):
        if line:
            parts = line.split('\t')
            if len(parts) >= 4 and i < len(counts):
                c = int(parts[3])
                counts[i] += (1 if c > 0 else 0)

# Write consensus BED (regions in >= min_overlap samples)
kept = 0
with open(out_bed, 'w') as out, open(out_ann, 'w') as ann:
    # Write header for annotated file
    ann.write('chr\tstart\tend\twidth\tn_samples\n')
    for i, (chrom, start, end) in enumerate(regions):
        n = counts[i]
        if n >= min_overlap:
            width = end - start
            out.write(f'{chrom}\t{start}\t{end}\n')
            ann.write(f'{chrom}\t{start}\t{end}\t{width}\t{n}\n')
            kept += 1

print(f"  Consensus peaks (>= {min_overlap} samples): {kept}")
print(f"  Removed (< {min_overlap} samples): {len(regions) - kept}")
PYEOF

if [[ -f "$CONSENSUS_BED" ]]; then
    N_CONSENSUS=$(wc -l < "$CONSENSUS_BED") || N_CONSENSUS=0
    success "  Consensus peak set: ${N_CONSENSUS} peaks → ${CONSENSUS_BED}"
else
    warn "  Consensus BED not created — check python3 output above"
fi
echo ""

# =============================================================================
# STEP 5 — Validate key loci (positive controls from paper)
# Check that H3K4me3 peaks are present at BDNF, DRD2, PENK in Controls
# These genes are published as having reduced H3K4me3 in HD
# =============================================================================
header "Step 5: Validate positive control loci"
echo ""

# hg38 coordinates for TSS ± 2 kb of published HD H3K4me3 loci
declare -A LOCI
LOCI["BDNF"]="chr11:27700000-27730000"
LOCI["DRD2"]="chr11:113410000-113440000"
LOCI["PENK"]="chr8:57094000-57110000"
LOCI["HTT"]="chr4:3074000-3085000"

info "Checking consensus peaks at published H3K4me3 loci (hg38)..."
echo ""

for GENE in BDNF DRD2 PENK HTT; do
    COORDS="${LOCI[$GENE]}"
    CHR=$(echo "$COORDS" | cut -d: -f1)
    START=$(echo "$COORDS" | cut -d- -f1 | cut -d: -f2)
    END=$(echo "$COORDS" | cut -d- -f2)

    if [[ -f "$CONSENSUS_BED" ]]; then
        NHITS=$(awk -v c="$CHR" -v s="$START" -v e="$END" \
            '$1==c && $3>s && $2<e {count++} END{print count+0}' \
            "$CONSENSUS_BED") || NHITS=0
        if [[ "$NHITS" -gt 0 ]]; then
            success "  ${GENE} (${COORDS}): ${NHITS} consensus peak(s) — PRESENT"
        else
            warn "  ${GENE} (${COORDS}): 0 peaks — check coordinates or peak calling"
        fi
    fi
done
echo ""

# =============================================================================
# STEP 6 — Final summary
# =============================================================================
header "Step 6: Stage 6 Summary"
echo ""
echo "  Per-sample peaks   : ${OUTDIR}/SRR*_*/  (*_peaks.broadPeak)"
echo "  MACS3 logs         : ${OUTDIR}/logs/"
echo "  Merged all peaks   : ${MERGED_ALL}"
echo "  Blacklist filtered : ${MERGED_NOBLACKLIST}"
echo "  Consensus peak set : ${CONSENSUS_BED}"
echo "  Annotated counts   : ${CONSENSUS_ANNOTATED}"
echo ""
echo "  Called this run    : ${CALL_COUNT}"
echo "  Already existed    : ${SKIP_COUNT}"

if [[ -f "$CONSENSUS_BED" ]]; then
    N_CONSENSUS=$(wc -l < "$CONSENSUS_BED") || N_CONSENSUS=0
    echo ""
    success "  Final consensus peaks : ${N_CONSENSUS}"
    info    "  Paper reports         : ~28,608 peaks total"
    if [[ "$N_CONSENSUS" -gt 15000 && "$N_CONSENSUS" -lt 50000 ]]; then
        success "  Peak count is in expected range."
    else
        warn "  Peak count outside expected range — review MACS3 parameters."
    fi
fi

echo ""
success "Stage 6 complete."
echo ""
info "Next step: run stage7_diffbind.R"
info "  Input  : ${INDIR}/*.markdup.bam  +  ${OUTDIR}/SRR*/*_peaks.broadPeak"
info "  Command: Rscript stage7_diffbind.R"
echo ""
info "$(date): Stage 6 finished."
echo ""

# =============================================================================
# NOTES
# =============================================================================
cat << 'NOTES'
-----------------------------------------------------------------------------
STAGE 6 NOTES — GSE68952 MACS3 Peak Calling
-----------------------------------------------------------------------------
WHY --broad FOR H3K4me3:
  H3K4me3 marks active promoters with domains spanning 1-3 kb.
  --broad merges nearby enriched regions into contiguous domains.
  Do NOT use for TF ChIP-seq (use narrow peaks for CTCF, NFKB etc.)

WHY --nomodel --extsize 200:
  Single-end 36 nt reads cannot be used to build a shift model
  (MACS3 model requires both +/- strand pileups at the same loci,
  which works best with longer reads or paired-end data).
  extsize 200 = estimated fragment size (nucleosome ~147 bp + ~50 bp linker)
  This is standard practice for short SE ChIP-seq reads.

WHY --keep-dup auto:
  Duplicates are MARKED in the BAM (not removed).
  MACS3 auto calculates the expected number of duplicate reads
  from a Poisson model and removes only the excess.
  This is more principled than a hard duplicate removal step.

CONSENSUS PEAK SET RATIONALE:
  Individual sample peaks are noisy. The consensus set:
    - Removes peaks unique to one sample (likely noise)
    - Retains peaks reproducible in >= 2 samples
    - Provides a fixed set of regions for DiffBind counting
  DiffBind will re-count reads in ALL consensus regions for all 12 samples,
  regardless of whether a given sample called a peak there or not.

C3 DEPTH NOTE (from Stage 5 PCA):
  C3 has the highest usable read depth (8.49M) and may call more peaks.
  This is expected and not a problem — check C3 peak count in Step 3.
  If C3 has > 2x more peaks than other Controls, this is a depth effect
  and DiffBind's TMM/DESeq2 normalisation corrects for it.

LOW-READ SAMPLES (C4: 2.25M, HD2: 2.47M, HD4: 2.09M):
  May call fewer peaks. This is fine — DiffBind uses the consensus
  peak set (counted in all samples) not the per-sample peak files
  for differential analysis.

EXPECTED PEAK COUNTS PER SAMPLE (H3K4me3, hg38):
  Typical range: 15,000 - 40,000 broad peaks per sample
  Low-read samples: 8,000 - 20,000 expected
  Paper total (all samples merged): ~28,608

VALIDATION LOCI:
  BDNF, DRD2, PENK, HTT should all have H3K4me3 peaks in Controls.
  Reduced signal in HD samples is detected in Stage 7 (DiffBind),
  not visible in raw peak presence/absence.
-----------------------------------------------------------------------------
NOTES
