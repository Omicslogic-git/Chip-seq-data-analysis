#!/usr/bin/env bash
# =============================================================================
# GSE68952 ChIP-seq Pipeline — STAGE 4: Alignment to hg38
# =============================================================================
# Study    : H3K4me3 ChIP-seq · Huntington Disease vs Control · hg38
#
# POST-TRIM QC SUMMARY (informs alignment decisions):
#   Quality trimming    : 0.00% bases removed — data was already clean
#   Reads passing filter: ~100% all samples (loss < 0.01%)
#   Avg read length     : 35.5 nt (from 36.0 nt) — soft 3' clipping only
#   Low-read samples    : C4 (2.84M), HD4 (2.72M) — clean quality, proceed
#   Adapter note        : 26-40% reads had adapter detected but 0 bases
#                         trimmed — Illumina read-through at read end; normal
#   Poly-G              : 0.00% detected — Genome Analyzer II uses 4-channel
#                         chemistry; poly-G artifact does not apply here
#
# ALIGNMENT STRATEGY:
#   Aligner  : Bowtie2 --very-sensitive (most sensitive preset for ChIP-seq)
#   Index    : GRCh38_noalt_as (no ALT contigs — avoids read splitting)
#   Filter   : MAPQ >= 30 (unique-enough mappers only)
#   Dups     : Picard MarkDuplicates (mark, not remove — DiffBind handles)
#   Blacklist: Apply ENCODE hg38 blacklist v2 to remove artefact regions
#
# PER-SAMPLE OUTPUT (in 04_aligned/):
#   {SRR}.raw.bam              — unfiltered sorted BAM (deleted after step)
#   {SRR}.filtered.bam         — MAPQ>=30, no blacklist, sorted + indexed
#   {SRR}.markdup.bam          — duplicates marked (used for DiffBind/MACS2)
#   {SRR}.markdup.bam.bai      — BAM index
#   {SRR}.flagstat.txt         — alignment statistics
#   {SRR}.markdup_metrics.txt  — Picard duplication metrics
#
# INPUT  : 03_trimmed/SRR*_trimmed.fq.gz
# OUTPUT : 04_aligned/
#
# REQUIREMENTS (chipseq_hg38 conda environment):
#   bowtie2 >= 2.5   samtools >= 1.19   picard >= 3.1   bedtools >= 2.31
#
# HOW TO RUN:
#   conda activate chipseq_hg38
#   chmod +x stage4_align.sh
#   ./stage4_align.sh 2>&1 | tee run_stage4.log
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

INDIR="03_trimmed"
OUTDIR="04_aligned"
GENOME_DIR="genome"
SAMPLESHEET="samplesheet.csv"

# Reference paths (created by stage1_download.sh)
BT2_INDEX="${GENOME_DIR}/GRCh38_noalt_as/GRCh38_noalt_as"
BLACKLIST="${GENOME_DIR}/hg38-blacklist.v2.bed"

# Filtering thresholds
MAPQ_THRESH=30          # minimum mapping quality (filters multimappers)

# Low-read samples from Stage 2/3 QC — monitor mapping more carefully
LOW_READ_SAMPLES="SRR2028281 SRR2028286 SRR2028288"
LOW_READ_MAPPED_MIN=2000000   # warn if post-filter reads < 2M

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

# =============================================================================
# HELPER: safe version string for any tool
# FIX: picard returns exit code 1 for --version; || true prevents set -e kill
# =============================================================================
get_version() {
    local tool="$1"
    local ver=""
    case "$tool" in
        picard)
            ver=$(picard --version 2>&1 | head -1) || \
            ver=$(picard MarkDuplicates --version 2>&1 | head -1) || \
            ver="installed (version query not supported)"
            ;;
        *)
            ver=$("$tool" --version 2>&1 | head -1) || \
            ver=$("$tool" version    2>&1 | head -1) || \
            ver="installed"
            ;;
    esac
    echo "$ver"
}

# =============================================================================
# STEP 0 — Preflight
# =============================================================================
header "GSE68952 Stage 4: Alignment to hg38"
info "Date     : $(date)"
info "Threads  : ${THREADS}"
info "Index    : ${BT2_INDEX}"
info "MAPQ     : >= ${MAPQ_THRESH}"
echo ""

info "Checking required tools..."
for tool in bowtie2 samtools picard bedtools; do
    if command -v "$tool" &>/dev/null; then
        VER=$(get_version "$tool")
        success "  $tool : $VER"
    else
        error "$tool not found. Run: conda activate chipseq_hg38"
    fi
done
echo ""

# Check Bowtie2 index
info "Checking Bowtie2 hg38 index..."
if [[ -f "${BT2_INDEX}.1.bt2" ]] || [[ -f "${BT2_INDEX}.1.bt2l" ]]; then
    success "Index found: ${BT2_INDEX}"
else
    error "Bowtie2 index not found at ${BT2_INDEX}\nRun stage1_download.sh to download the index."
fi

# Check blacklist
info "Checking ENCODE blacklist..."
if [[ -f "$BLACKLIST" ]]; then
    NR=$(wc -l < "$BLACKLIST")
    success "Blacklist: ${BLACKLIST} (${NR} regions)"
else
    error "Blacklist not found: ${BLACKLIST}\nRun stage1_download.sh to download it."
fi
echo ""

# Check input trimmed FASTQs
info "Checking trimmed FASTQ inputs..."
MISSING=0
for SRR in "${SRR_ORDERED[@]}"; do
    FQ="${INDIR}/${SRR}_trimmed.fq.gz"
    if [[ -f "$FQ" ]]; then
        SZ=$(du -sh "$FQ" | cut -f1)
        success "  ${SRR}_trimmed.fq.gz  [${SZ}]"
    else
        warn "  MISSING: $FQ"
        MISSING=$((MISSING + 1))
    fi
done
[[ $MISSING -gt 0 ]] && error "$MISSING trimmed FASTQ(s) missing. Run stage3_trim.sh first."
echo ""

# =============================================================================
# STEP 1 — Create output directories
# =============================================================================
header "Step 1: Create directories"
mkdir -p "${OUTDIR}/logs" "${OUTDIR}/flagstats" "${OUTDIR}/metrics"
success "Output directory: ${OUTDIR}/"
echo ""

# =============================================================================
# STEP 2 — Align, filter, mark duplicates (per sample)
# =============================================================================
header "Step 2: Align -> filter -> mark duplicates"

ALIGN_COUNT=0
SKIP_COUNT=0

for SRR in "${SRR_ORDERED[@]}"; do
    NAME="${SAMPLE_NAME[$SRR]}"
    FQ="${INDIR}/${SRR}_trimmed.fq.gz"

    RAW_BAM="${OUTDIR}/${SRR}.raw.bam"
    FILT_BAM="${OUTDIR}/${SRR}.filtered.bam"
    MARKDUP_BAM="${OUTDIR}/${SRR}.markdup.bam"
    MARKDUP_METRICS="${OUTDIR}/metrics/${SRR}.markdup_metrics.txt"
    FLAGSTAT="${OUTDIR}/flagstats/${SRR}.flagstat.txt"
    BT2_LOG="${OUTDIR}/logs/${SRR}.bowtie2.log"

    echo "────────────────────────────────────────────────────────────"
    info "Processing: ${SRR}  |  ${NAME}"

    # Skip if final markdup BAM already exists and is indexed
    if [[ -f "${MARKDUP_BAM}" && -f "${MARKDUP_BAM}.bai" ]]; then
        SZ=$(du -sh "$MARKDUP_BAM" | cut -f1)
        success "Already complete [${SZ}] — skipping: ${MARKDUP_BAM}"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    # ── 2a. Bowtie2 alignment ────────────────────────────────────────────
    # --very-sensitive : -D 20 -R 3 -N 1 -L 20 -i S,1,0.50
    # --no-unal        : suppress unaligned reads (smaller BAM)
    # Pipe directly into samtools sort to avoid large SAM on disk
    info "  [1/4] Bowtie2 alignment (--very-sensitive)..."
    bowtie2 \
        --very-sensitive \
        -x "${BT2_INDEX}" \
        -U "${FQ}" \
        --no-unal \
        -p "${THREADS}" \
        2> "${BT2_LOG}" \
    | samtools sort \
        -@ "${THREADS}" \
        -o "${RAW_BAM}" \
        -
    samtools index -@ "${THREADS}" "${RAW_BAM}"

    # Extract and display alignment rate (grep always succeeds here)
    ALIGN_RATE=$(grep "overall alignment rate" "${BT2_LOG}" | awk '{print $1}') || ALIGN_RATE="?"
    success "  Alignment rate: ${ALIGN_RATE}"

    # ── 2b. Filter: MAPQ >= 30 + remove blacklist regions ────────────────
    # -F 4  : skip unmapped reads
    # -q 30 : keep MAPQ >= 30 only
    # bedtools intersect -v : discard reads overlapping blacklist
    info "  [2/4] Filtering: MAPQ >= ${MAPQ_THRESH} + blacklist removal..."
    samtools view \
        -@ "${THREADS}" \
        -F 4 \
        -q "${MAPQ_THRESH}" \
        -b "${RAW_BAM}" \
    | bedtools intersect \
        -v \
        -abam stdin \
        -b "${BLACKLIST}" \
    | samtools sort \
        -@ "${THREADS}" \
        -o "${FILT_BAM}" \
        -
    samtools index -@ "${THREADS}" "${FILT_BAM}"

    FILT_READS=$(samtools view -c "${FILT_BAM}")
    success "  Filtered reads (MAPQ>=${MAPQ_THRESH}, no blacklist): ${FILT_READS}"

    # ── 2c. Picard MarkDuplicates ─────────────────────────────────────────
    # REMOVE_DUPLICATES false: mark only — MACS2 and DiffBind handle removal
    # OPTICAL_DUPLICATE_PIXEL_DISTANCE 100: correct for Genome Analyzer II
    #   (use 2500 for patterned flowcell instruments like NovaSeq)
    info "  [3/4] Picard MarkDuplicates..."
    picard MarkDuplicates \
        -I "${FILT_BAM}" \
        -O "${MARKDUP_BAM}" \
        -M "${MARKDUP_METRICS}" \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \
        --REMOVE_DUPLICATES false \
        --VALIDATION_STRINGENCY LENIENT \
        --QUIET true \
        2>/dev/null
    samtools index -@ "${THREADS}" "${MARKDUP_BAM}"
    success "  MarkDuplicates complete: ${MARKDUP_BAM}"

    # ── 2d. samtools flagstat ─────────────────────────────────────────────
    info "  [4/4] Generating flagstat..."
    samtools flagstat \
        -@ "${THREADS}" \
        "${MARKDUP_BAM}" \
        > "${FLAGSTAT}"
    success "  Flagstat: ${FLAGSTAT}"

    # ── Low-read check (string match, no boolean under set -e) ───────────
    # FIX: use grep-based string match instead of IS_LOW=false + boolean
    #      comparison, which is fragile under set -e
    if echo "${LOW_READ_SAMPLES}" | grep -qw "${SRR}"; then
        if [[ "${FILT_READS}" =~ ^[0-9]+$ ]] && \
           [[ "${FILT_READS}" -lt "${LOW_READ_MAPPED_MIN}" ]]; then
            warn "  WARNING: ${SRR} has only ${FILT_READS} usable reads"
            warn "  (below ${LOW_READ_MAPPED_MIN} threshold) — flag for sensitivity analysis"
        else
            success "  Low-read sample ${SRR} passes mapped threshold: ${FILT_READS} reads"
        fi
    fi

    # Clean up intermediate files to save disk space
    rm -f "${RAW_BAM}" "${RAW_BAM}.bai"
    info "  Removed intermediate raw BAM."

    ALIGN_COUNT=$((ALIGN_COUNT + 1))
    echo ""
done

# =============================================================================
# STEP 3 — Alignment summary table
# FIX: awk "BEGIN{exit !(condition)}" returns exit 1 when false — kills script
#      under set -e. Use awk with print and compare output string instead.
# =============================================================================
header "Step 3: Alignment summary"
echo ""
printf "%-14s %-8s %-14s %-14s %-10s %-10s\n" \
    "SRR" "Sample" "Total_reads" "Mapped" "%Mapped" "%Dup"
printf '%0.s-' {1..72}; echo ""

WARN_SAMPLES=()

for SRR in "${SRR_ORDERED[@]}"; do
    NAME="${SAMPLE_NAME[$SRR]}"
    FLAGSTAT="${OUTDIR}/flagstats/${SRR}.flagstat.txt"
    METRICS="${OUTDIR}/metrics/${SRR}.markdup_metrics.txt"

    if [[ ! -f "$FLAGSTAT" ]]; then
        printf "%-14s %-8s %-14s\n" "$SRR" "$NAME" "MISSING"
        continue
    fi

    TOTAL=$(grep  "in total"  "$FLAGSTAT" | awk '{print $1}') || TOTAL="?"
    MAPPED=$(grep "mapped ("  "$FLAGSTAT" | head -1 | awk '{print $1}') || MAPPED="?"
    PCT_MAP=$(grep "mapped (" "$FLAGSTAT" | head -1 | grep -oP '[0-9]+\.[0-9]+(?=%)') || PCT_MAP="?"
    PCT_MAP_DISPLAY="${PCT_MAP}%"

    PCT_DUP="?"
    if [[ -f "$METRICS" ]]; then
        # Picard metrics: PCT_DUPLICATION is column 9 in the data row after LIBRARY header
        PCT_DUP_RAW=$(grep -A2 "^LIBRARY" "$METRICS" 2>/dev/null | tail -1 | awk '{print $9}') || true
        if [[ -n "${PCT_DUP_RAW:-}" && "${PCT_DUP_RAW}" != "?" ]]; then
            PCT_DUP=$(awk "BEGIN{printf \"%.1f%%\", ${PCT_DUP_RAW} * 100}") || PCT_DUP="?"
        fi
    fi

    # FIX: use awk to print "LOW" or "" instead of using exit code for comparison
    LOW_FLAG=""
    if [[ "${PCT_MAP}" != "?" ]]; then
        LOW_FLAG=$(awk -v p="${PCT_MAP}" 'BEGIN{print (p+0 < 70) ? "LOW" : ""}')
        if [[ "${LOW_FLAG}" == "LOW" ]]; then
            WARN_SAMPLES+=("${SRR}(${NAME}): ${PCT_MAP_DISPLAY} mapping rate")
        fi
    fi

    STATUS=""
    [[ "${LOW_FLAG}" == "LOW" ]] && STATUS="<< LOW MAPPING"

    printf "%-14s %-8s %-14s %-14s %-10s %-10s %s\n" \
        "$SRR" "$NAME" "${TOTAL}" "${MAPPED}" \
        "${PCT_MAP_DISPLAY}" "${PCT_DUP}" "${STATUS}"
done

echo ""
if [[ ${#WARN_SAMPLES[@]} -gt 0 ]]; then
    warn "Samples with low mapping rate (< 70%):"
    for W in "${WARN_SAMPLES[@]}"; do
        warn "  WARNING: $W"
    done
    echo ""
fi

# =============================================================================
# STEP 4 — Final file verification
# =============================================================================
header "Step 4: Output file verification"
echo ""

MISSING_OUT=0
for SRR in "${SRR_ORDERED[@]}"; do
    MARKDUP_BAM="${OUTDIR}/${SRR}.markdup.bam"
    BAI="${MARKDUP_BAM}.bai"
    if [[ -f "$MARKDUP_BAM" && -f "$BAI" ]]; then
        SZ=$(du -sh "$MARKDUP_BAM" | cut -f1)
        success "  ${SRR}.markdup.bam  [${SZ}]"
    else
        warn "  MISSING: ${SRR}.markdup.bam"
        MISSING_OUT=$((MISSING_OUT + 1))
    fi
done

echo ""
info "============================================================"
info "  Stage 4 Summary"
info "============================================================"
echo "  Aligned this run  : ${ALIGN_COUNT}"
echo "  Already existed   : ${SKIP_COUNT}"
echo "  Missing outputs   : ${MISSING_OUT}"
echo ""
echo "  Output BAMs       : ${OUTDIR}/*.markdup.bam"
echo "  Flagstats         : ${OUTDIR}/flagstats/"
echo "  Picard metrics    : ${OUTDIR}/metrics/"
echo "  Bowtie2 logs      : ${OUTDIR}/logs/"
echo ""

if [[ $MISSING_OUT -eq 0 ]]; then
    success "Stage 4 complete. All 13 BAMs present."
else
    warn "Stage 4 incomplete — re-run to finish missing samples."
fi

echo ""
info "Next step: run stage5_postqc.sh"
info "  Input  : ${OUTDIR}/*.markdup.bam"
info "  Command: bash stage5_postqc.sh"
echo ""
info "$(date): Stage 4 finished."
echo ""

# =============================================================================
# NOTES
# =============================================================================
cat << 'NOTES'
-----------------------------------------------------------------------------
STAGE 4 NOTES — GSE68952 hg38 Alignment
-----------------------------------------------------------------------------
EXPECTED ALIGNMENT RATES (35.5 nt SE, human H3K4me3 ChIP, hg38):
  ChIP samples  : 85-95% overall alignment rate
  Input sample  : 90-97% (genomic DNA maps more uniformly)
  If any sample < 70%: check for contamination or index mismatch

EXPECTED MAPQ>=30 RETENTION:
  Typically 70-85% of aligned reads pass MAPQ >= 30
  Lower GC samples (e.g. HD2 at 54%) may have slightly lower retention

EXPECTED % DUPLICATION (post-filter, 36 nt SE):
  ChIP samples  : 15-50%
  Input sample  : 5-20%
  36 nt SE reads show lower apparent duplication than longer reads
  because fewer reads share identical start coordinates by chance

LOW-READ SAMPLES (C4 / HD2 / HD4):
  If usable reads after filtering < 2M:
    - Still proceed through peak calling
    - DiffBind DESeq2 normalises by library size automatically
    - Run a sensitivity analysis without these samples if key results change

PICARD DUPLICATE MARKING:
  Duplicates are MARKED but NOT removed (REMOVE_DUPLICATES=false)
  MACS2 uses --keep-dup auto (handles marked duplicates)
  DiffBind uses dba.count() which respects duplicate flags per config

OPTICAL DUPLICATE PIXEL DISTANCE:
  Set to 100 for Genome Analyzer II (non-patterned flowcell)
  Would be 2500 for NovaSeq/NextSeq (patterned flowcell)
  Using wrong value slightly over/under-counts optical duplicates

DISK USAGE:
  Raw BAMs deleted after filtering to save space
  markdup BAMs: ~2-8 GB each; 13 samples ~ 50 GB total
  filtered BAMs kept for now — delete with: rm 04_aligned/*.filtered.bam
  after confirming markdup BAMs are correct
-----------------------------------------------------------------------------
NOTES
