#!/usr/bin/env bash
# =============================================================================
# GSE68952 ChIP-seq Pipeline — STAGE 2: Quality Control
# =============================================================================
# Study    : H3K4me3 ChIP-seq in Huntington Disease vs Control human neurons
# Genome   : hg38 (GRCh38)
#
# What this stage does:
#   1. Runs FastQC on all 13 raw FASTQ files (parallel)
#   2. Aggregates all FastQC reports into a single MultiQC HTML report
#   3. Prints a per-sample summary table of key QC metrics
#   4. Flags any sample that fails key thresholds
#
# INPUT  : 01_raw_fastq/SRR202828{78..90}.fastq.gz
# OUTPUT : 02_fastqc/           — per-sample FastQC .html + .zip
#          02_fastqc/multiqc/   — aggregated MultiQC report
#
# REQUIREMENTS (all in chipseq_hg38 conda environment):
#   fastqc >= 0.12   multiqc >= 1.19
#
# HOW TO RUN:
#   conda activate chipseq_hg38
#   chmod +x stage2_qc.sh
#   ./stage2_qc.sh 2>&1 | tee run_stage2.log
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
THREADS=8                        # FastQC threads (one thread per sample in parallel)
INDIR="01_raw_fastq"
OUTDIR="02_fastqc"
MULTIQC_DIR="${OUTDIR}/multiqc"
SAMPLESHEET="samplesheet.csv"

# QC thresholds — flag samples that fall outside these
MIN_READS=5000000                # minimum acceptable read count (~5M)
MIN_LENGTH=30                    # minimum acceptable mean read length (nt)
MAX_DUPLICATION=80               # max % duplication before warning

# All 13 SRR accessions in order
SRR_ORDERED=(
    SRR2028278 SRR2028279 SRR2028280 SRR2028281 SRR2028282 SRR2028283
    SRR2028284
    SRR2028285 SRR2028286 SRR2028287 SRR2028288 SRR2028289 SRR2028290
)

# Sample name lookup (matches samplesheet)
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
# STEP 0 — Preflight checks
# =============================================================================
header "GSE68952 Stage 2: Quality Control"
info "Date    : $(date)"
info "Threads : ${THREADS}"
echo ""

info "Checking required tools..."
for tool in fastqc multiqc; do
    if command -v "$tool" &>/dev/null; then
        VER=$("$tool" --version 2>&1 | head -1)
        success "$tool : $VER"
    else
        error "$tool not found. Run: conda activate chipseq_hg38"
    fi
done
echo ""

# Verify input files exist
info "Checking input FASTQ files..."
MISSING_INPUT=0
for SRR in "${SRR_ORDERED[@]}"; do
    FQ="${INDIR}/${SRR}.fastq.gz"
    if [[ -f "$FQ" ]]; then
        SZ=$(du -sh "$FQ" | cut -f1)
        success "  ${SRR}  [${SZ}]"
    else
        warn "  MISSING: $FQ"
        MISSING_INPUT=$((MISSING_INPUT + 1))
    fi
done

if [[ $MISSING_INPUT -gt 0 ]]; then
    error "$MISSING_INPUT input file(s) missing. Run stage1_download.sh first."
fi
echo ""

# =============================================================================
# STEP 1 — Create output directories
# =============================================================================
header "Step 1: Create directories"
mkdir -p "${OUTDIR}" "${MULTIQC_DIR}"
success "Output dirs: ${OUTDIR}/  ${MULTIQC_DIR}/"
echo ""

# =============================================================================
# STEP 2 — Run FastQC on all samples (parallel)
# Uses GNU parallel if available, otherwise sequential fallback.
# FastQC --threads controls threads PER FILE; running multiple files
# simultaneously multiplies throughput on multi-core servers.
# =============================================================================
header "Step 2: FastQC"
info "Running FastQC on 13 samples..."
echo ""

# Check if any FastQC reports are missing (skip already-done samples)
FASTQC_TODO=()
for SRR in "${SRR_ORDERED[@]}"; do
    ZIP="${OUTDIR}/${SRR}_fastqc.zip"
    if [[ -f "$ZIP" ]]; then
        success "  Already done — skipping: ${SRR}"
    else
        FASTQC_TODO+=("${INDIR}/${SRR}.fastq.gz")
    fi
done
echo ""

if [[ ${#FASTQC_TODO[@]} -gt 0 ]]; then
    info "Running FastQC on ${#FASTQC_TODO[@]} sample(s)..."

    if command -v parallel &>/dev/null; then
        # GNU parallel: run up to THREADS FastQC jobs simultaneously
        info "GNU parallel detected — running ${THREADS} jobs in parallel"
        printf '%s\n' "${FASTQC_TODO[@]}" | \
            parallel -j "${THREADS}" \
            fastqc --outdir "${OUTDIR}" --threads 1 --quiet {}
    else
        # Sequential fallback
        warn "GNU parallel not found — running FastQC sequentially (slower)"
        for FQ in "${FASTQC_TODO[@]}"; do
            BASENAME=$(basename "$FQ")
            info "  Processing: $BASENAME"
            fastqc --outdir "${OUTDIR}" --threads 2 --quiet "$FQ"
        done
    fi

    success "FastQC complete."
else
    success "All FastQC reports already present — skipping."
fi
echo ""

# =============================================================================
# STEP 3 — Run MultiQC to aggregate all FastQC reports
# =============================================================================
header "Step 3: MultiQC"
info "Aggregating FastQC reports..."

multiqc \
    "${OUTDIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "GSE68952_multiqc_report" \
    --title "GSE68952 H3K4me3 ChIP-seq — Raw Read QC" \
    --comment "6 HD + 6 Control + 1 Input | hg38 | single-end 36 nt" \
    --force \
    --quiet

success "MultiQC report: ${MULTIQC_DIR}/GSE68952_multiqc_report.html"
echo ""

# =============================================================================
# STEP 4 — Parse FastQC data and print per-sample summary table
# Extracts key metrics from FastQC summary.txt files inside each .zip
# =============================================================================
header "Step 4: Per-sample QC summary"
echo ""

# Print table header
printf "%-12s %-8s %-12s %-10s %-12s %-10s %-8s\n" \
    "SRR" "Sample" "Condition" "Reads" "Length(nt)" "GC%" "Status"
printf '%0.s─' {1..78}; echo ""

FLAGS=()   # collect any flagged samples for the final report

for SRR in "${SRR_ORDERED[@]}"; do
    ZIP="${OUTDIR}/${SRR}_fastqc.zip"
    NAME="${SAMPLE_NAME[$SRR]}"

    # Determine condition from samplesheet
    COND=$(grep "^${SRR}," "${SAMPLESHEET}" 2>/dev/null | cut -d, -f3 || echo "?")

    if [[ ! -f "$ZIP" ]]; then
        printf "%-12s %-8s %-12s %-10s %-12s %-10s %-8s\n" \
            "$SRR" "$NAME" "$COND" "NO_ZIP" "—" "—" "ERROR"
        continue
    fi

    # Extract FastQC basic_statistics from inside the zip (no full unzip needed)
    STATS=$(unzip -p "$ZIP" "*/fastqc_data.txt" 2>/dev/null | \
            awk '/>>Basic Statistics/,/>>END_MODULE/')

    TOTAL_READS=$(echo "$STATS" | grep "Total Sequences" | awk '{print $NF}')
    SEQ_LENGTH=$(echo  "$STATS" | grep "Sequence length"  | awk '{print $NF}')
    GC_PCT=$(echo      "$STATS" | grep "%GC"             | awk '{print $NF}')

    # Extract per-module PASS/WARN/FAIL from summary.txt
    SUMMARY=$(unzip -p "$ZIP" "*/summary.txt" 2>/dev/null)
    PER_BASE_Q=$(echo   "$SUMMARY" | grep "Per base sequence quality" | awk '{print $1}')
    PER_SEQ_Q=$(echo    "$SUMMARY" | grep "Per sequence quality"      | awk '{print $1}')
    ADAPTER=$(echo      "$SUMMARY" | grep "Adapter Content"           | awk '{print $1}')
    DUPLICATION=$(echo  "$SUMMARY" | grep "Sequence Duplication"      | awk '{print $1}')

    # Build status string
    STATUS="OK"
    ISSUES=""
    [[ "$PER_BASE_Q"  == "FAIL" ]] && ISSUES+="LowQ "   && STATUS="WARN"
    [[ "$ADAPTER"     == "FAIL" ]] && ISSUES+="Adapt "  && STATUS="WARN"
    [[ "$DUPLICATION" == "FAIL" ]] && ISSUES+="Dup "    && STATUS="WARN"
    [[ -z "$TOTAL_READS"        ]] && STATUS="PARSE_ERR"

    # Check read count threshold
    if [[ "$TOTAL_READS" =~ ^[0-9]+$ ]] && [[ "$TOTAL_READS" -lt "$MIN_READS" ]]; then
        ISSUES+="LowN "
        STATUS="WARN"
    fi

    [[ "$STATUS" != "OK" ]] && FLAGS+=("${SRR}(${NAME}): ${ISSUES}")

    printf "%-12s %-8s %-12s %-10s %-12s %-10s %-8s\n" \
        "$SRR" "$NAME" "$COND" \
        "${TOTAL_READS:-?}" \
        "${SEQ_LENGTH:-?}" \
        "${GC_PCT:-?}" \
        "$STATUS"
done

echo ""

# =============================================================================
# STEP 5 — Key QC expectations for H3K4me3 ChIP-seq
# =============================================================================
header "Step 5: What to look for in the MultiQC report"

cat << 'NOTES'
  Open: 02_fastqc/multiqc/GSE68952_multiqc_report.html

  ┌─────────────────────────────────────────────────────────────────────────┐
  │  METRIC               EXPECTED FOR H3K4me3 ChIP-seq                    │
  ├─────────────────────────────────────────────────────────────────────────┤
  │  Total reads          ≥ 10M per sample (ENCODE minimum for H3K4me3)    │
  │  Read length          36 nt (SRA stored); paper reports 51 nt          │
  │  Per base quality     All green (Phred ≥ 28 across all positions)      │
  │  GC content           45–60% expected for human ChIP-seq               │
  │  Adapter content      Should be low/absent for 36 nt reads             │
  │  Sequence duplication High duplication is NORMAL for ChIP-seq          │
  │                       (enrichment = many reads from same fragment)     │
  │                       Do NOT fail samples based on duplication alone   │
  │  N content            Should be < 5% at any position                   │
  └─────────────────────────────────────────────────────────────────────────┘

  DUPLICATION NOTE (critical for ChIP-seq):
  FastQC will almost certainly flag "Sequence Duplication Levels" as WARN
  or FAIL for the ChIP samples — this is EXPECTED and NOT a problem.
  ChIP-seq intentionally enriches for specific genomic regions, producing
  many identical reads from the same fragments. Duplication levels of
  40–80% are normal. The Input sample may show lower duplication.

  WHAT WOULD ACTUALLY FAIL:
    - Per base quality drop in last few bases → trim with Trim Galore (Stage 3)
    - Adapter contamination → auto-handled by Trim Galore (Stage 3)
    - < 5M reads in any sample → flag for exclusion consideration
    - GC content bimodal peaks → possible contamination, investigate

NOTES

# =============================================================================
# STEP 6 — Final summary
# =============================================================================
header "Stage 2 Summary"
echo ""
echo "  FastQC reports : ${OUTDIR}/"
echo "  MultiQC report : ${MULTIQC_DIR}/GSE68952_multiqc_report.html"
echo ""

if [[ ${#FLAGS[@]} -gt 0 ]]; then
    warn "  Flagged samples (review in MultiQC before proceeding):"
    for F in "${FLAGS[@]}"; do
        warn "    ⚠  $F"
    done
    echo ""
    warn "  Review flagged samples in the MultiQC report."
    warn "  Duplication warnings are expected — see notes above."
else
    success "  All samples passed QC thresholds."
fi

echo ""
success "Stage 2 complete."
echo ""
info "Next step: run stage3_trim.sh"
info "  Input  : ${INDIR}/*.fastq.gz"
info "  Command: bash stage3_trim.sh"
echo ""
info "$(date): Stage 2 finished."
echo ""
