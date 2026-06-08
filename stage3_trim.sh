#!/usr/bin/env bash
# =============================================================================
# GSE68952 ChIP-seq Pipeline — STAGE 3: Adapter Trimming
# =============================================================================
# Study    : H3K4me3 ChIP-seq · Huntington Disease vs Control · hg38
#
# QC FINDINGS FROM STAGE 2 (informs trimming decisions):
#   Per-base quality  : ALL 13 samples PASS — no hard 3′ quality cliff
#   Adapter content   : ALL 13 samples PASS — no significant adapter carryover
#   Read length       : Uniform 36 nt across all samples (confirmed)
#   Duplication       : 4–12% (unusually low, consistent with 36 nt SE reads)
#   Low-read samples  : C4 (2.8M), HD2 (3.2M), HD4 (2.7M) — clean quality,
#                       proceed; revisit after alignment if mapped < 2M
#
# TRIMMING RATIONALE:
#   Even without visible adapter contamination or quality cliffs, Trim Galore
#   is still applied because:
#     (a) Soft 3′ quality trimming (--quality 20) removes marginal bases that
#         would otherwise slightly inflate mismatches in Bowtie2 alignment
#     (b) --length 20 ensures any accidental very-short reads are dropped
#     (c) The --fastqc re-run confirms trimming had the intended effect
#   DO NOT use hard length trimming (--hardtrim5/3) — reads are already 36 nt
#   and further shortening would reduce mappability.
#
# INPUT  : 01_raw_fastq/SRR*.fastq.gz
# OUTPUT : 03_trimmed/SRR*_trimmed.fq.gz      — trimmed FASTQs (→ alignment)
#          03_trimmed/reports/               — Trim Galore trim reports
#          03_trimmed/fastqc_post/           — post-trim FastQC
#          03_trimmed/multiqc_post/          — post-trim MultiQC report
#
# REQUIREMENTS (chipseq_hg38 conda environment):
#   trim-galore >= 0.6.10   cutadapt >= 4.6
#   fastqc >= 0.12          multiqc >= 1.19
#
# HOW TO RUN:
#   conda activate chipseq_hg38
#   chmod +x stage3_trim.sh
#   ./stage3_trim.sh 2>&1 | tee run_stage3.log
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
THREADS=8                    # Cutadapt threads (set to nproc for max speed)
INDIR="01_raw_fastq"
OUTDIR="03_trimmed"
REPORT_DIR="${OUTDIR}/reports"
FASTQC_POST_DIR="${OUTDIR}/fastqc_post"
MULTIQC_POST_DIR="${OUTDIR}/multiqc_post"
SAMPLESHEET="samplesheet.csv"

# ── Trim Galore parameters ────────────────────────────────────────────────────
# --quality 20    : trim 3′ bases with Phred < 20 (sliding window from 3′ end)
# --length 20     : discard reads shorter than 20 nt after trimming
# --dont_gzip     : write uncompressed FASTQ (we re-gzip with pigz for speed)
# --adapter       : not set — auto-detect (Illumina universal for SE data)
# --cores         : Cutadapt internal threads
# NO --hardtrim   : reads are 36 nt; hard trimming would reduce mappability
# NO --clip_R1    : no fixed 5′ adapter to clip (confirmed by FastQC)
TG_QUALITY=20
TG_MIN_LENGTH=20

# Samples ordered: Controls → Input → HD
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

# Low-read samples flagged in QC (monitor after alignment)
LOW_READ_SAMPLES=("SRR2028281" "SRR2028286" "SRR2028288")

# =============================================================================
# STEP 0 — Preflight
# =============================================================================
header "GSE68952 Stage 3: Adapter Trimming"
info "Date     : $(date)"
info "Threads  : ${THREADS}"
info "Quality  : Phred ≥ ${TG_QUALITY}"
info "Min len  : ${TG_MIN_LENGTH} nt"
echo ""

info "Checking required tools..."
for tool in trim_galore cutadapt fastqc multiqc; do
    if command -v "$tool" &>/dev/null; then
        VER=$("$tool" --version 2>&1 | head -1)
        success "  $tool : $VER"
    else
        error "$tool not found. Run: conda activate chipseq_hg38"
    fi
done
echo ""

# Verify pigz available (optional, faster compression)
if command -v pigz &>/dev/null; then
    GZIP_CMD="pigz --fast -p ${THREADS}"
    info "pigz detected — using parallel gzip"
else
    GZIP_CMD="gzip --fast"
    warn "pigz not found — using single-threaded gzip (slower)"
fi
echo ""

info "Checking input FASTQ files..."
MISSING=0
for SRR in "${SRR_ORDERED[@]}"; do
    FQ="${INDIR}/${SRR}.fastq.gz"
    if [[ -f "$FQ" ]]; then
        SZ=$(du -sh "$FQ" | cut -f1)
        success "  ${SRR}  [${SZ}]"
    else
        warn "  MISSING: $FQ"
        MISSING=$((MISSING + 1))
    fi
done
[[ $MISSING -gt 0 ]] && error "$MISSING input FASTQ(s) missing. Run stage1_download.sh first."
echo ""

# Low-read reminder
warn "Low-read samples flagged in QC (will still be trimmed):"
for SRR in "${LOW_READ_SAMPLES[@]}"; do
    warn "  ${SRR} (${SAMPLE_NAME[$SRR]}) — monitor mapped read count after alignment"
done
echo ""

# =============================================================================
# STEP 1 — Create output directories
# =============================================================================
header "Step 1: Create directories"
mkdir -p "${OUTDIR}" "${REPORT_DIR}" "${FASTQC_POST_DIR}" "${MULTIQC_POST_DIR}"
success "Directories created."
echo ""

# =============================================================================
# STEP 2 — Trim Galore (per sample)
# Processes each sample sequentially. On a server with many cores you can
# parallelise across samples using GNU parallel — see note at bottom.
# =============================================================================
header "Step 2: Trim Galore"

TRIM_COUNT=0
SKIP_COUNT=0

for SRR in "${SRR_ORDERED[@]}"; do
    NAME="${SAMPLE_NAME[$SRR]}"
    INFILE="${INDIR}/${SRR}.fastq.gz"
    # Trim Galore output name: <SRR>_trimmed.fq.gz
    OUTFILE="${OUTDIR}/${SRR}_trimmed.fq.gz"
    REPORT_FILE="${REPORT_DIR}/${SRR}_trimming_report.txt"

    echo "------------------------------------------------------------"
    info "Processing: ${SRR}  |  ${NAME}"

    # Skip if already trimmed
    if [[ -f "$OUTFILE" ]]; then
        SZ=$(du -sh "$OUTFILE" | cut -f1)
        success "Already trimmed (${SZ}) — skipping: ${OUTFILE}"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    # ── Run Trim Galore ────────────────────────────────────────────────────
    info "  Running Trim Galore..."
    trim_galore \
        --quality   "${TG_QUALITY}" \
        --length    "${TG_MIN_LENGTH}" \
        --cores     "${THREADS}" \
        --dont_gzip \
        --output_dir "${OUTDIR}" \
        "${INFILE}" \
        2>&1

    # Trim Galore writes: ${SRR}_trimmed.fq  (uncompressed, --dont_gzip)
    # Move trimming report to reports subdirectory
    TG_REPORT=$(ls "${OUTDIR}/${SRR}"*_trimming_report.txt 2>/dev/null | head -1 || true)
    if [[ -n "$TG_REPORT" ]]; then
        mv "$TG_REPORT" "${REPORT_FILE}"
        success "  Trim report: ${REPORT_FILE}"
    fi

    # ── Compress trimmed FASTQ ─────────────────────────────────────────────
    TRIMMED_FQ="${OUTDIR}/${SRR}_trimmed.fq"
    if [[ -f "$TRIMMED_FQ" ]]; then
        info "  Compressing with ${GZIP_CMD%% *}..."
        ${GZIP_CMD} "${TRIMMED_FQ}"
        SZ=$(du -sh "$OUTFILE" | cut -f1)
        success "  Done: ${OUTFILE}  [${SZ}]"
    else
        error "  Expected trimmed FASTQ not found: ${TRIMMED_FQ}"
    fi

    TRIM_COUNT=$((TRIM_COUNT + 1))
    echo ""
done

# =============================================================================
# STEP 3 — Post-trim FastQC
# Re-run FastQC on trimmed FASTQs to confirm quality improvement.
# Key things to check: Per base quality (should be flat/high), read length
# distribution (should peak at 36 nt with a small left tail, no 1-nt reads).
# =============================================================================
header "Step 3: Post-trim FastQC"
info "Running FastQC on trimmed FASTQs..."

FASTQC_TODO=()
for SRR in "${SRR_ORDERED[@]}"; do
    ZIP="${FASTQC_POST_DIR}/${SRR}_trimmed_fastqc.zip"
    TRIMMED="${OUTDIR}/${SRR}_trimmed.fq.gz"
    if [[ -f "$ZIP" ]]; then
        success "  Already done — skipping: ${SRR}_trimmed"
    elif [[ -f "$TRIMMED" ]]; then
        FASTQC_TODO+=("$TRIMMED")
    else
        warn "  Trimmed FASTQ missing — skipping FastQC for ${SRR}"
    fi
done

if [[ ${#FASTQC_TODO[@]} -gt 0 ]]; then
    if command -v parallel &>/dev/null; then
        info "Running ${#FASTQC_TODO[@]} FastQC jobs in parallel..."
        printf '%s\n' "${FASTQC_TODO[@]}" | \
            parallel -j "${THREADS}" \
            fastqc --outdir "${FASTQC_POST_DIR}" --threads 1 --quiet {}
    else
        for FQ in "${FASTQC_TODO[@]}"; do
            fastqc --outdir "${FASTQC_POST_DIR}" --threads 2 --quiet "$FQ"
        done
    fi
    success "Post-trim FastQC complete."
fi
echo ""

# =============================================================================
# STEP 4 — Post-trim MultiQC
# =============================================================================
header "Step 4: Post-trim MultiQC"
info "Generating post-trim MultiQC report..."

multiqc \
    "${OUTDIR}" \
    --outdir "${MULTIQC_POST_DIR}" \
    --filename "GSE68952_multiqc_posttrim" \
    --title "GSE68952 — Post-Trim QC" \
    --comment "Trim Galore: --quality ${TG_QUALITY} --length ${TG_MIN_LENGTH}" \
    --force \
    --quiet

success "Post-trim MultiQC: ${MULTIQC_POST_DIR}/GSE68952_multiqc_posttrim.html"
echo ""

# =============================================================================
# STEP 5 — Trim efficiency summary
# Parse each Trim Galore report to extract key stats
# =============================================================================
header "Step 5: Trimming efficiency summary"
echo ""
printf "%-12s %-8s %-14s %-14s %-14s\n" "SRR" "Sample" "Reads_in" "Reads_out" "%_Removed"
printf '%0.s─' {1..68}; echo ""

for SRR in "${SRR_ORDERED[@]}"; do
    NAME="${SAMPLE_NAME[$SRR]}"
    RFILE="${REPORT_DIR}/${SRR}_trimming_report.txt"

    if [[ ! -f "$RFILE" ]]; then
        printf "%-12s %-8s %-14s\n" "$SRR" "$NAME" "report missing"
        continue
    fi

    READS_IN=$(grep  "Total reads processed"    "$RFILE" | grep -o '[0-9,]*' | tr -d ',' | tail -1)
    READS_OUT=$(grep "Reads written"            "$RFILE" | grep -o '[0-9,]*' | tr -d ',' | tail -1)
    PCT_REMOVED=""
    if [[ -n "$READS_IN" && -n "$READS_OUT" && "$READS_IN" -gt 0 ]]; then
        PCT=$(awk "BEGIN{printf \"%.2f\", (1 - $READS_OUT/$READS_IN)*100}")
        PCT_REMOVED="${PCT}%"
    fi

    printf "%-12s %-8s %-14s %-14s %-14s\n" \
        "$SRR" "$NAME" \
        "${READS_IN:-?}" "${READS_OUT:-?}" "${PCT_REMOVED:-?}"
done
echo ""

# =============================================================================
# STEP 6 — Final summary
# =============================================================================
header "Stage 3 Summary"
echo ""
echo "  Trimmed this run   : ${TRIM_COUNT}"
echo "  Already existed    : ${SKIP_COUNT}"
echo "  Output directory   : ${OUTDIR}/"
echo "  Trim reports       : ${REPORT_DIR}/"
echo "  Post-trim FastQC   : ${FASTQC_POST_DIR}/"
echo "  Post-trim MultiQC  : ${MULTIQC_POST_DIR}/GSE68952_multiqc_posttrim.html"
echo ""

# Final file check
info "Verifying trimmed output files..."
MISSING_OUT=0
for SRR in "${SRR_ORDERED[@]}"; do
    OUTFILE="${OUTDIR}/${SRR}_trimmed.fq.gz"
    if [[ -f "$OUTFILE" ]]; then
        SZ=$(du -sh "$OUTFILE" | cut -f1)
        success "  ${SRR}_trimmed.fq.gz  [${SZ}]"
    else
        warn "  MISSING: ${SRR}_trimmed.fq.gz"
        MISSING_OUT=$((MISSING_OUT + 1))
    fi
done

echo ""
if [[ $MISSING_OUT -eq 0 ]]; then
    success "Stage 3 complete. All 13 trimmed FASTQs present."
else
    warn "Stage 3 incomplete — $MISSING_OUT file(s) missing. Re-run to finish."
fi

echo ""
info "Next step: run stage4_align.sh"
info "  Input  : ${OUTDIR}/*_trimmed.fq.gz"
info "  Command: bash stage4_align.sh"
echo ""
info "$(date): Stage 3 finished."
echo ""

# =============================================================================
# NOTES
# =============================================================================
cat << 'NOTES'
─────────────────────────────────────────────────────────────────────────────
STAGE 3 NOTES — GSE68952 ChIP-seq Trimming
─────────────────────────────────────────────────────────────────────────────
WHAT TO CHECK IN POST-TRIM MultiQC:
  Sequence length distribution : should peak at 36 nt, small left tail ≥ 20 nt
  Per base quality             : flat green across all 36 positions
  % reads removed              : expect 0–5% for clean data like this
  Adapter content              : should be 0% (was already 0% pre-trim)

EXPECTED TRIMMED FILE SIZES:
  Slightly smaller than raw (0–5% reduction for clean data)
  C4, HD2, HD4 will be smallest (low read count from Stage 2)

LOW-READ SAMPLES (C4 / HD2 / HD4):
  Still include in Stage 4. Bowtie2 alignment efficiency is the
  deciding metric — if usable mapped reads (MAPQ ≥ 30, non-dup) < 2M,
  consider flagging for sensitivity analysis but DO NOT remove from
  DiffBind; it uses read depth normalisation and handles imbalance.

PARALLELISING ACROSS SAMPLES (optional, faster):
  If you have ≥ 26 CPU cores, run all samples simultaneously:
    for SRR in SRR2028278 ... SRR2028290; do
        trim_galore --quality 20 --length 20 --cores 2 \
            01_raw_fastq/${SRR}.fastq.gz \
            --output_dir 03_trimmed &
    done
    wait
  Use --cores 2 per job when parallelising across samples.
─────────────────────────────────────────────────────────────────────────────
NOTES
