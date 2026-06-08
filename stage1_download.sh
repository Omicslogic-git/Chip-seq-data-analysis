#!/usr/bin/env bash
# =============================================================================
# GSE68952 ChIP-seq Pipeline — STAGE 1: Download
# =============================================================================
# Study    : H3K4me3 ChIP-seq in Huntington Disease vs Control human neurons
# Paper    : Dong, Tsuji et al. PLOS ONE 2015 (doi:10.1371/journal.pone.0144398)
# GEO      : GSE68952  |  SRA Study: SRP058390  |  BioProject: PRJNA284153
#
# Design   : 6 HD + 6 Control H3K4me3 ChIP-seq + 1 Input DNA = 13 samples
# Platform : Illumina Genome Analyzer II
# Layout   : Single-end, 36 nt stored spot length
# Genome   : Will align to hg38 (GRCh38)
#
# REQUIREMENTS:
#   - SRA Toolkit >= 3.0.0
#   - wget
#   - ~80 GB free disk space (raw FASTQs + reference genome)
#
# HOW TO RUN:
#   chmod +x stage1_download.sh
#   ./stage1_download.sh 2>&1 | tee run_stage1.log
# =============================================================================

set -euo pipefail

# ── Colour output helpers ──────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; NC='\033[0m'
info()    { echo -e "${CYAN}[INFO]${NC}  $*"; }
success() { echo -e "${GREEN}[OK]${NC}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error()   { echo -e "${RED}[ERROR]${NC} $*"; exit 1; }

# =============================================================================
# CONFIGURATION — edit THREADS to match your machine
# =============================================================================
THREADS=4
OUTDIR="01_raw_fastq"
GENOME_DIR="genome"
SAMPLESHEET="samplesheet.csv"

# =============================================================================
# VERIFIED SRR ACCESSIONS
# =============================================================================
declare -A SAMPLE_META
SAMPLE_META["SRR2028278"]="C1_GSM1688583:Control:M:55:NA"
SAMPLE_META["SRR2028279"]="C2_GSM1688584:Control:M:63:NA"
SAMPLE_META["SRR2028280"]="C3_GSM1688585:Control:M:64:NA"
SAMPLE_META["SRR2028281"]="C4_GSM1688586:Control:M:74:NA"
SAMPLE_META["SRR2028282"]="C5_GSM1688587:Control:M:81:NA"
SAMPLE_META["SRR2028283"]="C6_GSM1688588:Control:F:68:NA"
SAMPLE_META["SRR2028284"]="Input_GSM1688589:Input:F:69:NA"
SAMPLE_META["SRR2028285"]="HD1_GSM1688590:HD:M:55:45"
SAMPLE_META["SRR2028286"]="HD2_GSM1688591:HD:M:71:43"
SAMPLE_META["SRR2028287"]="HD3_GSM1688592:HD:M:43:49"
SAMPLE_META["SRR2028288"]="HD4_GSM1688593:HD:M:68:42"
SAMPLE_META["SRR2028289"]="HD5_GSM1688594:HD:F:56:45"
SAMPLE_META["SRR2028290"]="HD6_GSM1688595:HD:M:69:42"

SRR_ORDERED=(
    SRR2028278 SRR2028279 SRR2028280 SRR2028281 SRR2028282 SRR2028283
    SRR2028284
    SRR2028285 SRR2028286 SRR2028287 SRR2028288 SRR2028289 SRR2028290
)

# =============================================================================
# STEP 0 — Preflight checks
# =============================================================================
info "============================================================"
info "  GSE68952 Stage 1: Download"
info "  $(date)"
info "============================================================"
echo ""

info "Checking required tools..."
for tool in prefetch fasterq-dump vdb-config wget; do
    if command -v "$tool" &>/dev/null; then
        success "$tool found: $(command -v $tool)"
    else
        error "$tool not found. Install SRA Toolkit: https://github.com/ncbi/sra-tools/wiki"
    fi
done
echo ""

info "Verifying SRA Toolkit configuration..."
if vdb-config --query /repository/user/main/public/root &>/dev/null; then
    SRA_CACHE=$(vdb-config --query /repository/user/main/public/root 2>/dev/null | grep -o '".*"' | tr -d '"' || echo "~/.ncbi/public")
    success "SRA cache: $SRA_CACHE"
else
    warn "SRA Toolkit may not be configured. Running vdb-config --interactive is recommended."
    warn "Proceeding anyway — prefetch will use default cache location."
fi
echo ""

# =============================================================================
# STEP 1 — Create directory structure
# =============================================================================
info "Creating directory structure..."
mkdir -p "${OUTDIR}" "${OUTDIR}/sra_cache" "${OUTDIR}/tmp" "${GENOME_DIR}"
success "Directories created."
echo ""

# =============================================================================
# STEP 2 — Write samplesheet
# =============================================================================
info "Writing samplesheet.csv..."
cat > "$SAMPLESHEET" << 'EOF'
SRRAccession,SampleName,Condition,Sex,AgeAtDeath,CAGRepeat,GEOAccession,PMI,Notes
SRR2028278,C1,Control,M,55,NA,GSM1688583,17.0,rostral PFC
SRR2028279,C2,Control,M,63,NA,GSM1688584,24.5,rostral PFC
SRR2028280,C3,Control,M,64,NA,GSM1688585,28.2,rostral PFC
SRR2028281,C4,Control,M,74,NA,GSM1688586,12.0,rostral PFC
SRR2028282,C5,Control,M,81,NA,GSM1688587,8.0,rostral PFC
SRR2028283,C6,Control,F,68,NA,GSM1688588,7.0,rostral PFC
SRR2028284,Input,Input,F,69,NA,GSM1688589,40.2,shared input for all ChIP samples
SRR2028285,HD1,HD,M,55,45,GSM1688590,37.3,BA9 PFC; onset age 44
SRR2028286,HD2,HD,M,71,43,GSM1688591,20.5,BA9 PFC; onset age 52
SRR2028287,HD3,HD,M,43,49,GSM1688592,21.3,BA9 PFC; onset age 28
SRR2028288,HD4,HD,M,68,42,GSM1688593,3.7,BA9 PFC; onset age 45
SRR2028289,HD5,HD,F,56,45,GSM1688594,19.0,BA9 PFC; onset age unknown
SRR2028290,HD6,HD,M,69,42,GSM1688595,19.1,BA9 PFC; onset age 50
EOF
success "Samplesheet written to $SAMPLESHEET"
echo ""

# =============================================================================
# STEP 3 — Download FASTQ files
# FIX: use COUNTER+=1 instead of ((COUNTER++)) to avoid set -e exit on zero
# =============================================================================
info "============================================================"
info "  Starting FASTQ downloads (13 samples)"
info "  Output: ${OUTDIR}/"
info "============================================================"
echo ""

DOWNLOAD_COUNT=0
SKIP_COUNT=0
FAIL_COUNT=0
FAILED_SRRS=()

for SRR in "${SRR_ORDERED[@]}"; do
    META="${SAMPLE_META[$SRR]}"
    SAMPLE_NAME=$(echo "$META" | cut -d: -f1)
    CONDITION=$(echo "$META" | cut -d: -f2)
    OUTFILE="${OUTDIR}/${SRR}.fastq.gz"

    echo "------------------------------------------------------------"
    info "Processing: $SRR  |  Sample: $SAMPLE_NAME  |  Condition: $CONDITION"

    # Skip if already downloaded
    if [[ -f "$OUTFILE" ]]; then
        FILESIZE=$(du -sh "$OUTFILE" | cut -f1)
        success "Already exists ($FILESIZE) — skipping: $OUTFILE"
        SKIP_COUNT=$((SKIP_COUNT + 1))   # ← safe: no ((++)) under set -e
        continue
    fi

    # --- prefetch ---
    info "  Prefetching $SRR from NCBI SRA..."
    if prefetch "$SRR" \
        --output-directory "${OUTDIR}/sra_cache" \
        --max-size 50G \
        2>&1 | grep -v "^$"; then
        success "  Prefetch complete."
    else
        warn "  Prefetch failed for $SRR. Skipping."
        FAILED_SRRS+=("$SRR")
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    # --- fasterq-dump ---
    info "  Converting to FASTQ (${THREADS} threads)..."
    if fasterq-dump "${OUTDIR}/sra_cache/${SRR}" \
        --outdir "$OUTDIR" \
        --outfile "${SRR}.fastq" \
        --threads "$THREADS" \
        --temp "${OUTDIR}/tmp" \
        --skip-technical \
        --force \
        2>&1 | grep -v "^$"; then
        success "  FASTQ conversion complete."
    else
        warn "  fasterq-dump failed for $SRR. Skipping."
        FAILED_SRRS+=("$SRR")
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    # --- gzip ---
    info "  Compressing to .fastq.gz..."
    gzip --fast "${OUTDIR}/${SRR}.fastq"
    FILESIZE=$(du -sh "$OUTFILE" | cut -f1)
    success "  Done: ${OUTFILE} (${FILESIZE})"

    rm -rf "${OUTDIR}/sra_cache/${SRR}"
    DOWNLOAD_COUNT=$((DOWNLOAD_COUNT + 1))
    echo ""
done

# =============================================================================
# STEP 4 — Download hg38 reference resources
# =============================================================================
echo ""
info "============================================================"
info "  Downloading hg38 reference genome resources"
info "============================================================"
echo ""

BOWTIE2_INDEX_DIR="${GENOME_DIR}/GRCh38_noalt_as"
if [[ -d "$BOWTIE2_INDEX_DIR" ]]; then
    success "Bowtie2 hg38 index already exists — skipping."
else
    info "Downloading Bowtie2 hg38 index (~3.5 GB compressed)..."
    info "This is a one-time download. It may take 10-30 minutes."
    wget -q --show-progress \
        -O "${GENOME_DIR}/GRCh38_noalt_as.zip" \
        "https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip" \
        && unzip -q "${GENOME_DIR}/GRCh38_noalt_as.zip" -d "${GENOME_DIR}/" \
        && rm "${GENOME_DIR}/GRCh38_noalt_as.zip"
    success "Bowtie2 index saved to: ${BOWTIE2_INDEX_DIR}/"
fi
echo ""

BLACKLIST="${GENOME_DIR}/hg38-blacklist.v2.bed"
if [[ -f "$BLACKLIST" ]]; then
    success "hg38 blacklist already exists — skipping."
else
    info "Downloading ENCODE hg38 blacklist v2..."
    wget -q --show-progress \
        -O "${BLACKLIST}.gz" \
        "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz" \
        && gunzip "${BLACKLIST}.gz"
    NREGIONS=$(wc -l < "$BLACKLIST")
    success "Blacklist saved: ${BLACKLIST} (${NREGIONS} regions)"
fi
echo ""

CHROMSIZES="${GENOME_DIR}/hg38.chrom.sizes"
if [[ -f "$CHROMSIZES" ]]; then
    success "hg38 chrom.sizes already exists — skipping."
else
    info "Downloading hg38 chromosome sizes..."
    wget -q --show-progress \
        -O "$CHROMSIZES" \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
    success "Chrom sizes saved: ${CHROMSIZES}"
fi
echo ""

# =============================================================================
# STEP 5 — Verify outputs and print summary
# =============================================================================
info "============================================================"
info "  Verification — checking all expected files"
info "============================================================"
echo ""

MISSING=()
PRESENT=()

for SRR in "${SRR_ORDERED[@]}"; do
    OUTFILE="${OUTDIR}/${SRR}.fastq.gz"
    if [[ -f "$OUTFILE" ]]; then
        FILESIZE=$(du -sh "$OUTFILE" | cut -f1)
        success "${SRR}.fastq.gz  [${FILESIZE}]"
        PRESENT+=("$SRR")
    else
        warn "MISSING: ${SRR}.fastq.gz"
        MISSING+=("$SRR")
    fi
done

echo ""
info "============================================================"
info "  STAGE 1 SUMMARY"
info "============================================================"
echo ""
echo "  Downloaded this run : ${DOWNLOAD_COUNT}"
echo "  Already existed     : ${SKIP_COUNT}"
echo "  Failed              : ${FAIL_COUNT}"
echo "  Present (total)     : ${#PRESENT[@]} / 13"
echo "  Missing             : ${#MISSING[@]}"

if [[ ${#MISSING[@]} -gt 0 ]]; then
    warn ""
    warn "  Missing files: ${MISSING[*]}"
    warn "  Re-run this script to retry failed downloads."
fi

echo ""
info "  Reference resources:"
[[ -d "$BOWTIE2_INDEX_DIR" ]] && success "  Bowtie2 hg38 index : ${BOWTIE2_INDEX_DIR}/" \
    || warn "  Bowtie2 hg38 index : MISSING"
[[ -f "$BLACKLIST"         ]] && success "  hg38 blacklist     : ${BLACKLIST}" \
    || warn "  hg38 blacklist     : MISSING"
[[ -f "$CHROMSIZES"        ]] && success "  hg38 chrom.sizes   : ${CHROMSIZES}" \
    || warn "  hg38 chrom.sizes   : MISSING"

echo ""
if [[ ${#MISSING[@]} -eq 0 ]]; then
    success "Stage 1 complete. All 13 FASTQ files present."
    echo ""
    info "Next step: run stage2_qc.sh"
    info "  Input  : ${OUTDIR}/*.fastq.gz"
    info "  Command: bash stage2_qc.sh"
else
    warn "Stage 1 incomplete. Re-run to download missing files."
fi

echo ""
info "  $(date): Stage 1 finished."
echo ""

cat << 'NOTES'
─────────────────────────────────────────────────────────────────────────────
NOTES FOR GSE68952 STAGE 1
─────────────────────────────────────────────────────────────────────────────
DISK SPACE ESTIMATE:
  Raw FASTQs (13 samples × ~3-5 GB each) : ~50 GB
  hg38 Bowtie2 index (compressed)        : ~3.5 GB
  BAM files (Stage 3, all samples)       : ~50 GB
  Total pipeline workspace               : ~150 GB recommended

READ LENGTH NOTE:
  SRA metadata reports AvgSpotLen = 36 nt.
  The paper states 51 nt sequencing runs, but the stored spot length
  is 36 nt. Use --length 20 minimum in Trim Galore; no hard clipping.

SHARED INPUT NOTE:
  SRR2028284 (GSM1688589) is the ONLY input DNA sample.
  Used as control for ALL 12 ChIP-seq samples. Matches paper design.
  DiffBind handles shared input correctly.

SRA CACHE CLEANUP:
  .sra cache files are deleted after each successful conversion.
  Re-run the script at any time to retry failed downloads.
─────────────────────────────────────────────────────────────────────────────
NOTES
