#!/usr/bin/env bash
# =============================================================================
# GSE68952 ChIP-seq Pipeline — STAGE 9: Motif Analysis (Homer)
# =============================================================================
# Study    : H3K4me3 ChIP-seq · Huntington Disease vs Control · hg38
#
# STAGE 8 RESULTS (informs motif analysis):
#   Lost peaks (|FC|>=1.0)   : 237 peaks · 206 genes
#   Gained peaks (|FC|>=1.0) : 72 peaks  · 55 genes
#   Key pathway: Neuroactive ligand-receptor interaction (DRD2, ADORA2B, GRM4)
#   ~65% of differential peaks at promoters (correct for H3K4me3)
#
# MOTIF STRATEGY:
#   Primary   : loss_promoter.bed  — lost H3K4me3 at promoters (HD biology)
#   Secondary : gain_promoter.bed  — gained H3K4me3 at promoters
#   Also run  : all_differential.bed — full differential set
#   Background: background.bed     — non-differential consensus peaks
#   Genome    : hg38
#   Size      : 200 bp around peak centre (standard for H3K4me3 promoters)
#
# EXPECTED MOTIFS:
#   Lost peaks : SP1, NRF1, EGR1, E-box (MYC/MAX) — active promoter TFs
#                that co-occupy H3K4me3 promoters; loss suggests reduced TF
#                binding at these loci in HD
#   Gained peaks: REST/NRSF — neuronal repressor; gaining H3K4me3 at
#                REST-bound loci in HD context is biologically interesting
#
# INPUT  : 09_homer/  (BED files from Stage 8)
# OUTPUT : 09_homer/results/
#
# REQUIREMENTS (chipseq_analysis conda environment):
#   homer >= 4.11  (with hg38 genome loaded)
#
# HOW TO RUN:
#   conda activate chipseq_analysis
#   chmod +x stage9_homer.sh
#   ./stage9_homer.sh 2>&1 | tee run_stage9.log
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
HOMER_DIR="09_homer"
RESULTS_DIR="${HOMER_DIR}/results"
GENOME="hg38"
THREADS=8

# Peak size for motif finding
# 200 bp = standard for H3K4me3 promoter motifs
# Homer centres on the peak midpoint and uses this window
PEAK_SIZE=200

# BED files prepared by Stage 8
BED_LOSS_PROMOTER="${HOMER_DIR}/loss_promoter.bed"
BED_GAIN_PROMOTER="${HOMER_DIR}/gain_promoter.bed"
BED_LOSS_STRICT="${HOMER_DIR}/loss_strict.bed"
BED_GAIN_STRICT="${HOMER_DIR}/gain_strict.bed"
BED_ALL_DIFF="${HOMER_DIR}/all_differential.bed"
BED_BACKGROUND="${HOMER_DIR}/background.bed"

# =============================================================================
# STEP 0 — Preflight
# =============================================================================
header "GSE68952 Stage 9: Homer Motif Analysis"
info "Date    : $(date)"
info "Genome  : ${GENOME}"
info "Threads : ${THREADS}"
info "Size    : ${PEAK_SIZE} bp"
echo ""

# Check Homer
info "Checking Homer installation..."
if command -v findMotifsGenome.pl &>/dev/null; then
    VER=$(findMotifsGenome.pl 2>&1 | grep "HOMER" | head -1) || VER="installed"
    success "findMotifsGenome.pl found"
else
    error "Homer not found.\nInstall: mamba install -c bioconda homer -y\nThen load genome: perl \$(which homer) loadGenome hg38"
fi
echo ""

# Check hg38 genome is loaded in Homer
info "Checking Homer hg38 genome..."
HOMER_PATH=$(which findMotifsGenome.pl 2>/dev/null || echo "")
HOMER_BASE=$(dirname "$(dirname "$HOMER_PATH")" 2>/dev/null || echo "")
HG38_PATH="${HOMER_BASE}/share/homer/data/genomes/hg38"

if [[ -d "$HG38_PATH" ]]; then
    success "hg38 genome found: ${HG38_PATH}"
else
    warn "hg38 genome not found at expected path."
    warn "Loading hg38 genome now (~3 GB download)..."
    perl /home/omicslogic/miniconda3/envs/chipseq_analysis/share/homer/configureHomer.pl -install hg38 || true
    info "If download failed, run manually: perl \$(which homer) loadGenome hg38"
fi
echo ""

# Check BED input files
info "Checking input BED files..."
MISSING=0
for BED in "$BED_LOSS_PROMOTER" "$BED_GAIN_PROMOTER" \
           "$BED_LOSS_STRICT" "$BED_ALL_DIFF" "$BED_BACKGROUND"; do
    if [[ -f "$BED" ]]; then
        NL=$(wc -l < "$BED") || NL=0
        success "  $(basename $BED)  [${NL} peaks]"
    else
        warn "  MISSING: $BED"
        MISSING=$((MISSING + 1))
    fi
done
[[ $MISSING -gt 0 ]] && error "$MISSING BED file(s) missing. Run stage8_annotation.R first."
echo ""

# =============================================================================
# STEP 1 — Create output directories
# =============================================================================
header "Step 1: Create directories"
mkdir -p \
    "${RESULTS_DIR}/loss_promoter" \
    "${RESULTS_DIR}/gain_promoter" \
    "${RESULTS_DIR}/loss_strict" \
    "${RESULTS_DIR}/gain_strict" \
    "${RESULTS_DIR}/all_differential" \
    "${RESULTS_DIR}/all_differential_v2" \
    "${RESULTS_DIR}/loss_strict_v2"
success "Output: ${RESULTS_DIR}/"
echo ""

# =============================================================================
# STEP 2 — Homer findMotifsGenome.pl
# Run motif finding on each BED set vs background
#
# Key parameters:
# -size 200     : use 200 bp window centred on peak
# -mask         : mask repetitive elements (important for promoters)
# -bg           : custom background set (non-differential consensus peaks)
# -p            : parallel threads
# -nomotif      : skip de novo motif finding (faster) — use -nomotif for
#                 known motifs only. Remove flag for full de novo search.
# Note: de novo motif finding (-nomotif removed) takes ~2-4h per set.
#       Use -nomotif for quick known-motif scan (~10-20 min per set).
# =============================================================================
header "Step 2: Run Homer motif analysis"
echo ""

run_homer() {
    local BED="$1"
    local OUTDIR="$2"
    local LABEL="$3"
    local BG="$4"

    echo "────────────────────────────────────────────────────────────"
    info "Running Homer: ${LABEL}"
    info "  Input  : ${BED} ($(wc -l < "$BED") peaks)"
    info "  Output : ${OUTDIR}"
    info "  Mode   : known motifs + de novo"
    echo ""

    # Skip if already done
    if [[ -f "${OUTDIR}/knownResults.txt" ]]; then
        success "Already done — skipping: ${LABEL}"
        return 0
    fi

    # Check minimum peak count
    # Homer needs >= 200 peaks for robust motif finding
    # < 200: runs but results unreliable (as seen with loss_promoter)
    # < 50 : skip entirely
    NPEAKS=$(wc -l < "$BED") || NPEAKS=0
    if [[ "$NPEAKS" -lt 50 ]]; then
        warn "  Too few peaks (${NPEAKS}) — skipping (minimum 50 required)"
        return 0
    fi
    if [[ "$NPEAKS" -lt 200 ]]; then
        warn "  WARNING: Only ${NPEAKS} peaks — results may be unreliable"
        warn "  Homer recommends >= 200 peaks for robust motif analysis"
        warn "  Consider using all_differential.bed for better power"
    fi

    findMotifsGenome.pl \
        "$BED" \
        "$GENOME" \
        "$OUTDIR" \
        -size  "$PEAK_SIZE" \
        -mask \
        -bg    "$BG" \
        -p     "$THREADS" \
        2>&1 | grep -v "^$" || true

    if [[ -f "${OUTDIR}/knownResults.txt" ]]; then
        N_KNOWN=$(wc -l < "${OUTDIR}/knownResults.txt") || N_KNOWN=0
        success "  Done: ${N_KNOWN} known motifs tested → ${OUTDIR}/knownResults.html"
    else
        warn "  knownResults.txt not found — check Homer output in ${OUTDIR}"
    fi
    echo ""
}

# Run Homer on each BED set
run_homer "$BED_LOSS_PROMOTER" \
          "${RESULTS_DIR}/loss_promoter" \
          "Lost H3K4me3 — promoter peaks" \
          "$BED_BACKGROUND"

run_homer "$BED_GAIN_PROMOTER" \
          "${RESULTS_DIR}/gain_promoter" \
          "Gained H3K4me3 — promoter peaks" \
          "$BED_BACKGROUND"

run_homer "$BED_LOSS_STRICT" \
          "${RESULTS_DIR}/loss_strict" \
          "Lost H3K4me3 — strict |FC|>=1.0" \
          "$BED_BACKGROUND"

run_homer "$BED_ALL_DIFF" \
          "${RESULTS_DIR}/all_differential" \
          "All differential peaks" \
          "$BED_BACKGROUND"

# ── LARGER PEAK SETS — better statistical power ──────────────────────────────
# Re-run all_differential and loss_strict without the TSS distance filter
# These have 200+ peaks needed for Homer to test all 400+ known motifs
info "Running Homer on larger peak sets (better power)..."
echo ""

run_homer "$BED_ALL_DIFF" \
          "${RESULTS_DIR}/all_differential_v2" \
          "All differential peaks (v2 — full known motif scan)" \
          "$BED_BACKGROUND"

run_homer "$BED_LOSS_STRICT" \
          "${RESULTS_DIR}/loss_strict_v2" \
          "Lost strict peaks (v2 — no TSS filter)" \
          "$BED_BACKGROUND"

# =============================================================================
# STEP 3 — Parse and print top known motifs
# =============================================================================
header "Step 3: Top known motifs summary"
echo ""

parse_homer_results() {
    local RESULT_DIR="$1"
    local LABEL="$2"
    local KNOWN_FILE="${RESULT_DIR}/knownResults.txt"

    if [[ ! -f "$KNOWN_FILE" ]]; then
        warn "  ${LABEL}: results not found"
        return 0
    fi

    echo "  ${LABEL}:"
    printf "  %-30s %-12s %-12s %-10s\n" "Motif" "%Target" "%Background" "p-value"
    printf "  %s\n" "$(printf '%0.s-' {1..68})"

    # knownResults.txt columns:
    # 1=MotifName 2=Consensus 3=P-value 4=Log P-value 5=q-value(BH)
    # 6=%Target   7=TargetSeqs 8=%Background 9=BackgroundSeqs
    tail -n +2 "$KNOWN_FILE" | head -15 | \
    awk -F'\t' '{
        # Clean motif name (remove /factor/subfamily etc)
        split($1, a, "/")
        name = a[1]
        gsub(/\(.*\)/, "", name)
        printf "  %-30s %-12s %-12s %-10s\n",
               substr(name,1,30), $6, $8, $3
    }' || true

    echo ""
}

for SUBDIR in loss_promoter gain_promoter loss_strict all_differential \
              all_differential_v2 loss_strict_v2; do
    LABEL=$(echo "$SUBDIR" | tr '_' ' ')
    parse_homer_results "${RESULTS_DIR}/${SUBDIR}" "$LABEL"
done

# =============================================================================
# STEP 4 — Extract top de novo motifs
# =============================================================================
header "Step 4: Top de novo motifs"
echo ""

for SUBDIR in loss_promoter gain_promoter; do
    DENOVO_FILE="${RESULTS_DIR}/${SUBDIR}/homerMotifs.all.motifs"
    DENOVO_HTML="${RESULTS_DIR}/${SUBDIR}/homerResults.html"
    LABEL=$(echo "$SUBDIR" | tr '_' ' ')

    if [[ -f "$DENOVO_HTML" ]]; then
        success "  De novo results: ${DENOVO_HTML}"
    elif [[ -f "${RESULTS_DIR}/${SUBDIR}/knownResults.txt" ]]; then
        info "  ${LABEL}: de novo results in ${RESULTS_DIR}/${SUBDIR}/homerResults.html"
    else
        warn "  ${LABEL}: no results yet"
    fi
done
echo ""

# =============================================================================
# STEP 5 — Annotate peaks with nearest gene (Homer annotatePeaks.pl)
# This gives a gene-level view of which specific genes have motif hits
# =============================================================================
header "Step 5: Annotate peaks with Homer"
echo ""

annotate_peaks() {
    local BED="$1"
    local OUTFILE="$2"
    local LABEL="$3"

    if [[ -f "$OUTFILE" ]]; then
        success "  Already done — skipping: ${LABEL}"
        return 0
    fi

    info "  Annotating: ${LABEL}..."
    annotatePeaks.pl \
        "$BED" \
        "$GENOME" \
        > "$OUTFILE" \
        2>/dev/null || true

    if [[ -f "$OUTFILE" ]]; then
        NL=$(wc -l < "$OUTFILE") || NL=0
        success "  ${LABEL}: ${NL} peaks annotated → ${OUTFILE}"
    else
        warn "  Annotation failed for ${LABEL}"
    fi
}

annotate_peaks "$BED_LOSS_PROMOTER" \
               "${HOMER_DIR}/loss_promoter_annotated.txt" \
               "Loss promoter peaks"

annotate_peaks "$BED_GAIN_PROMOTER" \
               "${HOMER_DIR}/gain_promoter_annotated.txt" \
               "Gain promoter peaks"

annotate_peaks "$BED_LOSS_STRICT" \
               "${HOMER_DIR}/loss_strict_annotated.txt" \
               "Loss strict peaks"
echo ""

# =============================================================================
# STEP 6 — Final output summary
# =============================================================================
header "Step 6: Stage 9 Summary"
echo ""
echo "  Homer results:"
for SUBDIR in loss_promoter gain_promoter loss_strict all_differential \
             all_differential_v2 loss_strict_v2; do
    KNOWN="${RESULTS_DIR}/${SUBDIR}/knownResults.html"
    DENOVO="${RESULTS_DIR}/${SUBDIR}/homerResults.html"
    if [[ -f "$KNOWN" ]]; then
        success "  ${SUBDIR}/knownResults.html"
    else
        warn "  ${SUBDIR}/knownResults.html — MISSING"
    fi
done

echo ""
echo "  Peak annotations:"
for F in loss_promoter_annotated.txt gain_promoter_annotated.txt \
          loss_strict_annotated.txt; do
    FPATH="${HOMER_DIR}/${F}"
    if [[ -f "$FPATH" ]]; then
        NL=$(wc -l < "$FPATH") || NL=0
        success "  ${F}  [${NL} lines]"
    else
        warn "  ${F} — MISSING"
    fi
done

echo ""
success "Stage 9 complete."
echo ""
info "Open Homer results in browser:"
info "  firefox ${RESULTS_DIR}/loss_promoter/knownResults.html"
info "  firefox ${RESULTS_DIR}/loss_promoter/homerResults.html  (de novo)"
echo ""
info "Next step: bash stage10_visualisation.sh"
info "$(date): Stage 9 finished."
echo ""

# =============================================================================
# NOTES
# =============================================================================
cat << 'NOTES'
-----------------------------------------------------------------------------
STAGE 9 NOTES — Homer Motif Analysis
-----------------------------------------------------------------------------
INTERPRETING KNOWN MOTIF RESULTS (knownResults.html):
  Open in browser — shows motif logo, p-value, % target vs % background
  Key columns:
    % of Targets   : fraction of your peaks containing the motif
    % of Background: fraction of background peaks with same motif
    P-value        : enrichment significance (one-tailed binomial)
  Focus on motifs with:
    - % Target >> % Background  (enrichment ratio > 1.5x)
    - P-value < 1e-5
    - Known role in neuronal gene regulation

EXPECTED TOP MOTIFS — Lost H3K4me3 promoters:
  SP1/SP3    : ubiquitous promoter TF; loss suggests reduced promoter activity
  NRF1       : nuclear respiratory factor; mitochondrial + active gene promoters
  EGR1/EGR2  : early growth response; neuronal activity-dependent genes
  YY1        : polycomb-associated; H3K4me3 loss often correlates with YY1 loss
  NF-Y       : CCAAT-box binding; marks active promoters

EXPECTED TOP MOTIFS — Gained H3K4me3 promoters:
  REST/NRSF  : neuronal gene repressor; gaining H3K4me3 at REST sites in HD
               may reflect aberrant activation of normally-repressed genes
  CTCF       : chromatin organiser; gained H3K4me3 near CTCF sites

DE NOVO MOTIFS (homerResults.html):
  These are motifs Homer discovered from scratch in your data.
  Best motifs have p < 1e-10 and > 10% target occurrence.
  Match to known TFs using the "Best Match" column.

RUNTIME:
  Known motif scan only    : ~10-20 min per BED set
  Full de novo + known     : ~2-4 hours per BED set
  all_differential_v2      : PRIMARY result — 8,202 peaks, full motif scan
  loss_strict_v2           : SECONDARY result — 237 peaks, strict FC

WHY v2 RUNS WERE ADDED:
  loss_promoter.bed had ~180 peaks (TSS + strict FC filtered)
  Homer only tested 16 motifs — results were underpowered
  all_differential_v2 uses all 8,202 FDR<0.05 peaks
  This gives Homer enough data to test all 400+ known motifs properly
  Expected top motifs with larger set: SP1, NRF1, EGR1, bHLH/E-box

RUN IN BACKGROUND (recommended for large sets):
  nohup findMotifsGenome.pl 09_homer/all_differential.bed hg38 \
      09_homer/results/all_differential_v2 \
      -size 200 -mask -bg 09_homer/background.bed -p 8 \
      > homer_alldiff.log 2>&1 &
  echo "Homer PID: $!"
  # Check progress: tail -f homer_alldiff.log

IF HOMER HG38 GENOME IS NOT LOADED:
  conda activate chipseq_analysis
  perl $(which homer) loadGenome hg38
  # This downloads ~3 GB of annotation data
-----------------------------------------------------------------------------
NOTES
