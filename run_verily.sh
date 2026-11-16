#!/bin/bash
# ASAP snRNA-seq Pipeline — Verily Workbench (Google Batch)
#
# Prerequisites:
#   1. wb auth login && gcloud auth application-default login
#   2. wb workspace set --id=asap-crn-reference-workspace-v5-copy-5986
#   3. Samplesheet uploaded to ws_files (gs://...) — see SAMPLESHEET below
#   4. 1000G VCF uploaded to ws_files refs/ — see CELLSNP_VCF below
#   5. nextflow_verily.config: fill in CHANGE_ME values
#
# Usage:
#   bash run_verily.sh                    # SoupX + Harmony, all samples
#   bash run_verily.sh -resume            # resume after failure
#   bash run_verily.sh --skip-ambient     # skip ambient correction
#   bash run_verily.sh --scvi             # scVI integration (GPU)
#   bash run_verily.sh --males            # males only
#   bash run_verily.sh --females          # females only

set -e

# ── Fill these in ──────────────────────────────────────────────────────────────
# Get the ws_files bucket root with:  wb resource resolve ws_files
WS_FILES_BUCKET="gs://cloned-ws-files-wb-lambent-papaya-979"

SAMPLESHEET="${WS_FILES_BUCKET}/samplesheets/ASAP_samplesheet.xlsx"
OUTDIR="${WS_FILES_BUCKET}/results_ASAP_verily"
WORKDIR="${WS_FILES_BUCKET}/nextflow_work"
CELLSNP_VCF="${WS_FILES_BUCKET}/refs/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
# ──────────────────────────────────────────────────────────────────────────────

# ── Defaults ──────────────────────────────────────────────────────────────────
AMBIENT_TOOL="soupx"
INTEGRATION_TOOL="harmony"
SEX_FILTER="all"
EXTRA_FILTER=""
RESUME=""

# ── Parse flags ───────────────────────────────────────────────────────────────
for arg in "$@"; do
    case "$arg" in
        --skip-ambient)  AMBIENT_TOOL="skip"      ;;
        --scvi)          INTEGRATION_TOOL="scvi"  ;;
        --males)         SEX_FILTER="M"           ;;
        --females)       SEX_FILTER="F"           ;;
        -resume)         RESUME="-resume"         ;;
    esac
done

echo "ambient_tool     : ${AMBIENT_TOOL}"
echo "integration_tool : ${INTEGRATION_TOOL}"
echo "sex_filter       : ${SEX_FILTER}"
echo "extra_filter     : ${EXTRA_FILTER:-'(none)'}"
echo "workdir          : ${WORKDIR}"
echo ""

wb nextflow run "${PIPELINE_DIR}/main.nf" \
    --samplesheet      "${SAMPLESHEET}" \
    --outdir           "${OUTDIR}" \
    --ambient_tool     "${AMBIENT_TOOL}" \
    --integration_tool "${INTEGRATION_TOOL}" \
    --sex_filter       "${SEX_FILTER}" \
    --extra_filter     "${EXTRA_FILTER}" \
    --cellsnp_vcf      "${CELLSNP_VCF}" \
    -c "${PIPELINE_DIR}/nextflow_verily.config" \
    -w "${WORKDIR}" \
    ${RESUME}
