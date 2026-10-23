#!/usr/bin/env bash
set -euo pipefail

# Run quality check with FastQC and MultiQC
# Input: FASTQ files in ../fastq/
# Output: FastQC reports in ../outputs/qc/fastqc/, MultiQC report in ../outputs/qc/

# Resolve paths relative to this script, so it works from any working directory
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

FASTQ_DIR="${PROJECT_DIR}/fastq"
QC_DIR="${PROJECT_DIR}/outputs/qc"
FASTQC_DIR="${QC_DIR}/fastqc"
THREADS=${THREADS:-8}

mkdir -p "$FASTQC_DIR"

# Collect FASTQ files up front so a missing/empty directory fails loudly
shopt -s nullglob
fastq_files=("${FASTQ_DIR}"/*.fastq.gz)
shopt -u nullglob

if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo "ERROR: no *.fastq.gz files found in ${FASTQ_DIR}" >&2
    exit 1
fi

echo "Running FastQC on ${#fastq_files[@]} files with ${THREADS} threads..."
fastqc --threads "$THREADS" -o "$FASTQC_DIR" "${fastq_files[@]}"

echo "Running MultiQC..."
multiqc --force "$FASTQC_DIR" -o "$QC_DIR"

echo "Quality check complete! Reports saved in ${QC_DIR}"
