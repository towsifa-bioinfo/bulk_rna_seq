#!/usr/bin/env bash
# Download FASTQ files using fastq-dl
# Reads accessions from ../inputs/SRR_Acc_List.txt
# Saves all FASTQ files in the root-level fastq/ directory
# Tools: https://github.com/rpetit3/fastq-dl

# Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778

mkdir -p ../fastq

while read -r ACC; do
    [[ -z "$ACC" ]] && continue
    echo "Downloading $ACC ..."
    fastq-dl \
      --accession "$ACC" \
      --cpus 4 \
      --outdir ../fastq
done < ../inputs/SRR_Acc_List.txt

echo "All FASTQ files are saved in ../fastq/"