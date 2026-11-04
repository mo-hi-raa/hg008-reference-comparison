#!/usr/bin/env bash
#SBATCH --job-name=bwa_index_grch38
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/%x_%j.out

# Build the bwa index for the GRCh38 (GIAB v3) reference.
# Submit from the project folder: sbatch scripts/s02_bwa_index_grch38.sh
set -euo pipefail

mkdir -p data/processed/reference/grch38

pixi run --frozen -e pipeline bwa index \
    -p data/processed/reference/grch38/GRCh38_GIABv3 \
    data/raw/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz

echo "done: index files are in data/processed/reference/grch38/"
