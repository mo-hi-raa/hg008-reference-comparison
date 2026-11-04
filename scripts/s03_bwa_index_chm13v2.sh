#!/usr/bin/env bash
#SBATCH --job-name=bwa_index_chm13v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/%x_%j.out

# Build the bwa index for the T2T-CHM13v2.0 (maskedY, rCRS) reference.
# Submit from the project folder: sbatch scripts/s03_bwa_index_chm13v2.sh
set -euo pipefail

mkdir -p data/processed/reference/chm13v2

pixi run --frozen -e pipeline bwa index \
    -p data/processed/reference/chm13v2/CHM13v2_maskedY_rCRS \
    data/raw/chm13v2.0_maskedY_rCRS.fa.gz

echo "done: index files are in data/processed/reference/chm13v2/"
