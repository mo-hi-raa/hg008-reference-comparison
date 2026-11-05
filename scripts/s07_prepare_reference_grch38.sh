#!/usr/bin/env bash
# Prepare GRCh38 for the callers: plain FASTA, .fai index (samtools)
# and .dict (GATK; Mutect2 needs it).
# Run from the project folder with: sbatch scripts/s07_prepare_reference_grch38.sh
#SBATCH --job-name=prep_ref_grch38
#SBATCH --partition=ws-ia
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail

RUN="pixi run --frozen -e pipeline"
IN=data/raw/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
DIR=data/processed/reference/grch38

zcat $IN > $DIR/GRCh38_GIABv3.fasta
$RUN samtools faidx $DIR/GRCh38_GIABv3.fasta
$RUN gatk CreateSequenceDictionary -R $DIR/GRCh38_GIABv3.fasta -O $DIR/GRCh38_GIABv3.dict
echo "done: FASTA, .fai and .dict in $DIR"
