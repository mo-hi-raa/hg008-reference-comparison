#!/usr/bin/env bash
# Prepare T2T-CHM13v2.0 for the callers: plain FASTA, .fai index (samtools)
# and .dict (GATK; Mutect2 needs it).
# Run from the project folder with: sbatch scripts/s13_prepare_reference_grch38.sh
#SBATCH --job-name=prep_ref_chm13v2
#SBATCH --partition=ws-ia
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail

RUN="pixi run --frozen -e pipeline"
IN=data/raw/chm13v2.0_maskedY_rCRS.fa.gz
DIR=data/processed/reference/chm13v2

zcat $IN > $DIR/CHM13v2_maskedY_rCRS.fasta
$RUN samtools faidx $DIR/CHM13v2_maskedY_rCRS.fasta
$RUN gatk CreateSequenceDictionary -R $DIR/CHM13v2_maskedY_rCRS.fasta -O $DIR/CHM13v2_maskedY_rCRS.dict
echo "done: FASTA, .fai and .dict in $DIR"
