#!/usr/bin/env bash
# DSA chr10 pilot: extract chr10 read pairs for one sample from the T2T (CHM13) BAM.
# Usage: sbatch scripts/s25a_pilot_dsa_reads.sh HG008-T
#SBATCH --job-name=pilot_reads
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
SM=$1
RUN="pixi run --frozen -e pipeline"
P=data/processed/dsa_pilot_chr10
mkdir -p $P/reads $P/bams $P/tmp
$RUN samtools view -u -F 0x900 data/processed/alignments/chm13v2/$SM.chm13v2.bam chr10 \
  | $RUN samtools collate -u -O -T $P/tmp/$SM - \
  | $RUN samtools fastq -n -1 $P/reads/${SM}_R1.fastq.gz -2 $P/reads/${SM}_R2.fastq.gz -0 /dev/null -s /dev/null -
echo "read pairs: $(( $(zcat $P/reads/${SM}_R1.fastq.gz | wc -l) / 4 ))"
echo "done: $SM"
