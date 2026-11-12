#!/usr/bin/env bash
# Full DSA run: sort pieces, merge, mark duplicates, index for one sample and one reference.
# Usage: sbatch scripts/s31_finish_dsa.sh HG008-T combined_mq60
#SBATCH --job-name=finish_dsa
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
SM=$1; LEVEL=$2
RUN="pixi run --frozen -e pipeline"
A=data/processed/alignments/dsa_$LEVEL
mkdir -p $A/sorted_pieces $A/tmp
for bam in $A/unsorted/${SM}_*.bam; do
  name=$(basename $bam .bam)
  [ -s $A/sorted_pieces/$name.sorted.bam ] || $RUN samtools sort -@ 16 -m 2G -T $A/tmp/$name -o $A/sorted_pieces/$name.sorted.bam $bam
  echo "sorted: $name"
done
$RUN samtools merge -@ 16 -u -o - $A/sorted_pieces/${SM}_*.sorted.bam | $RUN samtools markdup -@ 16 -f $A/${SM}.markdup_stats.txt - $A/${SM}.dsa_$LEVEL.bam
$RUN samtools index -@ 16 $A/${SM}.dsa_$LEVEL.bam
$RUN samtools flagstat -@ 16 $A/${SM}.dsa_$LEVEL.bam > $A/${SM}.dsa_$LEVEL.flagstat.txt
echo "done: $A/${SM}.dsa_$LEVEL.bam"
