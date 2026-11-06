#!/usr/bin/env bash
# Finish the GRCh38 alignment for ONE sample: sort each piece by genome
# position, merge the pieces, mark duplicates, index.
# Usage: sbatch scripts/s09_finish_bam_grch38.sh HG008-T
#SBATCH --job-name=finish_bam_grch38
#SBATCH --partition=ws-ia
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=23:59:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail

SM=$1
RUN="pixi run --frozen -e pipeline"
IN=data/processed/alignments/grch38/unsorted
SORTED=data/processed/alignments/grch38/sorted_pieces
OUT=data/processed/alignments/grch38
T=$SLURM_CPUS_PER_TASK
mkdir -p $SORTED $OUT/tmp

# 1. Sort each piece by genome position (skip pieces already sorted)
for bam in $IN/${SM}_*.bam; do
  name=$(basename $bam .bam)
  if [ ! -s $SORTED/$name.sorted.bam ]; then
    $RUN samtools sort -@ $T -m 2G -T $OUT/tmp/$name -o $SORTED/$name.sorted.bam $bam
  fi
  echo "sorted: $name"
done

# 2. Merge the sorted pieces and mark duplicates in one stream
$RUN samtools merge -@ $T -u -o - $SORTED/${SM}_*.sorted.bam \
  | $RUN samtools markdup -@ $T -f $OUT/${SM}.markdup_stats.txt - $OUT/${SM}.grch38.bam

# 3. Index and summarise
$RUN samtools index -@ $T $OUT/${SM}.grch38.bam
$RUN samtools flagstat -@ $T $OUT/${SM}.grch38.bam > $OUT/${SM}.grch38.flagstat.txt
echo "done: $OUT/${SM}.grch38.bam"
