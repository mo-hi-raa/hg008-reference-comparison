#!/usr/bin/env bash
# Finish the graph-arm alignment for ONE sample: rename chromosomes
# (CHM13#0#chr1 -> chr1, to match the CHM13 reference), sort each piece,
# merge the pieces, mark duplicates, index.
# Usage: sbatch scripts/s22_finish_bam_graph.sh HG008-T
#SBATCH --job-name=finish_bam_graph
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail

SM=$1
RUN="pixi run --frozen -e pipeline"
IN=data/processed/alignments/graph/unsorted
SORTED=data/processed/alignments/graph/sorted_pieces
OUT=data/processed/alignments/graph
T=$SLURM_CPUS_PER_TASK
mkdir -p $SORTED $OUT/tmp

# 1. Rename chromosomes and sort each piece (skip pieces already sorted)
for bam in $IN/${SM}_*.bam; do
  name=$(basename $bam .bam)
  if [ ! -s $SORTED/$name.sorted.bam ]; then
    $RUN samtools view -H $bam | sed 's/SN:CHM13#0#/SN:/' > $OUT/tmp/$name.header.sam
    $RUN samtools reheader $OUT/tmp/$name.header.sam $bam \
      | $RUN samtools sort -@ $T -m 2G -T $OUT/tmp/$name -o $SORTED/$name.sorted.bam -
  fi
  echo "sorted: $name"
done

# 2. Merge the sorted pieces and mark duplicates in one stream
$RUN samtools merge -@ $T -u -o - $SORTED/${SM}_*.sorted.bam \
  | $RUN samtools markdup -@ $T -f $OUT/${SM}.markdup_stats.txt - $OUT/${SM}.graph.bam

# 3. Index, summarise, and show the first chromosome names
$RUN samtools index -@ $T $OUT/${SM}.graph.bam
$RUN samtools flagstat -@ $T $OUT/${SM}.graph.bam > $OUT/${SM}.graph.flagstat.txt
$RUN samtools view -H $OUT/${SM}.graph.bam | grep "^@SQ" | head -3 || true
echo "done: $OUT/${SM}.graph.bam"
