#!/usr/bin/env bash
# Full DSA run: align one read input (line N of metadata/dsa_align_inputs.tsv) to one DSA reference.
# combined_mq60 = V2 (MAPQ reset to 60, original kept in tag mq; Sohn et al. 2026); hap1/hap2 = V3.
# Usage: sbatch --array=1-12%1 scripts/s30_align_dsa.sh combined_mq60   (or hap1, hap2)
#SBATCH --job-name=align_dsa
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=results/logs/%x_%A_%a.out
set -euo pipefail
LEVEL=$1
RUN="pixi run --frozen -e pipeline"
D=data/processed/reference/dsa_v6.3
case $LEVEL in
  combined_mq60) REF=$D/HG008N_v6.3.fasta ;;
  hap1)          REF=$D/HG008N_v6.3.hap1.fasta ;;
  hap2)          REF=$D/HG008N_v6.3.hap2.fasta ;;
  *) echo "unknown reference: $LEVEL"; exit 1 ;;
esac
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/dsa_align_inputs.tsv)
SM=$(echo "$LINE" | cut -f1); R1=$(echo "$LINE" | cut -f2)
R2=${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}; R2=${R2/_R1.fastq.gz/_R2.fastq.gz}
B=$(basename $R1); B=${B%_R1_001.fastq.gz}; B=${B%_R1.fastq.gz}; ID=${B%_part??}
OUT=data/processed/alignments/dsa_$LEVEL/unsorted; mkdir -p $OUT
F=$OUT/${SM}_${B}.bam
if [ -s $F ]; then echo "already done: $F"; exit 0; fi
echo "task $SLURM_ARRAY_TASK_ID: $SM $B (read group $ID) -> $LEVEL"
RG="@RG\tID:${ID}\tSM:${SM}\tPL:ILLUMINA\tPU:${ID}"
if [ $LEVEL = combined_mq60 ]; then
  $RUN bwa mem -t 30 -R "$RG" $REF $R1 $R2 \
    | awk 'BEGIN{OFS="\t"} /^@/ {print; next} { if (int($2/4)%2==0) { o=$5; $5=60; $0=$0 "\tmq:i:" o } print }' \
    | $RUN samtools fixmate -m -O bam -@ 2 - $F.part
else
  $RUN bwa mem -t 30 -R "$RG" $REF $R1 $R2 | $RUN samtools fixmate -m -O bam -@ 2 - $F.part
fi
mv $F.part $F
echo "done: $F"
