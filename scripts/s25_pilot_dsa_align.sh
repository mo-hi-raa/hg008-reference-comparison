#!/usr/bin/env bash
# DSA chr22 pilot: align one sample's chr22 reads to one DSA reference
# (combined / hap1 / hap2), fixmate, sort, mark duplicates, index,
# then count mapped and MAPQ-0 reads on chr22.
# Usage: sbatch scripts/s25_pilot_dsa_align.sh HG008-T combined
#SBATCH --job-name=pilot_dsa_align
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
SM=$1; LEVEL=$2
RUN="pixi run --frozen -e pipeline"
P=data/processed/dsa_pilot
D=data/processed/reference/dsa_v6.3
case $LEVEL in
  combined) REF=$D/HG008N_v6.3.fasta;      REGIONS="chr22_hap1 chr22_hap2" ;;
  hap1)     REF=$D/HG008N_v6.3.hap1.fasta; REGIONS="chr22_hap1" ;;
  hap2)     REF=$D/HG008N_v6.3.hap2.fasta; REGIONS="chr22_hap2" ;;
esac
OUT=$P/bams/${SM}.dsa_${LEVEL}.bam
RG="@RG\tID:${SM}_pilot\tSM:${SM}\tPL:ILLUMINA"
$RUN bwa mem -t 14 -R "$RG" $REF $P/reads/${SM}_R1.fastq.gz $P/reads/${SM}_R2.fastq.gz \
  | $RUN samtools fixmate -m -u - - \
  | $RUN samtools sort -u -@ 2 -m 2G -T $P/tmp/${SM}_${LEVEL} - \
  | $RUN samtools markdup -@ 2 - $OUT
$RUN samtools index $OUT
ALL=$($RUN samtools view -c -F 0x904 $OUT $REGIONS)
Q1=$($RUN samtools view -c -F 0x904 -q 1 $OUT $REGIONS)
echo "check a: $SM $LEVEL mapped on chr22 = $ALL, with MAPQ>=1 = $Q1, MAPQ 0 = $((ALL - Q1))"
echo "done: $OUT"
