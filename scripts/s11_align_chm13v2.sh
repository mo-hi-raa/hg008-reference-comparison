#!/usr/bin/env bash
# Align the BCM 2024 Illumina reads to T2T-CHM13v2.0 (maskedY_rCRS) (real run, step C4).
# Array task 1 = normal (1 FASTQ pair); tasks 2-5 = the 4 tumour FASTQ pairs.
# Each pair -> one BAM with mate tags (fixmate -m), ready for sorting,
# merging and duplicate marking in the next step.
# Run from the project folder with: sbatch scripts/s08_align_grch38.sh
#SBATCH --job-name=align_chm13v2
#SBATCH --partition=ws-ia
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=23:59:00
#SBATCH --array=1-5%2
#SBATCH --output=results/logs/%x_%A_%a.out
set -euo pipefail

RUN="pixi run --frozen -e pipeline"
REF=data/processed/reference/chm13v2/CHM13v2_maskedY_rCRS
READS=data/raw/reads/BCM_Illumina-WGS_20240313
MD5=metadata/BCM_Illumina-WGS_20240313_checksum.md5
OUT=data/processed/alignments/chm13v2/unsorted
mkdir -p $OUT

r1=$(grep "_R1_001.fastq.gz" $MD5 | sed -n "${SLURM_ARRAY_TASK_ID}p" | awk '{print $2}')
r2=${r1/_R1_001/_R2_001}
case $r1 in
  HG008-N-D_fastqs/*) SM=HG008-N-D ;;
  HG008-T_fastqs/*)   SM=HG008-T ;;
esac
ID=$(basename $r1 _R1_001.fastq.gz)
echo "task $SLURM_ARRAY_TASK_ID: sample $SM, read group $ID"

$RUN bwa mem -t $((SLURM_CPUS_PER_TASK - 2)) -R "@RG\tID:$ID\tSM:$SM\tPL:ILLUMINA\tPU:$ID" $REF $READS/$r1 $READS/$r2 \
  | $RUN samtools fixmate -m -@ 2 - $OUT/${SM}_${ID}.bam
echo "done: $OUT/${SM}_${ID}.bam"
