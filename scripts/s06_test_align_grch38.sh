#!/usr/bin/env bash
# PIPELINE TEST (decision C3, option 2). These numbers are NOT results.
# First 2 million read pairs of one tumour FASTQ pair and of the normal pair,
# read straight from the BCM server, aligned to the full GRCh38; chr22 kept
# for the calling and scoring tests.
# Run from the project folder with: sbatch scripts/s06_test_align_grch38.sh
#SBATCH --job-name=test_align_grch38
#SBATCH --partition=ws-ia
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail

BASE=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313
REF=data/processed/reference/grch38/GRCh38_GIABv3
OUT=results/test_grch38
LINES=8000000
RUN="pixi run --frozen -e pipeline"
mkdir -p "$OUT"

# 1. Test reads: first 8,000,000 lines = 2,000,000 reads (each read is 4 lines)
(curl -sL "$BASE/HG008-T_fastqs/HJVHNDSX7-1-IDUDI0031_S12_L001_R1_001.fastq.gz" | zcat | head -n $LINES || true) | gzip > $OUT/T_R1.fastq.gz
(curl -sL "$BASE/HG008-T_fastqs/HJVHNDSX7-1-IDUDI0031_S12_L001_R2_001.fastq.gz" | zcat | head -n $LINES || true) | gzip > $OUT/T_R2.fastq.gz
(curl -sL "$BASE/HG008-N-D_fastqs/HV5TMDSX7-1-IDUDI0034_S1_L001_R1_001.fastq.gz" | zcat | head -n $LINES || true) | gzip > $OUT/N_R1.fastq.gz
(curl -sL "$BASE/HG008-N-D_fastqs/HV5TMDSX7-1-IDUDI0034_S1_L001_R2_001.fastq.gz" | zcat | head -n $LINES || true) | gzip > $OUT/N_R2.fastq.gz

for f in T_R1 T_R2 N_R1 N_R2; do
  n=$(zcat $OUT/$f.fastq.gz | wc -l)
  echo "$f: $n lines"
  [ "$n" -eq "$LINES" ]
done

# 2. Align to the full GRCh38, sort, index
$RUN bwa mem -t 16 -R '@RG\tID:HG008-T_test\tSM:HG008-T\tPL:ILLUMINA' $REF $OUT/T_R1.fastq.gz $OUT/T_R2.fastq.gz | $RUN samtools sort -@ 4 -o $OUT/HG008-T_test.bam -
$RUN bwa mem -t 16 -R '@RG\tID:HG008-N-D_test\tSM:HG008-N-D\tPL:ILLUMINA' $REF $OUT/N_R1.fastq.gz $OUT/N_R2.fastq.gz | $RUN samtools sort -@ 4 -o $OUT/HG008-N-D_test.bam -

for s in HG008-T HG008-N-D; do
  $RUN samtools index $OUT/${s}_test.bam
  $RUN samtools flagstat $OUT/${s}_test.bam > $OUT/${s}_test.flagstat.txt
  # 3. Keep chr22 only, for the calling and scoring tests
  $RUN samtools view -b -o $OUT/${s}_test_chr22.bam $OUT/${s}_test.bam chr22
  $RUN samtools index $OUT/${s}_test_chr22.bam
done
echo "done: test BAMs in $OUT"
