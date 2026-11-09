#!/usr/bin/env bash
# Align the BCM 2024 Illumina reads to the HPRC v2.1 CHM13-based graph (giraffe).
# Same structure as s08: task 1 = normal; tasks 2-5 = the 4 tumour FASTQ pairs.
# giraffe surjects each read onto CHM13 (--ref-paths metadata/chm13_ref_paths.txt -o BAM); fixmate -m prepares for markdup.
# Run from the project folder with: sbatch scripts/s21_align_graph.sh
#SBATCH --job-name=align_graph
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=23:59:00
#SBATCH --array=1-5%2
#SBATCH --output=results/logs/%x_%A_%a.out
set -euo pipefail

RUN="pixi run --frozen -e pipeline"
G=data/raw/graph/hprc-v2.1-mc-chm13/hprc-v2.1-mc-chm13.d46
READS=data/raw/reads/BCM_Illumina-WGS_20240313
MD5=metadata/BCM_Illumina-WGS_20240313_checksum.md5
OUT=data/processed/alignments/graph/unsorted
mkdir -p $OUT

r1=$(grep "_R1_001.fastq.gz" $MD5 | sed -n "${SLURM_ARRAY_TASK_ID}p" | awk '{print $2}')
r2=${r1/_R1_001/_R2_001}
case $r1 in
  HG008-N-D_fastqs/*) SM=HG008-N-D ;;
  HG008-T_fastqs/*)   SM=HG008-T ;;
esac
ID=$(basename $r1 _R1_001.fastq.gz)
echo "task $SLURM_ARRAY_TASK_ID: sample $SM, read group $ID"

$RUN vg giraffe -t $((SLURM_CPUS_PER_TASK - 2)) -Z $G.gbz -d $G.dist \
  -m $G.shortread.withzip.min -z $G.shortread.zipcodes \
  -f $READS/$r1 -f $READS/$r2 -N $SM -R $ID --ref-paths metadata/chm13_ref_paths.txt -o BAM \
  | $RUN samtools fixmate -m -@ 2 - $OUT/${SM}_${ID}.bam
echo "done: $OUT/${SM}_${ID}.bam"
