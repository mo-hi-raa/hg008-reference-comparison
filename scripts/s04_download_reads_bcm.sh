#!/usr/bin/env bash
#SBATCH --job-name=dl_bcm_fastq
#SBATCH --partition=ws-ia
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=23:00:00
#SBATCH --array=1-10%5
#SBATCH --output=results/logs/%x_%A_%a.out
set -euo pipefail

BASE=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313
MD5=metadata/BCM_Illumina-WGS_20240313_checksum.md5
OUT=data/raw/reads/BCM_Illumina-WGS_20240313

line=$(grep "fastq.gz" "$MD5" | sed -n "${SLURM_ARRAY_TASK_ID}p")
path=$(echo "$line" | awk '{print $2}')

mkdir -p "$OUT/$(dirname "$path")"
curl -fL --retry 5 --retry-delay 60 -C - -o "$OUT/$path.part" "$BASE/$path"
mv "$OUT/$path.part" "$OUT/$path"

cd "$OUT"
echo "$line" | md5sum -c -
echo "done: $path"
