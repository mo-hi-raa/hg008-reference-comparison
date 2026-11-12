#!/usr/bin/env bash
# DSA chr10 pilot prep: .fai + .dict for the 3 DSA references, chr10 BED files,
# and V2 BAMs (combined BAMs with MAPQ reset to 60, original kept in tag mq; Sohn et al. 2026).
#SBATCH --job-name=pilot_dsa_prep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
RUN="pixi run --frozen -e pipeline"
P=data/processed/dsa_pilot_chr10
D=data/processed/reference/dsa_v6.3
for L in "" .hap1 .hap2; do
  F=$D/HG008N_v6.3$L.fasta
  [ -s $F.fai ] || $RUN samtools faidx $F
  [ -s ${F%.fasta}.dict ] || $RUN gatk CreateSequenceDictionary -R $F
done
mkdir -p $P/beds
awk -v OFS='\t' '$1=="chr10_hap1" || $1=="chr10_hap2" {print $1, 0, $2}' $D/HG008N_v6.3.fasta.fai > $P/beds/combined.chr10.bed
awk -v OFS='\t' '$1=="chr10_hap1" {print $1, 0, $2}' $D/HG008N_v6.3.fasta.fai > $P/beds/hap1.chr10.bed
awk -v OFS='\t' '$1=="chr10_hap2" {print $1, 0, $2}' $D/HG008N_v6.3.fasta.fai > $P/beds/hap2.chr10.bed
for SM in HG008-T HG008-N-D; do
  IN=$P/bams/$SM.dsa_combined.bam
  OUT=$P/bams/$SM.dsa_combined_mq60.bam
  $RUN samtools view -h $IN \
    | awk 'BEGIN{OFS="\t"} /^@/ {print; next} { if (int($2/4)%2==0) { o=$5; $5=60; for (i=12; i<=NF; i++) if ($i ~ /^MQ:i:/) $i="MQ:i:60"; $0=$0 "\tmq:i:" o } print }' \
    | $RUN samtools view -b -@ 4 -o $OUT -
  $RUN samtools index $OUT
done
cat $P/beds/*.bed
echo done
