#!/usr/bin/env bash
# DSA chr10 pilot 2b prep: Mutect2 filtering, PASS-only calls for both callers,
# chr10 subset of the native v6.3 truth, and chr10 hap2<->hap1 alignments (for twins).
#SBATCH --job-name=pilot_dsa_2b_prep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
RUN="pixi run --frozen -e pipeline"
P=$(pwd -P)
D=$P/data/processed/reference/dsa_v6.3
C=$P/results/dsa_pilot_chr10/calls
E=$P/results/dsa_pilot_chr10/eval
mkdir -p $E
# 1. Mutect2 filtering, then PASS-only calls for both callers
for V in V1 V2 V3hap1 V3hap2; do
  case $V in
    V3hap1) REF=$D/HG008N_v6.3.hap1.fasta ;;
    V3hap2) REF=$D/HG008N_v6.3.hap2.fasta ;;
    *)      REF=$D/HG008N_v6.3.fasta ;;
  esac
  M=$C/$V/mutect2
  $RUN gatk FilterMutectCalls -R $REF -V $M/raw.vcf.gz --stats $M/raw.vcf.gz.stats -O $M/calls.vcf.gz
  for CALLER in deepsomatic mutect2; do
    $RUN gatk SelectVariants -R $REF -V $C/$V/$CALLER/calls.vcf.gz --exclude-filtered -O $C/$V/$CALLER/pass.vcf.gz
    echo "PASS calls: $V $CALLER $(zcat $C/$V/$CALLER/pass.vcf.gz | grep -vc '^#')"
  done
done
# 2. chr10 part of the native v6.3 truth and its regions
T=$P/data/raw/benchmark/HG008Nv6.3/HG008-T_somatic_smvar_benchmark_v0.4_HG008Nv6.3
zcat ${T}_tumorvariants.vcf.gz | awk '/^#/ || $1=="chr10_hap1" || $1=="chr10_hap2"' > $E/truth.chr10.vcf
awk '$1=="chr10_hap1" || $1=="chr10_hap2"' ${T}.bed > $E/truth.chr10.bed
grep -v '^#' $E/truth.chr10.vcf | awk '{t=(length($4)==1 && length($5)==1) ? "SNV" : "indel"; n[$1" "t]++} END {for (k in n) print "truth:", k, n[k]}'
# 3. line chr10_hap2 up against chr10_hap1, and the reverse
$RUN samtools faidx $D/HG008N_v6.3.fasta chr10_hap1 > $E/chr10_hap1.fa
$RUN samtools faidx $D/HG008N_v6.3.fasta chr10_hap2 > $E/chr10_hap2.fa
$RUN minimap2 -cx asm5 --cs -t 8 $E/chr10_hap1.fa $E/chr10_hap2.fa > $E/hap2_to_hap1.paf
$RUN minimap2 -cx asm5 --cs -t 8 $E/chr10_hap2.fa $E/chr10_hap1.fa > $E/hap1_to_hap2.paf
echo "largest hap2->hap1 alignment blocks (query, target, strand, length, mapq):"
awk '{print $1, $6, $5, $4-$3, $12}' $E/hap2_to_hap1.paf | sort -k4,4nr | head -5
echo done
