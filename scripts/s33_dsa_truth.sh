#!/usr/bin/env bash
# Full DSA run: genome-wide haplotype-agnostic ("twin") truth for scoring (Way B, as in the pilot).
# Per chromosome: split truth + regions, align chrN_hap2 <-> chrN_hap1 (minimap2), make twins (pilot_lift.py).
#SBATCH --job-name=dsa_truth
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
RUN="pixi run --frozen -e pipeline"
P=$(pwd -P)
D=$P/data/processed/reference/dsa_v6.3
T=$P/data/raw/benchmark/HG008Nv6.3/HG008-T_somatic_smvar_benchmark_v0.4_HG008Nv6.3
F=$P/results/dsa_eval; E=$F/truth; mkdir -p $E
for CHR in chr{1..22} chrX; do
  C=$E/$CHR; mkdir -p $C
  zcat ${T}_tumorvariants.vcf.gz | awk -v a=${CHR}_hap1 -v b=${CHR}_hap2 '/^#/ || $1==a || $1==b' > $C/truth.$CHR.vcf
  awk -v a=${CHR}_hap1 -v b=${CHR}_hap2 '$1==a || $1==b' ${T}.bed > $C/truth.$CHR.bed
  if [ ! -s $C/hap1_to_hap2.paf ]; then
    $RUN samtools faidx $D/HG008N_v6.3.fasta ${CHR}_hap1 > $C/${CHR}_hap1.fa
    $RUN samtools faidx $D/HG008N_v6.3.fasta ${CHR}_hap2 > $C/${CHR}_hap2.fa
    $RUN minimap2 -cx asm5 --cs -t 8 $C/${CHR}_hap1.fa $C/${CHR}_hap2.fa > $C/hap2_to_hap1.paf 2> $C/minimap2.log
    $RUN minimap2 -cx asm5 --cs -t 8 $C/${CHR}_hap2.fa $C/${CHR}_hap1.fa > $C/hap1_to_hap2.paf 2>> $C/minimap2.log
  fi
  python3 -I $P/scripts/pilot_lift.py $C $CHR
done
{ grep '^#' $E/chr1/truth_twins.vcf
  for CHR in chr{1..22} chrX; do grep -v '^#' $E/$CHR/truth_twins.vcf || true; done; } > $F/truth_twins.vcf
for CHR in chr{1..22} chrX; do cat $E/$CHR/eval.bed; done > $F/eval.bed
rm -f $F/truth_twins.vcf.gz $F/truth_twins.vcf.gz.tbi
$RUN bgzip -f $F/truth_twins.vcf
$RUN tabix -p vcf $F/truth_twins.vcf.gz
echo "twin truth records: $(zcat $F/truth_twins.vcf.gz | grep -vc '^#')"
echo "done: twin truth in $F"
