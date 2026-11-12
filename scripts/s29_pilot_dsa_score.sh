#!/usr/bin/env bash
# DSA pilot scoring (Way B, haplotype-agnostic): twin truth, V3 hap1+hap2 merge,
# rtg vcfeval per version and caller, counted with each mutation once.
# Usage: sbatch scripts/s29_pilot_dsa_score.sh dsa_pilot_chr10 chr10
#SBATCH --job-name=pilot_dsa_score
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
PILOT=$1; CHR=$2
RUN="pixi run --frozen -e pipeline"
P=$(pwd -P)
D=$P/data/processed/reference/dsa_v6.3
E=$P/results/$PILOT/eval
C=$P/results/$PILOT/calls
SDF=$D/HG008N_v6.3.sdf
[ -d $SDF ] || $RUN rtg format -o $SDF $D/HG008N_v6.3.fasta
python3 -I $P/scripts/pilot_lift.py $E $CHR
$RUN gatk SortVcf -I $E/truth_twins.vcf -O $E/truth_twins.vcf.gz -SD $D/HG008N_v6.3.dict --CREATE_INDEX true
for V in V1 V2 V3; do
  for CALLER in deepsomatic mutect2; do
    if [ $V = V3 ]; then
      mkdir -p $C/V3/$CALLER
      { zcat $C/V3hap1/$CALLER/pass.vcf.gz; zcat $C/V3hap2/$CALLER/pass.vcf.gz | grep -v '^#' || true; } > $C/V3/$CALLER/pass_merged.vcf
      $RUN gatk SortVcf -I $C/V3/$CALLER/pass_merged.vcf -O $C/V3/$CALLER/pass.vcf.gz -SD $D/HG008N_v6.3.dict --CREATE_INDEX true
    fi
    O=$E/vcfeval_${V}_$CALLER; rm -rf $O
    if $RUN rtg vcfeval -b $E/truth_twins.vcf.gz -c $C/$V/$CALLER/pass.vcf.gz -t $SDF \
         --bed-regions $E/eval.bed --squash-ploidy --sample TRUTH,HG008-T -o $O > $O.log 2>&1; then
      python3 -I $P/scripts/pilot_count.py $O "$CHR $V $CALLER"
    else
      echo "vcfeval FAILED: $CHR $V $CALLER"; tail -n 5 $O.log
    fi
  done
done
echo "done: $CHR"
