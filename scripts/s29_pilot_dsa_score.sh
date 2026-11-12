#!/usr/bin/env bash
# DSA pilot scoring (Way B): twin truth, V3 hap1+hap2 merge, vcfeval, each mutation counted once.
# VCFs compressed with bgzip and indexed with tabix (GATK SortVcf refuses the NIST truth header).
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
rm -f $E/truth_twins.vcf.gz $E/truth_twins.vcf.gz.tbi
$RUN bgzip -f $E/truth_twins.vcf
$RUN tabix -p vcf $E/truth_twins.vcf.gz
for V in V1 V2 V3; do
  for CALLER in deepsomatic mutect2; do
    if [ $V = V3 ]; then
      M=$C/V3/$CALLER; mkdir -p $M; rm -f $M/pass_merged.vcf.gz $M/pass_merged.vcf.gz.tbi
      A=$C/V3hap1/$CALLER/pass.vcf.gz; B=$C/V3hap2/$CALLER/pass.vcf.gz
      { zcat $A | grep '^##' || true; zcat $B | grep '^##contig' || true; zcat $A | grep '^#CHROM' || true
        zcat $A | grep -v '^#' || true; zcat $B | grep -v '^#' || true; } > $M/pass_merged.vcf
      $RUN bgzip -f $M/pass_merged.vcf
      $RUN tabix -p vcf $M/pass_merged.vcf.gz
      CALLS=$M/pass_merged.vcf.gz
    else
      CALLS=$C/$V/$CALLER/pass.vcf.gz
    fi
    O=$E/vcfeval_${V}_$CALLER; rm -rf $O
    if $RUN rtg vcfeval -b $E/truth_twins.vcf.gz -c $CALLS -t $SDF --bed-regions $E/eval.bed --squash-ploidy --sample TRUTH,HG008-T -o $O > $O.log 2>&1; then
      python3 -I $P/scripts/pilot_count.py $O "$CHR $V $CALLER"
    else
      echo "vcfeval FAILED: $CHR $V $CALLER"; tail -n 5 $O.log
    fi
  done
done
echo "done: $CHR"
