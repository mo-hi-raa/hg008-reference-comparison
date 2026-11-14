#!/usr/bin/env bash
# Full DSA run: join per-chromosome calls, Mutect2 filtering, PASS only, score against the twin truth (Way B).
# V2 = combined_mq60 calls; V3 = hap1 + hap2 calls joined.
# Usage: sbatch scripts/s34_dsa_score.sh V2   (or V3)
#SBATCH --job-name=dsa_score
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
V=$1
RUN="pixi run --frozen -e pipeline"
P=$(pwd -P)
D=$P/data/processed/reference/dsa_v6.3
F=$P/results/dsa_eval
case $V in V2) LEVELS="combined_mq60" ;; V3) LEVELS="hap1 hap2" ;; *) echo "unknown version: $V"; exit 1 ;; esac
pass_calls() {
  L=$1; CALLER=$2; C=$P/results/calls_dsa_$L/$CALLER
  case $L in combined_mq60) REF=$D/HG008N_v6.3.fasta ;; *) REF=$D/HG008N_v6.3.$L.fasta ;; esac
  VCFS=""; STATS=""
  for CHR in chr{1..22} chrX; do
    if [ $CALLER = deepsomatic ]; then VCFS="$VCFS -I $C/by_chrom/$CHR.vcf.gz"
    else VCFS="$VCFS -I $C/by_chrom/$CHR.raw.vcf.gz"; STATS="$STATS -stats $C/by_chrom/$CHR.raw.vcf.gz.stats"; fi
  done
  if [ $CALLER = deepsomatic ]; then
    $RUN gatk MergeVcfs $VCFS -O $C/calls.vcf.gz
  else
    $RUN gatk MergeVcfs $VCFS -O $C/raw.vcf.gz
    $RUN gatk MergeMutectStats $STATS -O $C/raw.vcf.gz.stats
    $RUN gatk FilterMutectCalls -R $REF -V $C/raw.vcf.gz --stats $C/raw.vcf.gz.stats -O $C/calls.vcf.gz
  fi
  $RUN gatk SelectVariants -R $REF -V $C/calls.vcf.gz --exclude-filtered -O $C/pass.vcf.gz
  echo "PASS calls: $L $CALLER $(zcat $C/pass.vcf.gz | grep -vc '^#')"
}
for CALLER in deepsomatic mutect2; do
  for L in $LEVELS; do pass_calls $L $CALLER; done
  if [ $V = V2 ]; then
    CALLS=$P/results/calls_dsa_combined_mq60/$CALLER/pass.vcf.gz
  else
    M=$F/calls_V3_$CALLER.vcf; A=$P/results/calls_dsa_hap1/$CALLER/pass.vcf.gz; B=$P/results/calls_dsa_hap2/$CALLER/pass.vcf.gz
    { zcat $A | grep '^##' || true; zcat $B | grep '^##contig' || true; zcat $A | grep '^#CHROM' || true
      zcat $A | grep -v '^#' || true; zcat $B | grep -v '^#' || true; } > $M
    rm -f $M.gz $M.gz.tbi; $RUN bgzip -f $M; $RUN tabix -p vcf $M.gz; CALLS=$M.gz
  fi
  O=$F/vcfeval_${V}_$CALLER; rm -rf $O
  if $RUN rtg vcfeval -b $F/truth_twins.vcf.gz -c $CALLS -t $D/HG008N_v6.3.sdf --bed-regions $F/eval.bed \
       --squash-ploidy --sample TRUTH,HG008-T --threads 8 -o $O > $O.log 2>&1; then
    python3 -I $P/scripts/pilot_count.py $O "genome $V $CALLER"
  else
    echo "vcfeval FAILED: $V $CALLER"; tail -n 5 $O.log
  fi
done
echo "done: $V"
