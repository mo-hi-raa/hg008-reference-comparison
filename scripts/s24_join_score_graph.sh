#!/usr/bin/env bash
# JOIN + FILTER + SCORE (C6), graph arm: join s23's 24 per-chromosome VCFs,
# filter Mutect2 once genome-wide, keep PASS calls, score against NIST v0.4 (all.bed).
# Usage: sbatch scripts/s24_join_score_graph.sh deepsomatic   (or: mutect2)
#SBATCH --job-name=join_score_graph
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
ulimit -n $(ulimit -Hn)

CALLER=$1
P=$(pwd -P)
RUN="pixi run --frozen -e pipeline"
T=$SLURM_CPUS_PER_TASK
REF=$P/data/processed/reference/chm13v2/CHM13v2_maskedY_rCRS.fasta
SDF=$P/data/processed/reference/chm13v2/CHM13v2_maskedY_rCRS.sdf
TRUTH=$P/data/raw/HG008-T_somatic_smvar_benchmark_v0.4_CHM13v2_tumorvariants.vcf.gz
TRUTH_BED=$P/data/raw/HG008-T_somatic_smvar_benchmark_v0.4_CHM13v2_all.bed
TRUTH_SAMPLE=dipHG8Tv3.4chr19fixvCHM13v2z200k.hap1.bam
IN=$P/results/calls_graph/$CALLER/by_chrom
OUT=$P/results/calls_graph/$CALLER
CHROMS=(chr{1..22} chrX chrY)

# 1. Join the 24 chromosome files in order (Mutect2: also join stats and filter once)
VCFS=""
if [ "$CALLER" = deepsomatic ]; then
  for c in "${CHROMS[@]}"; do VCFS="$VCFS -I $IN/$c.vcf.gz"; done
  $RUN gatk MergeVcfs $VCFS -O $OUT/calls.vcf.gz
elif [ "$CALLER" = mutect2 ]; then
  STATS=""
  for c in "${CHROMS[@]}"; do
    VCFS="$VCFS -I $IN/$c.raw.vcf.gz"; STATS="$STATS --stats $IN/$c.raw.vcf.gz.stats"
  done
  $RUN gatk MergeVcfs $VCFS -O $OUT/raw.vcf.gz
  $RUN gatk MergeMutectStats $STATS -O $OUT/raw.vcf.gz.stats
  $RUN gatk FilterMutectCalls -R $REF -V $OUT/raw.vcf.gz --stats $OUT/raw.vcf.gz.stats -O $OUT/calls.vcf.gz
else
  echo "unknown caller: $CALLER"; exit 1
fi

# 2. Keep only PASS calls (the fix for aardvark's crash on rejected multi-allele records)
$RUN gatk SelectVariants -V $OUT/calls.vcf.gz --exclude-filtered -O $OUT/calls.pass.vcf.gz
echo "== $CALLER: PASS calls, whole genome: $(zcat $OUT/calls.pass.vcf.gz | grep -vc '^#')"

# 3. Score with vcfeval (main) and show it right away
$RUN rtg vcfeval -b $TRUTH -c $OUT/calls.pass.vcf.gz -t $SDF \
  --bed-regions $TRUTH_BED --sample $TRUTH_SAMPLE,HG008-T \
  --squash-ploidy --threads $T -o $OUT/vcfeval
echo "== $CALLER: vcfeval, whole genome"; cat $OUT/vcfeval/summary.txt

# 4. Score with aardvark (second check); a failure here does not undo the vcfeval result
if $RUN aardvark compare -r $REF -t $TRUTH -q $OUT/calls.pass.vcf.gz \
     -b $TRUTH_BED --truth-sample $TRUTH_SAMPLE --query-sample HG008-T \
     --threads $T -o $OUT/aardvark; then
  echo "== $CALLER: aardvark"; cat $OUT/aardvark/summary.tsv
else
  echo "== $CALLER: aardvark FAILED (vcfeval result above still stands)"
fi
echo "done: graph $CALLER result in $OUT"
