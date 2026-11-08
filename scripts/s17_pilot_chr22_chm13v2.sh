#!/usr/bin/env bash
# PILOT, T2T-CHM13v2.0: call and score somatic variants on chr22 using the real
# finished CHM13 BAMs. Same settings as the GRCh38 pilot (s14), plus PASS-only scoring.
# Usage: sbatch scripts/s17_pilot_chr22_chm13v2.sh deepsomatic   (or: mutect2)
#SBATCH --job-name=pilot_chr22_chm13v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
ulimit -n $(ulimit -Hn)

CALLER=$1
P=$(pwd -P)
RUN="pixi run --frozen -e pipeline"
T=$SLURM_CPUS_PER_TASK
REF=$P/data/processed/reference/chm13v2/CHM13v2_maskedY_rCRS.fasta
SDF=$P/data/processed/reference/chm13v2/CHM13v2_maskedY_rCRS.sdf
T_BAM=$P/data/processed/alignments/chm13v2/HG008-T.chm13v2.bam
N_BAM=$P/data/processed/alignments/chm13v2/HG008-N-D.chm13v2.bam
TRUTH=$P/data/raw/HG008-T_somatic_smvar_benchmark_v0.4_CHM13v2_tumorvariants.vcf.gz
TRUTH_BED=$P/data/raw/HG008-T_somatic_smvar_benchmark_v0.4_CHM13v2_all.bed
TRUTH_SAMPLE=dipHG8Tv3.4chr19fixvCHM13v2z200k.hap1.bam
OUT=$P/results/pilot_chm13v2_chr22/$CALLER
mkdir -p $OUT
awk '$1=="chr22"' $TRUTH_BED > $OUT/benchmark_chr22.bed

# 1. Call on chr22, and time it
start=$SECONDS
if [ "$CALLER" = deepsomatic ]; then
  mkdir -p $OUT/tmp
  singularity exec -B /storage containers/deepsomatic_1.8.0.sif run_deepsomatic \
    --model_type=WGS --ref=$REF --reads_tumor=$T_BAM --reads_normal=$N_BAM \
    --sample_name_tumor=HG008-T --sample_name_normal=HG008-N-D \
    --regions=chr22 --num_shards=$T --intermediate_results_dir=$OUT/tmp \
    --output_vcf=$OUT/calls.vcf.gz
elif [ "$CALLER" = mutect2 ]; then
  $RUN gatk Mutect2 -R $REF -I $T_BAM -I $N_BAM -normal HG008-N-D -L chr22 \
    --native-pair-hmm-threads $T -O $OUT/raw.vcf.gz
  $RUN gatk FilterMutectCalls -R $REF -V $OUT/raw.vcf.gz -O $OUT/calls.vcf.gz
else
  echo "unknown caller: $CALLER"; exit 1
fi
echo "TIME: $CALLER on chr22 took $(( (SECONDS - start) / 60 )) minutes"
$RUN gatk SelectVariants -V $OUT/calls.vcf.gz --exclude-filtered -O $OUT/calls.pass.vcf.gz
echo "== $CALLER: PASS calls on chr22: $(zcat $OUT/calls.pass.vcf.gz | grep -vc '^#')"

# 2. Score with vcfeval (main), then aardvark (second check)
$RUN rtg vcfeval -b $TRUTH -c $OUT/calls.pass.vcf.gz -t $SDF \
  --bed-regions $OUT/benchmark_chr22.bed --sample $TRUTH_SAMPLE,HG008-T \
  --squash-ploidy --threads $T -o $OUT/vcfeval
echo "== $CALLER: vcfeval, chr22"; cat $OUT/vcfeval/summary.txt
if $RUN aardvark compare -r $REF -t $TRUTH -q $OUT/calls.pass.vcf.gz \
     -b $OUT/benchmark_chr22.bed --truth-sample $TRUTH_SAMPLE --query-sample HG008-T \
     --threads $T -o $OUT/aardvark; then
  echo "== $CALLER: aardvark"; cat $OUT/aardvark/summary.tsv
else
  echo "== $CALLER: aardvark FAILED (vcfeval result above still stands)"
fi
echo "done: T2T pilot $CALLER in $OUT"
