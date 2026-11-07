#!/usr/bin/env bash
# PILOT (C5b): call and score somatic variants on chr22 using the REAL
# finished GRCh38 BAMs. Main setting: Mutect2 without gnomAD/panel of normals.
# Usage: sbatch scripts/s14_pilot_chr22_grch38.sh deepsomatic   (or: mutect2)
#SBATCH --job-name=pilot_chr22_grch38
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
REF=$P/data/processed/reference/grch38/GRCh38_GIABv3.fasta
SDF=$P/data/processed/reference/grch38/GRCh38_GIABv3.sdf
T_BAM=$P/data/processed/alignments/grch38/HG008-T.grch38.bam
N_BAM=$P/data/processed/alignments/grch38/HG008-N-D.grch38.bam
TRUTH=$P/data/raw/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz
TRUTH_BED=$P/data/raw/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed
TRUTH_SAMPLE=dipHG8Tv3.4chr19fixvGRCh38z200k.hap1.bam
OUT=$P/results/pilot_grch38_chr22/$CALLER
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

# 2. Score with vcfeval and aardvark
$RUN rtg vcfeval -b $TRUTH -c $OUT/calls.vcf.gz -t $SDF \
  --bed-regions $OUT/benchmark_chr22.bed --sample $TRUTH_SAMPLE,HG008-T \
  --squash-ploidy --threads $T -o $OUT/vcfeval
$RUN aardvark compare -r $REF -t $TRUTH -q $OUT/calls.vcf.gz \
  -b $OUT/benchmark_chr22.bed --truth-sample $TRUTH_SAMPLE --query-sample HG008-T \
  --threads $T -o $OUT/aardvark

# 3. Show the results
echo "== $CALLER: PASS calls on chr22: $(zcat $OUT/calls.vcf.gz | grep -v '^#' | awk '$7=="PASS"' | wc -l)"
echo "== $CALLER: vcfeval"; cat $OUT/vcfeval/summary.txt
echo "== $CALLER: aardvark"; cat $OUT/aardvark/summary.tsv
echo "done: pilot $CALLER in $OUT"
