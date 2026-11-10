#!/usr/bin/env bash
# WHOLE-GENOME CALLING (C5c), graph arm (reads on CHM13): one array task per chromosome.
# Joining, Mutect2 filtering and scoring come next in s16.
# Usage: sbatch --array=1-24%N scripts/s23_call_graph.sh deepsomatic   (or: mutect2)
#SBATCH --job-name=call_graph
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/%x_%A_%a.out
set -euo pipefail
ulimit -n $(ulimit -Hn)

CALLER=$1
CHROMS=(chr{1..22} chrX chrY)
CHR=${CHROMS[$((SLURM_ARRAY_TASK_ID - 1))]}
P=$(pwd -P)
RUN="pixi run --frozen -e pipeline"
T=$SLURM_CPUS_PER_TASK
REF=$P/data/processed/reference/chm13v2/CHM13v2_maskedY_rCRS.fasta
T_BAM=$P/data/processed/alignments/graph/HG008-T.graph.bam
N_BAM=$P/data/processed/alignments/graph/HG008-N-D.graph.bam
OUT=$P/results/calls_graph/$CALLER/by_chrom
mkdir -p $OUT

start=$SECONDS
if [ "$CALLER" = deepsomatic ]; then
  mkdir -p $OUT/tmp_$CHR
  singularity exec -B /storage containers/deepsomatic_1.8.0.sif run_deepsomatic \
    --model_type=WGS --ref=$REF --reads_tumor=$T_BAM --reads_normal=$N_BAM \
    --sample_name_tumor=HG008-T --sample_name_normal=HG008-N-D \
    --regions=$CHR --num_shards=$T --intermediate_results_dir=$OUT/tmp_$CHR \
    --output_vcf=$OUT/$CHR.vcf.gz
elif [ "$CALLER" = mutect2 ]; then
  $RUN gatk Mutect2 -R $REF -I $T_BAM -I $N_BAM -normal HG008-N-D -L $CHR \
    --native-pair-hmm-threads $T -O $OUT/$CHR.raw.vcf.gz
else
  echo "unknown caller: $CALLER"; exit 1
fi
echo "TIME: $CALLER on $CHR took $(( (SECONDS - start) / 60 )) minutes"
echo "done: $CALLER $CHR"
