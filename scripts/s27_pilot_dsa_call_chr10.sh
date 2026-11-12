#!/usr/bin/env bash
# DSA chr10 pilot calling: same DeepSomatic / Mutect2 commands as s15, on chr10 of the DSA.
# Usage: sbatch --cpus-per-task=N scripts/s27_pilot_dsa_call.sh deepsomatic|mutect2 V1|V2|V3hap1|V3hap2
#SBATCH --job-name=pilot_dsa_call
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail
ulimit -n $(ulimit -Hn)
CALLER=$1; V=$2
P=$(pwd -P)
RUN="pixi run --frozen -e pipeline"
T=$SLURM_CPUS_PER_TASK
D=$P/data/processed/reference/dsa_v6.3
B=$P/data/processed/dsa_pilot_chr10
case $V in
  V1)     S=dsa_combined;      REF=$D/HG008N_v6.3.fasta;      BED=$B/beds/combined.chr10.bed ;;
  V2)     S=dsa_combined_mq60; REF=$D/HG008N_v6.3.fasta;      BED=$B/beds/combined.chr10.bed ;;
  V3hap1) S=dsa_hap1;          REF=$D/HG008N_v6.3.hap1.fasta; BED=$B/beds/hap1.chr10.bed ;;
  V3hap2) S=dsa_hap2;          REF=$D/HG008N_v6.3.hap2.fasta; BED=$B/beds/hap2.chr10.bed ;;
  *) echo "unknown version: $V"; exit 1 ;;
esac
T_BAM=$B/bams/HG008-T.$S.bam
N_BAM=$B/bams/HG008-N-D.$S.bam
OUT=$P/results/dsa_pilot_chr10/calls/$V/$CALLER
mkdir -p $OUT
start=$SECONDS
if [ "$CALLER" = deepsomatic ]; then
  mkdir -p $OUT/tmp
  singularity exec -B /storage $P/containers/deepsomatic_1.8.0.sif run_deepsomatic \
    --model_type=WGS --ref=$REF --reads_tumor=$T_BAM --reads_normal=$N_BAM \
    --sample_name_tumor=HG008-T --sample_name_normal=HG008-N-D \
    --regions=$BED --num_shards=$T --intermediate_results_dir=$OUT/tmp \
    --output_vcf=$OUT/calls.vcf.gz
elif [ "$CALLER" = mutect2 ]; then
  $RUN gatk Mutect2 -R $REF -I $T_BAM -I $N_BAM -normal HG008-N-D -L $BED \
    --native-pair-hmm-threads $T -O $OUT/raw.vcf.gz
else
  echo "unknown caller: $CALLER"; exit 1
fi
echo "TIME: $CALLER $V took $(( (SECONDS - start) / 60 )) minutes"
echo "done: $CALLER $V"
