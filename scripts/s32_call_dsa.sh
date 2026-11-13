#!/usr/bin/env bash
# Full DSA run: call one chromosome (chr1-22, X) with DeepSomatic or Mutect2, same settings as s15/s27.
# combined_mq60 = V2 (both copies called together); hap1 / hap2 = V3 (joined before scoring).
# Usage: sbatch --array=1-23%2 scripts/s32_call_dsa.sh combined_mq60 deepsomatic
#SBATCH --job-name=call_dsa
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/%x_%A_%a.out
set -euo pipefail
ulimit -n $(ulimit -Hn)
LEVEL=$1; CALLER=$2
CHROMS=(chr{1..22} chrX)
CHR=${CHROMS[$((SLURM_ARRAY_TASK_ID - 1))]}
P=$(pwd -P)
RUN="pixi run --frozen -e pipeline"
T=$SLURM_CPUS_PER_TASK
D=$P/data/processed/reference/dsa_v6.3
case $LEVEL in
  combined_mq60) REF=$D/HG008N_v6.3.fasta;      CONTIGS="${CHR}_hap1 ${CHR}_hap2" ;;
  hap1)          REF=$D/HG008N_v6.3.hap1.fasta; CONTIGS="${CHR}_hap1" ;;
  hap2)          REF=$D/HG008N_v6.3.hap2.fasta; CONTIGS="${CHR}_hap2" ;;
  *) echo "unknown reference: $LEVEL"; exit 1 ;;
esac
A=$P/data/processed/alignments/dsa_$LEVEL
T_BAM=$A/HG008-T.dsa_$LEVEL.bam; N_BAM=$A/HG008-N-D.dsa_$LEVEL.bam
OUT=$P/results/calls_dsa_$LEVEL/$CALLER/by_chrom; mkdir -p $OUT
BED=$OUT/$CHR.bed
awk -v c="$CONTIGS" 'BEGIN{n=split(c,a," "); for(i=1;i<=n;i++) w[a[i]]=1} ($1 in w){print $1"\t0\t"$2}' $REF.fai > $BED
start=$SECONDS
if [ "$CALLER" = deepsomatic ]; then
  mkdir -p $OUT/tmp_$CHR
  singularity exec -B /storage $P/containers/deepsomatic_1.8.0.sif run_deepsomatic \
    --model_type=WGS --ref=$REF --reads_tumor=$T_BAM --reads_normal=$N_BAM \
    --sample_name_tumor=HG008-T --sample_name_normal=HG008-N-D \
    --regions=$BED --num_shards=$T --intermediate_results_dir=$OUT/tmp_$CHR \
    --output_vcf=$OUT/$CHR.vcf.gz
elif [ "$CALLER" = mutect2 ]; then
  $RUN gatk Mutect2 -R $REF -I $T_BAM -I $N_BAM -normal HG008-N-D -L $BED \
    --native-pair-hmm-threads $T -O $OUT/$CHR.raw.vcf.gz
else
  echo "unknown caller: $CALLER"; exit 1
fi
echo "TIME: $CALLER $LEVEL $CHR took $(( (SECONDS - start) / 60 )) minutes"
echo "done: $CALLER $LEVEL $CHR"
