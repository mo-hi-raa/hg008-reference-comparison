#!/usr/bin/env bash
# PIPELINE TEST (C3c): call and score somatic variants on chr22 of the
# 2M-read-pair test BAMs. Scores will be near zero (tiny data); this only
# checks that every command runs and every file is valid.
# Run on Karakoram from the project folder: sbatch scripts/s10_test_call_score_chr22.sh
#SBATCH --job-name=test_call_chr22
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail

P=$(pwd -P)
RUN="pixi run --frozen -e pipeline"
REF=$P/data/processed/reference/grch38/GRCh38_GIABv3.fasta
SDF=$P/data/processed/reference/grch38/GRCh38_GIABv3.sdf
T_BAM=$P/results/test_grch38/HG008-T_test_chr22.bam
N_BAM=$P/results/test_grch38/HG008-N-D_test_chr22.bam
TRUTH=$P/data/raw/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz
TRUTH_BED=$P/data/raw/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed
TRUTH_SAMPLE=dipHG8Tv3.4chr19fixvGRCh38z200k.hap1.bam
OUT=$P/results/test_grch38/calls_chr22
mkdir -p $OUT/deepsomatic_tmp

# 1. DeepSomatic (inside its container)
singularity exec -B /storage containers/deepsomatic_1.8.0.sif run_deepsomatic \
  --model_type=WGS --ref=$REF \
  --reads_tumor=$T_BAM --reads_normal=$N_BAM \
  --sample_name_tumor=HG008-T --sample_name_normal=HG008-N-D \
  --regions=chr22 --num_shards=16 \
  --intermediate_results_dir=$OUT/deepsomatic_tmp \
  --output_vcf=$OUT/deepsomatic.vcf.gz

# 2. Mutect2, then its filtering step
$RUN gatk Mutect2 -R $REF -I $T_BAM -I $N_BAM -normal HG008-N-D -L chr22 -O $OUT/mutect2_raw.vcf.gz
$RUN gatk FilterMutectCalls -R $REF -V $OUT/mutect2_raw.vcf.gz -O $OUT/mutect2.vcf.gz

# 3. Benchmark regions on chr22 only
awk '$1=="chr22"' $TRUTH_BED > $OUT/benchmark_chr22.bed
echo "chr22 benchmark regions: $(wc -l < $OUT/benchmark_chr22.bed)"
echo "chr22 truth variants: $(zcat $TRUTH | grep -c '^chr22' || true)"
echo "first truth records:"; zcat $TRUTH | grep -v '^#' | head -3

# 4. RTG's reference format (made once, reused later)
[ -d $SDF ] || $RUN rtg format -o $SDF $REF

# 5. Score each caller with vcfeval and aardvark
for caller in deepsomatic mutect2; do
  $RUN rtg vcfeval -b $TRUTH -c $OUT/$caller.vcf.gz -t $SDF \
    --bed-regions $OUT/benchmark_chr22.bed --sample $TRUTH_SAMPLE,HG008-T \
    --squash-ploidy --threads 16 -o $OUT/vcfeval_$caller
  $RUN aardvark compare -r $REF -t $TRUTH -q $OUT/$caller.vcf.gz \
    -b $OUT/benchmark_chr22.bed --truth-sample $TRUTH_SAMPLE --query-sample HG008-T \
    --threads 16 -o $OUT/aardvark_$caller
done

# 6. Show the results
for caller in deepsomatic mutect2; do
  echo "== $caller: calls in VCF: $(zcat $OUT/$caller.vcf.gz | grep -vc '^#' || true)"
  echo "== $caller: vcfeval"; cat $OUT/vcfeval_$caller/summary.txt
  echo "== $caller: aardvark"; ls $OUT/aardvark_$caller
done
echo "done: C3c test in $OUT"
