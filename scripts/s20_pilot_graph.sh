#!/bin/bash
#SBATCH -p cpu
#SBATCH --job-name=pilot_graph
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH -o results/logs/%x_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
RUN="pixi run --frozen -e pipeline"
G=data/raw/graph/hprc-v2.1-mc-chm13/hprc-v2.1-mc-chm13.d46
READS=data/raw/reads/BCM_Illumina-WGS_20240313
IN=results/test_grch38
OUT=results/pilot_graph
mkdir -p "$IN" "$OUT"

# 0. Make the test reads if they are not here (same first 2M pairs as s06)
LINES=8000000
make_test() { (zcat "$READS/$1" | head -n $LINES || true) | gzip > "$IN/$2"; }
if [ ! -f $IN/T_R1.fastq.gz ]; then
  make_test HG008-T_fastqs/HJVHNDSX7-1-IDUDI0031_S12_L001_R1_001.fastq.gz T_R1.fastq.gz
  make_test HG008-T_fastqs/HJVHNDSX7-1-IDUDI0031_S12_L001_R2_001.fastq.gz T_R2.fastq.gz
  make_test HG008-N-D_fastqs/HV5TMDSX7-1-IDUDI0034_S1_L001_R1_001.fastq.gz N_R1.fastq.gz
  make_test HG008-N-D_fastqs/HV5TMDSX7-1-IDUDI0034_S1_L001_R2_001.fastq.gz N_R2.fastq.gz
fi
echo "test reads ready: $(zcat $IN/T_R1.fastq.gz | wc -l) lines in T_R1"

# 1. Who is in the graph? Is HG008 there?
$RUN vg paths -M -x $G.gbz | cut -f3 | sort -u > $OUT/graph_samples.txt
echo "samples in graph: $(wc -l < $OUT/graph_samples.txt)"
if grep -qiw "HG008" $OUT/graph_samples.txt; then echo "WARNING: HG008 IS in the graph"; else echo "OK: HG008 not in the graph"; fi

# 2-3. giraffe on the test reads, then sort, index, count
for s in T N; do
  if [ $s = T ]; then SM=HG008-T; else SM=HG008-N-D; fi
  $RUN vg giraffe -t 16 -Z $G.gbz -d $G.dist -m $G.shortread.withzip.min -z $G.shortread.zipcodes \
    -f $IN/${s}_R1.fastq.gz -f $IN/${s}_R2.fastq.gz -N $SM -R ${SM}_test --ref-paths metadata/chm13_ref_paths.txt -o BAM > $OUT/${s}.unsorted.bam
  $RUN samtools sort -@ 4 -o $OUT/${SM}_graph_test.bam $OUT/${s}.unsorted.bam
  $RUN samtools index $OUT/${SM}_graph_test.bam
  $RUN samtools flagstat $OUT/${SM}_graph_test.bam > $OUT/${SM}_graph_test.flagstat.txt
  rm $OUT/${s}.unsorted.bam
done

# 4. What are the chromosomes called?
$RUN samtools view -H $OUT/HG008-T_graph_test.bam | grep "^@SQ" | head -5 || true
echo "done: graph pilot in $OUT"
