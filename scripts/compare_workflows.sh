#!/usr/bin/env bash
# Lesson 10 check: Snakemake and Nextflow must give exactly the same PASS calls
# and the same comparison summary. Headers are skipped (they hold dates and paths).
set -euo pipefail

snake=results/capstone/calls/pass.vcf.gz
nf=results/capstone_nextflow/pass.vcf.gz

if diff <(gzip -dc "$snake" | grep -v '^#') <(gzip -dc "$nf" | grep -v '^#') \
   && diff results/capstone/comparison/summary.tsv results/capstone_nextflow/comparison/summary.tsv; then
    echo "Same PASS calls and same summary in Snakemake and Nextflow."
else
    echo "Snakemake and Nextflow results differ." >&2
    exit 1
fi
