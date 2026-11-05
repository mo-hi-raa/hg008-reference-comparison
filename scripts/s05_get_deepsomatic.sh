#!/usr/bin/env bash
# Get the DeepSomatic 1.8.0 container, the version used in the HG008 paper
# (bioRxiv 2026.05.01.722316, p. 53). Saves containers/deepsomatic_1.8.0.sif.
# Run from the project folder with: sbatch scripts/s05_get_deepsomatic.sh
#SBATCH --job-name=get_deepsomatic
#SBATCH --partition=ws-ia
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=results/logs/%x_%j.out
set -euo pipefail

mkdir -p containers tmp_singularity
export SINGULARITY_TMPDIR=$PWD/tmp_singularity
export SINGULARITY_CACHEDIR=$PWD/tmp_singularity

singularity pull containers/deepsomatic_1.8.0.sif docker://google/deepsomatic:1.8.0
singularity exec containers/deepsomatic_1.8.0.sif run_deepsomatic --version
rm -rf tmp_singularity
echo "done: containers/deepsomatic_1.8.0.sif"
