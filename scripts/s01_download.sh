#!/usr/bin/env bash
# Download every raw file listed in metadata/sources.tsv into data/raw/.
# Run from the project folder: bash scripts/s01_download.sh
set -euo pipefail

mkdir -p data/raw

tail -n +2 metadata/sources.tsv | while IFS=$'\t' read -r file description source url checksum; do
    if [ -s "data/raw/$file" ]; then
        echo "already have $file"
        continue
    fi
    echo "downloading $file"
    curl -fL --retry 3 -o "data/raw/$file.part" "$url"
    mv "data/raw/$file.part" "data/raw/$file"
done

echo "done: files are in data/raw/"
