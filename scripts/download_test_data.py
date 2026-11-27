"""Step 1 of the CBIO 7501 capstone: fetch the test data from NIST.

Only our stretch of chr22 (the region in config/config.yaml) is taken
from NIST's public files, so the whole genome is never downloaded:
  - tumour and normal reads, from NIST's GRCh38 BAMs, saved as FASTQ
  - the matching piece of the GRCh38 reference
  - the known mutations (benchmark VCF) and trusted regions (BED)
Each file is then checked against metadata/capstone_checksums.sha256.

Run: pixi run -e capstone python scripts/download_test_data.py
"""

import gzip
import hashlib
import subprocess
import sys
import time
from pathlib import Path

import yaml

# ---- Settings, read from config/config.yaml ----
ROOT = Path(__file__).resolve().parent.parent      # the project folder
with open(ROOT / "config" / "config.yaml") as handle:
    CONFIG = yaml.safe_load(handle)
REGION = CONFIG["region"]                          # chr22:34127223-34838003
NIST = CONFIG["nist"]                              # the five NIST addresses
OUT = ROOT / CONFIG["data_dir"]                    # data/raw/test
CHECKSUMS = ROOT / "metadata" / "capstone_checksums.sha256"
TRIES = 3                                          # attempts before giving up

FILES = [
    "HG008-T_R1.fq.gz", "HG008-T_R2.fq.gz",        # tumour reads
    "HG008-N-D_R1.fq.gz", "HG008-N-D_R2.fq.gz",    # normal reads
    "reference.fa",                                # reference piece
    "known_mutations.vcf",                         # NIST benchmark in the stretch
    "trusted_regions.bed",                         # NIST trusted regions in the stretch
]


def run(command):
    """Run a shell command inside OUT; try again if NIST does not answer."""
    for attempt in range(1, TRIES + 1):
        result = subprocess.run(["bash", "-c", "set -euo pipefail; " + command], cwd=OUT)
        if result.returncode == 0:
            return
        print(f"Attempt {attempt} of {TRIES} failed.", file=sys.stderr)
        if attempt < TRIES:
            time.sleep(10)
    sys.exit("Could not fetch the data from NIST. The server may be down; try again later.")


def fetch_reads(url, sample):
    """Reads in REGION from one NIST BAM, saved as a pair of FASTQ files."""
    run(
        f"samtools view -u {url} {REGION}"
        f" | samtools sort -n -u -T sort_{sample} -"
        f" | samtools fastq -n -0 /dev/null -s /dev/null"
        f" -1 {sample}_R1.fq.gz -2 {sample}_R2.fq.gz -"
    )


def fetch_reference():
    """The piece of the GRCh38 reference that matches REGION."""
    run(f"samtools faidx {NIST['reference']} {REGION} > reference.fa")


def fetch_known_mutations():
    """NIST's known somatic mutations inside REGION (with the VCF header)."""
    run(f"tabix -h {NIST['benchmark_vcf']} {REGION} > known_mutations.vcf")


def fetch_trusted_regions():
    """NIST's trusted regions that overlap REGION."""
    chrom, span = REGION.split(":")
    start, end = span.split("-")
    run(
        f"curl -fsS {NIST['benchmark_bed']}"
        f" | awk -v c={chrom} -v s={start} -v e={end} '$1==c && $2<e && $3>s-1'"
        f" > trusted_regions.bed"
    )


def content_checksum(path):
    """sha256 of a file's content. .gz files are unzipped first, so the
    result does not depend on how a computer compresses."""
    opener = gzip.open if path.suffix == ".gz" else open
    digest = hashlib.sha256()
    with opener(path, "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def check_checksums():
    """Compare every file with metadata/capstone_checksums.sha256.
    If that file does not exist yet (our own first run), write it."""
    found = {name: content_checksum(OUT / name) for name in FILES}
    if not CHECKSUMS.exists():
        with open(CHECKSUMS, "w") as handle:
            for name, checksum in found.items():
                handle.write(f"{checksum}  {name}\n")
        print(f"First run: checksums saved to {CHECKSUMS.relative_to(ROOT)}")
        return
    expected = {}
    for line in CHECKSUMS.read_text().splitlines():
        checksum, name = line.split()
        expected[name] = checksum
    wrong = [name for name in FILES if found[name] != expected.get(name)]
    if wrong:
        sys.exit("Checksum does not match for: " + ", ".join(wrong))
    print("All checksums match.")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name in FILES:                     # allow overwriting a previous fetch
        if (OUT / name).exists():
            (OUT / name).chmod(0o644)
    print("Fetching tumour reads...")
    fetch_reads(NIST["tumour_bam"], "HG008-T")
    print("Fetching normal reads...")
    fetch_reads(NIST["normal_bam"], "HG008-N-D")
    print("Fetching the reference piece...")
    fetch_reference()
    print("Fetching known mutations and trusted regions...")
    fetch_known_mutations()
    fetch_trusted_regions()
    check_checksums()
    for name in FILES:                     # raw data is read-only (lesson 02)
        (OUT / name).chmod(0o444)
    print(f"Done: test data is in {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
