"""Step 6 of the CBIO 7501 capstone: compare our PASS calls with NIST's
known mutations, inside NIST's trusted regions.

Our calls are numbered from the start of the reference piece; NIST's are
numbered from the start of chr22. The piece starts at chr22:34,127,223, so
position on chr22 = position on the piece + 34,127,222.

A call matches a known mutation when the position, the reference letters
and the new letters are all the same.

Usage: python scripts/compare_calls.py CALLS.vcf.gz KNOWN.vcf TRUSTED.bed REGION OUT_DIR
"""

import gzip
import sys
from pathlib import Path

calls_vcf, known_vcf, trusted_bed, region, out_dir = sys.argv[1:6]

chrom, span = region.split(":")                # chr22, 34127223-34838003
offset = int(span.split("-")[0]) - 1           # 34127222


def read_vcf(path, shift=0):
    """Set of (position on chr22, ref, alt), one entry per new allele."""
    opener = gzip.open if path.endswith(".gz") else open
    variants = set()
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            pos, ref, alts = int(fields[1]) + shift, fields[3], fields[4]
            for alt in alts.split(","):
                variants.add((pos, ref, alt))
    return variants


def read_trusted(path):
    """Trusted regions as (start, end); BED counts start from 0."""
    regions = []
    with open(path) as handle:
        for line in handle:
            fields = line.split("\t")
            regions.append((int(fields[1]), int(fields[2])))
    return regions


def is_trusted(pos, regions):
    return any(start < pos <= end for start, end in regions)


all_calls = read_vcf(calls_vcf, shift=offset)
known = read_vcf(known_vcf)
trusted = read_trusted(trusted_bed)

calls = {v for v in all_calls if is_trusted(v[0], trusted)}
known = {v for v in known if is_trusted(v[0], trusted)}
outside = len(all_calls) - len(calls)

found = calls & known
missed = known - calls
extra = calls - known

out = Path(out_dir)
out.mkdir(parents=True, exist_ok=True)

with open(out / "comparison.tsv", "w") as table:
    table.write("chrom\tpos\tref\talt\tstatus\n")
    for status, group in [("found", found), ("missed", missed), ("extra", extra)]:
        for pos, ref, alt in sorted(group):
            table.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{status}\n")

recall = len(found) / len(known) if known else 0
precision = len(found) / len(calls) if calls else 0
with open(out / "summary.tsv", "w") as summary:
    summary.write("known\tcalls\tfound\tmissed\textra\tcalls_outside_trusted\trecall\tprecision\n")
    summary.write(f"{len(known)}\t{len(calls)}\t{len(found)}\t{len(missed)}\t{len(extra)}"
                  f"\t{outside}\t{recall:.3f}\t{precision:.3f}\n")

print(f"Known: {len(known)}  Calls: {len(calls)}  Found: {len(found)}  "
      f"Missed: {len(missed)}  Extra: {len(extra)}  Outside trusted: {outside}")
