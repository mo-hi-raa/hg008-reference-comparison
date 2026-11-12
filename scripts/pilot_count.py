# Count DSA pilot results from one vcfeval output, Way B: a truth variant counts as found
# if it or its twin was matched; found on both copies counts once ("echo").
# Usage: python3 -I pilot_count.py VCFEVAL_DIR "label"
import sys, gzip
d, label = sys.argv[1], sys.argv[2]
def recs(fn):
    with gzip.open(d + "/" + fn, "rt") as f:
        for l in f:
            if not l.startswith("#"): yield l.split("\t")
def cls(ref, alt):
    return "SNV" if len(ref) == 1 and all(len(a) == 1 for a in alt.split(",")) else "indel"
found, missed = {}, set()
for r in recs("tp-baseline.vcf.gz"): found[r[2]] = found.get(r[2], 0) + 1
for r in recs("fn.vcf.gz"): missed.add(r[2])
fp = [cls(r[3], r[4]) for r in recs("fp.vcf.gz")]
for c in ("SNV", "indel"):
    f = [k for k in found if cls(k.split(":")[2], k.split(":")[3]) == c]
    m = [k for k in missed if k not in found and cls(k.split(":")[2], k.split(":")[3]) == c]
    w = fp.count(c); echo = sum(1 for k in f if found[k] > 1)
    rec = len(f) / (len(f) + len(m)) if f or m else 0
    pre = len(f) / (len(f) + w) if f or w else 0
    print("RESULT %s %s: found %d, missed %d, recall %.3f | wrong calls %d, precision %.3f | found on both copies %d"
          % (label, c, len(f), len(m), rec, w, pre, echo))
