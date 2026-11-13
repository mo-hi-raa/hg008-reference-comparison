# Check d for the DSA pilot: where do the wrong calls (vcfeval fp.vcf.gz) fall?
# Usage: python3 -I pilot_fp_check.py EVAL_DIR CHR RUN   (RUN e.g. V2_deepsomatic)
import sys, gzip, re, bisect
E, CHR, RUN = sys.argv[1], sys.argv[2], sys.argv[3]
H1, H2 = CHR + "_hap1", CHR + "_hap2"
paf = {H1: E + "/hap1_to_hap2.paf", H2: E + "/hap2_to_hap1.paf"}
OP = re.compile(r":(\d+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)")
diff, cov = {}, {}
for c in (H1, H2):
    d, iv = [], []
    for line in open(paf[c]):
        f = line.rstrip("\n").split("\t")
        if f[4] != "+" or "tp:A:P" not in f[12:]: continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if cs is None: continue
        q = int(f[2]); iv.append((q, int(f[3])))
        for m in OP.finditer(cs):
            if m.group(1): q += int(m.group(1))
            elif m.group(2): d.append(q); q += 1
            elif m.group(4): d.extend(range(q, q + len(m.group(4)))); q += len(m.group(4))
            else: d.append(q)
    diff[c] = sorted(d); cov[c] = sorted(iv)
def near(lst, p, w):
    i = bisect.bisect_left(lst, p - w); return i < len(lst) and lst[i] <= p + w
def covered(c, p):
    iv = cov[c]; i = bisect.bisect_right(iv, (p, float("inf"))) - 1
    return i >= 0 and iv[i][0] <= p < iv[i][1]
ids, untw = {}, {H1: [], H2: []}
for l in gzip.open(E + "/truth_twins.vcf.gz", "rt"):
    if l.startswith("#"): continue
    r = l.split("\t"); ids.setdefault(r[2], []).append((r[0], int(r[1]) - 1))
for k, v in ids.items():
    if len(v) == 1: untw[v[0][0]].append(v[0][1])
for c in untw: untw[c].sort()
fps = []
with gzip.open(E + "/vcfeval_" + RUN + "/fp.vcf.gz", "rt") as f:
    for l in f:
        if l.startswith("#"): continue
        r = l.split("\t"); fps.append((r[0], int(r[1]) - 1, len(r[3]) == 1 and len(r[4]) == 1))
pos = {H1: sorted(p for c, p, s in fps if c == H1), H2: sorted(p for c, p, s in fps if c == H2)}
n = {}
def add(k): n[k] = n.get(k, 0) + 1
for c, p, snv in fps:
    t = "SNV" if snv else "indel"; add(t + " total")
    if not covered(c, p): add(t + " a) no matching place on other copy")
    elif near(diff[c], p, 5): add(t + " b) at a germline hap1/hap2 difference")
    elif near(untw[c], p, 50): add(t + " c) near a truth variant without twin")
    else: add(t + " d) other")
    lst = pos[c]; i = bisect.bisect_left(lst, p - 500); j = bisect.bisect_right(lst, p + 500)
    if j - i >= 3: add(t + " (in a cluster of 3+ within 1 kb)")
for k in sorted(n): print("FPCHECK", CHR, RUN, k, n[k])
