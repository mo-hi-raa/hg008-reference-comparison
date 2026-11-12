# Make haplotype-agnostic ("twin") truth for the DSA pilot (Way B, as SMaHT/Sohn 2026):
# each truth variant is copied to the matching position on the other haplotype,
# using the minimap2 hap-to-hap alignments (PAF with cs tags). A twin is made only if
# the whole REF lifts contiguously and the other haplotype has the same REF bases.
# Also lifts the benchmark BED to the other haplotype, so twins are inside the scored region.
# Usage: python3 -I pilot_lift.py EVAL_DIR CHR
import sys, re, bisect
E, CHR = sys.argv[1], sys.argv[2]
H1, H2 = CHR + "_hap1", CHR + "_hap2"
other = {H1: H2, H2: H1}
paf = {H1: E + "/hap1_to_hap2.paf", H2: E + "/hap2_to_hap1.paf"}
order = {H1: 0, H2: 1}

def fasta(path):
    with open(path) as f:
        return "".join(l.strip() for l in f if not l.startswith(">")).upper()
seq = {H1: fasta(E + "/" + H1 + ".fa"), H2: fasta(E + "/" + H2 + ".fa")}

OP = re.compile(r":(\d+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)")
def lift(path, positions):
    pos = sorted(set(positions)); out = {}
    if not pos: return out
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[4] != "+" or "tp:A:P" not in f[12:]: continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if cs is None: continue
        q, t = int(f[2]), int(f[7])
        if bisect.bisect_left(pos, q) == bisect.bisect_left(pos, int(f[3])): continue
        for m in OP.finditer(cs):
            if m.group(1):
                n = int(m.group(1))
                for p in pos[bisect.bisect_left(pos, q):bisect.bisect_left(pos, q + n)]:
                    out.setdefault(p, t + (p - q))
                q += n; t += n
            elif m.group(2):
                i = bisect.bisect_left(pos, q)
                if i < len(pos) and pos[i] == q: out.setdefault(q, t)
                q += 1; t += 1
            elif m.group(4): q += len(m.group(4))
            else: t += len(m.group(5))
    return out

hdr, recs = [], []
for l in open(E + "/truth." + CHR + ".vcf"):
    if l.startswith("##"): hdr.append(l)
    elif l.startswith("#"):
        c = l.rstrip("\n").split("\t"); hdr.append("\t".join(c[:9] + ["TRUTH"]) + "\n")
    else: recs.append(l.rstrip("\n").split("\t")[:10])
bed = [l.rstrip("\n").split("\t")[:3] for l in open(E + "/truth." + CHR + ".bed") if l.strip()]

need = {H1: [], H2: []}
for r in recs:
    p0 = int(r[1]) - 1; need[r[0]] += [p0, p0 + len(r[3]) - 1]
for c, s, e in bed: need[c] += [int(s), int(e) - 1]
maps = {c: lift(paf[c], need[c]) for c in (H1, H2)}

out, st = [], {}
def add(k): st[k] = st.get(k, 0) + 1
for r in recs:
    c, pos, ref, alt = r[0], int(r[1]), r[3], r[4]
    r[2] = c + ":" + str(pos) + ":" + ref + ":" + alt
    out.append(r); add(c + " variants")
    m = maps[c]; p0 = pos - 1; p1 = p0 + len(ref) - 1
    if p0 in m and p1 in m and m[p1] - m[p0] == p1 - p0:
        if seq[other[c]][m[p0]:m[p0] + len(ref)] == ref.upper():
            tw = list(r); tw[0] = other[c]; tw[1] = str(m[p0] + 1); out.append(tw); add(c + " twins made")
        else: add(c + " no twin: different bases on other copy")
    else: add(c + " no twin: position not on other copy")

lb = []
for c, s, e in bed:
    s, e = int(s), int(e); lb.append((c, s, e)); m = maps[c]
    if s in m and (e - 1) in m and m[e - 1] >= m[s]: lb.append((other[c], m[s], m[e - 1] + 1))
lb.sort(key=lambda x: (order[x[0]], x[1], x[2])); merged = []
for c, s, e in lb:
    if merged and merged[-1][0] == c and s <= merged[-1][2]: merged[-1][2] = max(merged[-1][2], e)
    else: merged.append([c, s, e])
with open(E + "/eval.bed", "w") as f:
    for c, s, e in merged: f.write(c + "\t" + str(s) + "\t" + str(e) + "\n")
out.sort(key=lambda r: (order[r[0]], int(r[1])))
with open(E + "/truth_twins.vcf", "w") as f:
    f.writelines(hdr)
    for r in out: f.write("\t".join(r) + "\n")
for k in sorted(st): print("twins:", CHR, k, st[k])
for h in (H1, H2): print("eval region:", h, sum(e - s for c, s, e in merged if c == h), "bases")
