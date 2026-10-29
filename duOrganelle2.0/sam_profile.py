#!/usr/bin/env python3
"""
sam_profile.py -- Stage A2 helper: turn a minimap2 SAM stream (paired or
single-end reads mapped to the small pt__/mt__ profiling reference) into the
same three files coverage_profile.py reads from BBMap:

  covstats  BBMap-covstats-like table: #ID Length Plus_reads Minus_reads
            Avg_fold Median_fold (per reference sequence)
  ihist     '#Median <insert>' (median TLEN of proper pairs, read 1 only)
  log       'Reads Used:' and 'mapped: xx.xx%' lines

Why: BBMap 40.05 stalled indefinitely in index construction for this tiny
reference on svrlinux01 (two Stage A2 runs: >19 h and >10 min stuck at
"Indexing threads started"). minimap2 -ax sr maps the same reads in minutes.

Counting rules: primary records only (no 0x100 secondary / 0x800
supplementary); a read counts as mapped if aligned with identity >= --minid
(1 - NM / aligned length), mimicking BBMap's minid=0.90. Depth = per-base
count of mapped reads over reference-consuming CIGAR ops (M, =, X, D).
Reads with equally good hits (e.g. plastid inverted repeat) are placed by
minimap2 on one of them, comparable to BBMap ambiguous=random.

Usage: minimap2 -ax sr -t 16 --secondary=no ref.fa R1.fq.gz R2.fq.gz \
         | sam_profile.py --ref ref.fa --covstats c.txt --ihist i.txt --log l.txt
"""
import argparse
import re
import sys

CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def ref_lengths(path):
    lens, name, n = {}, None, 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--covstats", required=True)
    ap.add_argument("--ihist", required=True)
    ap.add_argument("--log", required=True)
    ap.add_argument("--minid", type=float, default=0.90)
    a = ap.parse_args()

    lens = ref_lengths(a.ref)
    diff = {r: [0] * (L + 2) for r, L in lens.items()}
    plus = {r: 0 for r in lens}
    minus = {r: 0 for r in lens}
    total = mapped = 0
    tlens = []

    for line in sys.stdin:
        if line[0] == "@":
            continue
        f = line.split("\t", 11)
        flag = int(f[1])
        if flag & 0x900:
            continue
        total += 1
        if flag & 0x4 or f[2] not in diff or f[5] == "*":
            continue
        ops = CIG.findall(f[5])
        alen = sum(int(n) for n, o in ops if o in "MI=X")
        refspan = sum(int(n) for n, o in ops if o in "MD=XN")
        nm = 0
        tags = f[11] if len(f) > 11 else ""
        m = re.search(r"\tNM:i:(\d+)|^NM:i:(\d+)", tags)
        if m:
            nm = int(m.group(1) or m.group(2))
        if alen == 0 or 1 - nm / alen < a.minid:
            continue
        mapped += 1
        r, pos = f[2], int(f[3])
        d = diff[r]
        start = max(1, pos)
        end = min(lens[r], pos + refspan - 1)
        if end >= start:
            d[start] += 1
            d[end + 1] -= 1
        if flag & 0x10:
            minus[r] += 1
        else:
            plus[r] += 1
        if flag & 0x2 and flag & 0x40:
            t = abs(int(f[8]))
            if t > 0:
                tlens.append(t)

    with open(a.covstats, "w") as out:
        out.write("#ID\tLength\tPlus_reads\tMinus_reads\tAvg_fold\tMedian_fold\n")
        for r, L in lens.items():
            d, cov, run = diff[r], [], 0
            for i in range(1, L + 1):
                run += d[i]
                cov.append(run)
            avg = sum(cov) / L if L else 0
            med = sorted(cov)[L // 2] if L else 0
            out.write(f"{r}\t{L}\t{plus[r]}\t{minus[r]}\t{avg:.2f}\t{med}\n")
    with open(a.ihist, "w") as out:
        if tlens:
            tlens.sort()
            out.write(f"#Median\t{tlens[len(tlens) // 2]}\n")
    with open(a.log, "w") as out:
        out.write(f"Reads Used:\t{total}\n")
        pct = 100.0 * mapped / total if total else 0.0
        out.write(f"mapped:\t{pct:.4f}%\t{mapped}\n")
    print(f"[sam_profile] {total:,} primary reads, {mapped:,} mapped "
          f"({100.0 * mapped / max(total, 1):.2f}%), {len(tlens):,} proper pairs for insert size",
          file=sys.stderr)


if __name__ == "__main__":
    main()
