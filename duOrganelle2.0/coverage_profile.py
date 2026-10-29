#!/usr/bin/env python3
"""
coverage_profile.py -- Stage A2: per-compartment coverage profile and Stage B
route decision for DuOrganelle.

Reads the outputs of ONE BBMap pass of the raw reads against a combined
reference whose headers are prefixed "pt__" (plastid) and "mt__"
(mitochondrial), and reports:

  COV_PT / COV_MT   length-weighted median of BBMap's per-sequence
                    Median_fold, per compartment (medians, not means, so
                    MTPTs, NUMTs and multicopy repeats do not inflate them)
  ORG_FRAC          fraction of input reads mapped to any organellar
                    reference (dominated by the plastome)
  INSERT_MEDIAN     median insert size (BBMap ihist), paired data only
  SAMPLERATE        fraction of plastid-classified read pairs to KEEP so that
                    plastid depth ends up ~TARGET_RATIO x mitochondrial depth:
                        samplerate = target_ratio * COV_MT / COV_PT  (capped at 1)
  STAGEB_ROUTE      "recruit" (total-DNA library: GetOrganelle seed-and-extend
                    is selective and worth using) or "graph" (organelle-
                    enriched library: recruitment captures nearly everything,
                    so assemble globally with SPAdes after partial plastid
                    depletion and extract the mt subgraph instead)

WHY THIS STEP EXISTS: on the Fuchsia magellanica enriched library,
get_organelle_from_reads.py -F embplant_mt grew its word dictionary from 2 M to
210 M words in three rounds (reads recruited: ~22 % of all pairs by round 3) --
recruitment was selecting almost nothing while costing hours and RAM. The
decision between recruitment and a global-graph route can be made up front
from the organellar read fraction instead of discovered by a runaway run.

CAVEAT: ORG_FRAC is only meaningful if the pt reference is a (near-)complete
plastome (-r_pt, or Stage A's pt accumulation). With seed genes only, ORG_FRAC
is grossly underestimated and the route will default to "recruit"; force the
route with --route if you know the library type.

Writes KEY=value lines (numbers/words only) to --out-env for the shell to
source, and a human-readable summary to stdout.
"""
import argparse
import re
import sys


def weighted_median(pairs):
    """pairs: list of (value, weight). Returns weighted median or None."""
    pairs = [(v, w) for v, w in pairs if w > 0]
    if not pairs:
        return None
    pairs.sort()
    total = sum(w for _, w in pairs)
    acc = 0.0
    for v, w in pairs:
        acc += w
        if acc >= total / 2:
            return v
    return pairs[-1][0]


def parse_covstats(path):
    comp = {"pt": [], "mt": []}
    mapped_reads = 0
    with open(path) as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header = line.lstrip("#").split("\t")
                continue
            if header is None:
                continue
            row = dict(zip(header, line.split("\t")))
            name = row.get("ID", "").split()[0]
            try:
                length = int(row["Length"])
                med = float(row.get("Median_fold", row.get("Avg_fold", 0)))
                mapped_reads += int(row.get("Plus_reads", 0)) + int(row.get("Minus_reads", 0))
            except (KeyError, ValueError):
                continue
            for c in ("pt", "mt"):
                if name.startswith(c + "__"):
                    comp[c].append((med, length))
    return weighted_median(comp["pt"]), weighted_median(comp["mt"]), mapped_reads


def parse_ihist(path):
    try:
        with open(path) as f:
            for line in f:
                if line.startswith("#Median"):
                    return int(float(line.split("\t")[1]))
    except (FileNotFoundError, IndexError, ValueError):
        pass
    return None


def parse_reads_used(path):
    try:
        with open(path, errors="replace") as f:
            for line in f:
                m = re.match(r"\s*Reads Used:\s+(\d+)", line)
                if m:
                    return int(m.group(1))
    except FileNotFoundError:
        pass
    return None


def parse_mapped_pct(path):
    """Mean of BBMap's per-read-end 'mapped: xx.x%' lines (fraction of reads
    actually processed, so it stays valid when BBMap subsamples)."""
    vals = []
    try:
        with open(path, errors="replace") as f:
            for line in f:
                m = re.match(r"\s*mapped:\s+([\d.]+)%", line)
                if m:
                    vals.append(float(m.group(1)))
    except FileNotFoundError:
        pass
    return (sum(vals) / len(vals) / 100.0) if vals else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--covstats", required=True)
    ap.add_argument("--ihist", default="")
    ap.add_argument("--bbmap-log", default="", help="BBMap stderr, for 'Reads Used:' (total input reads)")
    ap.add_argument("--target-ratio", type=float, default=2.5,
                    help="target plastid:mitochondrial depth after depletion (default 2.5)")
    ap.add_argument("--enrich-threshold", type=float, default=0.5,
                    help="ORG_FRAC at or above which the library is treated as organelle-enriched (default 0.5)")
    ap.add_argument("--route", choices=["auto", "recruit", "graph"], default="auto")
    ap.add_argument("--sample-fraction", type=float, default=1.0,
                    help="fraction of reads BBMap mapped (its samplerate=); depths are divided by it")
    ap.add_argument("--out-env", required=True)
    a = ap.parse_args()

    cov_pt, cov_mt, mapped = parse_covstats(a.covstats)
    sf = a.sample_fraction if 0 < a.sample_fraction <= 1 else 1.0
    cov_pt = cov_pt / sf if cov_pt is not None else None
    cov_mt = cov_mt / sf if cov_mt is not None else None
    insert = parse_ihist(a.ihist) if a.ihist else None
    org_frac = parse_mapped_pct(a.bbmap_log) if a.bbmap_log else None
    if org_frac is None:
        total = parse_reads_used(a.bbmap_log) if a.bbmap_log else None
        org_frac = min(1.0, mapped / total) if (total and total > 0) else None

    if cov_pt and cov_mt and cov_mt > 0:
        samplerate = min(1.0, a.target_ratio * cov_mt / cov_pt)
    else:
        samplerate = 1.0

    if a.route != "auto":
        route, why = a.route, "forced with --route"
    elif org_frac is None:
        route, why = "recruit", "organellar read fraction could not be computed; defaulting to recruitment"
    elif org_frac >= a.enrich_threshold:
        route, why = "graph", f"organellar fraction {org_frac:.2f} >= {a.enrich_threshold} (enriched library)"
    else:
        route, why = "recruit", f"organellar fraction {org_frac:.2f} < {a.enrich_threshold} (total-DNA library)"

    fmt = lambda x, d=1: ("NA" if x is None else f"{x:.{d}f}")
    with open(a.out_env, "w") as f:
        f.write(f"COV_PT={fmt(cov_pt) if cov_pt is not None else 0}\n")
        f.write(f"COV_MT={fmt(cov_mt) if cov_mt is not None else 0}\n")
        f.write(f"ORG_FRAC={fmt(org_frac, 4) if org_frac is not None else 0}\n")
        f.write(f"INSERT_MEDIAN={insert if insert is not None else 0}\n")
        f.write(f"SAMPLERATE={samplerate:.5f}\n")
        f.write(f"STAGEB_ROUTE={route}\n")

    ratio = (cov_pt / cov_mt) if (cov_pt and cov_mt) else None
    print(f"[profile] median depth pt={fmt(cov_pt)}x  mt={fmt(cov_mt)}x  pt:mt={fmt(ratio)}")
    print(f"[profile] organellar read fraction={fmt(org_frac, 3)}  median insert={insert if insert else 'NA'} bp")
    print(f"[profile] plastid-pair samplerate={samplerate:.4f} (target pt:mt = {a.target_ratio})")
    print(f"[profile] Stage B route: {route.upper()} -- {why}")
    if cov_mt is None or cov_mt == 0:
        print("[profile] WARNING: no mitochondrial depth estimated (no reads on mt__ references); "
              "depletion will be skipped (samplerate=1).")


if __name__ == "__main__":
    main()
