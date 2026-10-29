#!/usr/bin/env python3
"""
graph_depth_screen.py -- Stage C: per-node depth + cross-organelle screening
and classification of an organelle assembly graph (.gfa).

For every node of the graph of ONE organelle (--organelle mt|pt) this script
reports length, depth, connectivity to the main component, how much of it is
already known (Stage A accumulation) and how much matches the OTHER
organelle (--anti-reference), and then CLASSIFIES it by depth + homology.

Classes, mt graph (anti-reference = plastome):
  mt              single-copy mitochondrial (depth ~ 1 x D)
  mt_repeat       depth ~ n x D (n >= 2); `copies` = n
  mtpt            plastid homology at MITOCHONDRIAL depth, linked to >= 1
                  non-plastid node -> plastid-derived insert (MTPT): keep
  mtpt_isolated   as above but no non-plastid neighbour: inspect manually
  plastid         plastid homology at PLASTID depth -> contamination: exclude
  low_depth       depth < --low-ratio x D: NUMT / nuclear repeat / noise
Classes, pt graph (anti-reference = mitochondrial set):
  pt, pt_repeat (IR expected at ~2 x D), pt_shared_with_mt (mt homology at
  plastid depth: genuine plastid sequence that also exists as an MTPT in the
  mitogenome: keep), mitochondrial (mt homology at mt depth -> contamination),
  low_depth.
Flags: `copy_uncertain` (depth ratio >= 1.3 but not near an integer: short
node with few k-mers, substoichiometric repeat, or inflated depth) and
`partial_other_depth_inflated` (5-50 % of the node is other-organelle
sequence, e.g. an MTPT identical to the plastome, whose reads inflate the
node's mean depth -- seen on Fuchsia magellanica node 7344871, 1.43x with a
3.4 kb plastid-IR-identical MTPT, which GetOrganelle scored as 2 copies).
`span` (repeats only): "reads" (< read length) / "pairs" (< median insert) /
"exceeds_insert" -> the repeat cannot be resolved with these short reads;
alternative paths through it should be reported, not forced into one circle.

Baseline depth D = length-weighted median depth of nodes >=
--min-baseline-len bp with < 20 % other-organelle homology (falls back to all
such nodes). The plastid/mito depth comparison uses D_other = length-weighted
median depth of nodes >= 80 % other-organelle homology, when such nodes exist.

Outputs in --outdir:
  repeat_candidates.tsv  original columns (unchanged order) + new columns
  bandage_colours.csv    load in Bandage via File > Load CSV data
  summary on stdout      incl. estimated genome size = sum(length x copies)
                         over target-organelle classes

CHANGES (Oct 2026, merged with the former classify_graph_nodes.py):
  * depth is read from DP:f (GetOrganelle) or dp, KC:i/len (SPAdes) or RC:i/len;
    the previous version read RC:i only, so GetOrganelle graphs (DP:f only)
    got coverage 0 for every node and no repeat could ever be flagged.
  * fold_over_median is now relative to the robust baseline D above instead of
    the unweighted median of ALL nodes, which plastid-contaminant nodes and
    many tiny nodes could shift arbitrarily.
  * pct_matches_known_content now uses merged-interval coverage (as the
    other-organelle screen already did) instead of the single longest hit.
  * BLAST databases are built inside --outdir, not next to the reference
    file (which may live in a read-only directory).
"""
import argparse
import math
import os
import subprocess
from collections import defaultdict

COLOURS = {
    "mt": "#6b8fb3", "pt": "#6b8fb3",
    "mt_repeat": "#e08a1e", "pt_repeat": "#e08a1e",
    "mtpt": "#2e9d5b", "pt_shared_with_mt": "#2e9d5b",
    "mtpt_isolated": "#9acb5a",
    "plastid": "#d1342f", "mitochondrial": "#d1342f",
    "low_depth": "#8e6bb3",
}


def sh(cmd):
    """Run cmd, and if it fails, SHOW why (stdout+stderr). See the matching
    comment in redundancy_and_consensus.py's sh(): stderr is captured, never
    sent to /dev/null, and decoded leniently (errors='replace') only when
    building the error message, so a stray non-UTF-8 byte cannot crash the run."""
    result = subprocess.run(cmd, shell=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"\nCommand failed (exit {result.returncode}): {cmd}\n"
            f"--- stdout ---\n{result.stdout.decode('utf-8', errors='replace')}"
            f"--- stderr ---\n{result.stderr.decode('utf-8', errors='replace')}"
        )
    return result


def parse_gfa(gfa_path):
    nodes, edges = {}, []
    with open(gfa_path) as f:
        for line in f:
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t")
                name, seq = parts[1], parts[2]
                tags = {}
                for t in parts[3:]:
                    bits = t.split(":", 2)
                    if len(bits) == 3:
                        tags[bits[0]] = bits[2]
                ln = int(tags.get("LN", len(seq) if seq != "*" else 0))
                if "DP" in tags:
                    cov = float(tags["DP"])
                elif "dp" in tags:
                    cov = float(tags["dp"])
                elif "KC" in tags and ln:
                    cov = float(tags["KC"]) / ln
                elif "RC" in tags and ln:
                    cov = float(tags["RC"]) / ln
                else:
                    cov = 0.0
                nodes[name] = {"seq": seq, "len": ln, "cov": cov}
            elif line.startswith("L\t"):
                parts = line.rstrip("\n").split("\t")
                if parts[1] != parts[3]:
                    edges.append((parts[1], parts[3]))
    return nodes, edges


def connected_components(nodes, adj):
    seen, comps = set(), []
    for n in nodes:
        if n in seen:
            continue
        stack, comp = [n], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(adj[x] - comp)
        seen |= comp
        comps.append(comp)
    return comps


def merge_intervals(ivs):
    merged = []
    for lo, hi in sorted(ivs):
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def weighted_median(pairs):
    pairs = sorted((v, w) for v, w in pairs if w > 0 and v > 0)
    if not pairs:
        return None
    half, acc = sum(w for _, w in pairs) / 2, 0
    for v, w in pairs:
        acc += w
        if acc >= half:
            return v


def pct_covered_by(node_fasta, ref, tag, outdir, nodes, threads, identity, hits_override=None):
    """% of each node's length covered by merged BLASTn hits (>= identity) to ref."""
    if hits_override:
        tsv = hits_override
    else:
        db = os.path.join(outdir, f"{tag}_blastdb")
        tsv = os.path.join(outdir, f"nodes_vs_{tag}.tsv")
        sh(f"makeblastdb -in '{ref}' -dbtype nucl -out '{db}'")
        sh(f"blastn -query '{node_fasta}' -db '{db}' "
           f'-outfmt "6 qseqid sseqid pident length qstart qend" '
           f"-evalue 1e-20 -word_size 20 -num_threads {threads} -out '{tsv}'")
    ivs = defaultdict(list)
    with open(tsv) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 6 or float(p[2]) < identity:
                continue
            lo, hi = sorted([int(p[4]), int(p[5])])
            ivs[p[0]].append((lo, hi))
    out = {}
    for name, v in nodes.items():
        bp = sum(hi - lo + 1 for lo, hi in merge_intervals(ivs.get(name, [])))
        out[name] = min(100.0, 100.0 * bp / v["len"]) if v["len"] else 0.0
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--accumulated", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--organelle", choices=["mt", "pt"], default="mt",
                    help="organelle this graph represents (sets class names and the contamination logic)")
    ap.add_argument("--depth-fold-threshold", type=float, default=2.0,
                    help="flag_repeat_candidate: coverage >= this x baseline depth D")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--anti-reference", "--chloroplast", dest="anti_reference",
                    help="FASTA of the OTHER organelle (preferably a complete plastome for an mt graph). "
                         "Drives both the cross-organelle contamination flag and the classification; "
                         "without it, only depth-based classes are assigned.")
    ap.add_argument("--other-frac-threshold", "--plastid-frac-threshold", dest="other_frac_threshold",
                    type=float, default=0.5)
    ap.add_argument("--insert", type=int, default=0, help="median insert size (bp); 0 = unknown")
    ap.add_argument("--read-len", type=int, default=0)
    ap.add_argument("--baseline-depth", type=float, help="override baseline depth D")
    ap.add_argument("--low-ratio", type=float, default=0.3)
    ap.add_argument("--repeat-tol", type=float, default=0.35, help="max |r - round(r)| to call n copies")
    ap.add_argument("--other-ratio-min", type=float, default=5.0,
                    help="without an other-organelle depth baseline: an other-homologous node at >= this "
                         "x D (mt graph) or <= 1/this x D (pt graph) is called contamination")
    ap.add_argument("--min-baseline-len", type=int, default=500)
    ap.add_argument("--identity", type=float, default=90)
    ap.add_argument("--other-hits", help=argparse.SUPPRESS)   # precomputed outfmt6, for testing
    ap.add_argument("--known-hits", help=argparse.SUPPRESS)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    nodes, edges = parse_gfa(args.gfa)
    nodes = {k: v for k, v in nodes.items() if v["len"] > 0}
    if not nodes:
        print("No nodes found in GFA.")
        return
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    comps = sorted(connected_components(nodes, adj), key=lambda c: -sum(nodes[n]["len"] for n in c))
    main_comp = comps[0] if comps else set()

    node_fasta = os.path.join(args.outdir, "all_graph_nodes.fasta")
    with open(node_fasta, "w") as f:
        for name, v in nodes.items():
            if v["seq"] != "*":
                f.write(f">{name}\n{v['seq']}\n")

    known_cov = {}
    if args.known_hits or (os.path.exists(args.accumulated) and os.path.getsize(args.accumulated) > 0):
        known_cov = pct_covered_by(node_fasta, args.accumulated, "accumulated", args.outdir, nodes,
                                   args.threads, args.identity, args.known_hits)
    other_cov = {}
    if args.other_hits or (args.anti_reference and os.path.exists(args.anti_reference)
                           and os.path.getsize(args.anti_reference) > 0):
        other_cov = pct_covered_by(node_fasta, args.anti_reference, "other_organelle", args.outdir, nodes,
                                   args.threads, args.identity, args.other_hits)

    oth = lambda n: other_cov.get(n, 0.0) / 100.0
    tgt, other_name = args.organelle, ("plastid" if args.organelle == "mt" else "mitochondrial")
    D = args.baseline_depth or weighted_median(
        [(v["cov"], v["len"]) for k, v in nodes.items() if oth(k) < 0.2 and v["len"] >= args.min_baseline_len]
    ) or weighted_median([(v["cov"], v["len"]) for k, v in nodes.items() if oth(k) < 0.2]) or 0.0
    D_other = weighted_median([(v["cov"], v["len"]) for k, v in nodes.items() if oth(k) >= 0.8])
    if D <= 0:
        print("WARNING: no usable baseline depth (no DP/KC/RC tags, or every node matches the other "
              "organelle); classification skipped. Use --baseline-depth to force one.")

    rows = []
    for name, v in nodes.items():
        r = v["cov"] / D if D > 0 else 0.0
        fo = oth(name)
        cls, copies, span, flag = "", 1, "", ""
        if D > 0:
            if other_cov and fo >= args.other_frac_threshold:
                if D_other and (D_other > 1.5 * D or D_other < D / 1.5):
                    other_level = abs(math.log(max(v["cov"], 1e-9) / D_other)) < abs(math.log(max(v["cov"], 1e-9) / D))
                elif tgt == "mt":
                    other_level = r >= args.other_ratio_min
                else:
                    other_level = r <= 1.0 / args.other_ratio_min
                if other_level:
                    cls = other_name
                elif tgt == "mt":
                    has_tgt_nbr = any(oth(n) < args.other_frac_threshold for n in adj.get(name, ()) if n in nodes)
                    cls = "mtpt" if has_tgt_nbr else "mtpt_isolated"
                    copies = max(1, round(r))
                else:
                    cls = "pt_shared_with_mt"
                    copies = max(1, round(r))
            elif r < args.low_ratio:
                cls = "low_depth"
            else:
                n = round(r)
                if n >= 2 and abs(r - n) <= args.repeat_tol:
                    cls, copies = f"{tgt}_repeat", n
                elif r >= 1.3:
                    # between integer copy numbers: a repeat whose depth is
                    # biased (short node, few k-mers), substoichiometric, or
                    # inflated by reads of the other organelle -- report as a
                    # repeat with an explicit uncertainty flag, not as single-copy
                    cls, copies = f"{tgt}_repeat", max(2, n)
                    flag = "copy_uncertain"
                else:
                    cls = tgt
                if 0.05 <= fo < args.other_frac_threshold:
                    # >= 5 % of the node is other-organelle sequence (e.g. an MTPT
                    # identical to the plastome): its mean depth is inflated by
                    # plastid reads, so a high depth ratio is NOT evidence of copy number
                    flag = (flag + ";" if flag else "") + "partial_other_depth_inflated"
            if copies >= 2:
                if args.read_len and v["len"] < args.read_len:
                    span = "reads"
                elif args.insert and v["len"] < args.insert:
                    span = "pairs"
                elif args.insert:
                    span = "exceeds_insert"
                else:
                    span = "unknown_insert"
        rows.append((name, v, r, cls, copies, span, flag))

    out_tsv = os.path.join(args.outdir, "repeat_candidates.tsv")
    with open(out_tsv, "w") as out:
        out.write("node\tlength_bp\tcoverage\tfold_over_median\tconnected_to_main_graph\t"
                  "pct_matches_known_content\tpct_matches_other_organelle\tflag_repeat_candidate\t"
                  "flag_likely_cross_organelle_contamination\tclass\tcopies\tspan\tflag\n")
        for name, v, r, cls, copies, span, flag in rows:
            connected = name in main_comp
            pk, po = known_cov.get(name), other_cov.get(name)
            flag_rep = (r >= args.depth_fold_threshold) or (not connected and v["len"] > 1000)
            flag_other = cls == other_name if cls else (po is not None and po >= args.other_frac_threshold * 100)
            out.write(f"{name}\t{v['len']}\t{v['cov']:.2f}\t{r:.2f}\t{connected}\t"
                      f"{'' if pk is None else f'{pk:.1f}'}\t{'' if po is None else f'{po:.1f}'}\t"
                      f"{flag_rep}\t{flag_other}\t{cls}\t{copies}\t{span}\t{flag}\n")
            if flag_other:
                print(f"WARNING: node {name} ({v['len']} bp, coverage {v['cov']:.2f}, {r:.2f}x baseline) is "
                      f"{po:.1f}% {other_name} at {other_name}-level depth -- likely cross-organelle "
                      f"contamination. Recommend excluding it before annotation.")

    # curated graph: drop nodes classified as the OTHER organelle (e.g. plastid
    # edges that entered an mt graph at MTPT boundaries after plastid depletion)
    excluded = [name for name, v, r, cls, *_ in rows if cls == other_name]
    with open(os.path.join(args.outdir, "excluded_nodes.txt"), "w") as f:
        f.write("\n".join(excluded) + ("\n" if excluded else ""))
    with open(args.gfa) as fin, open(os.path.join(args.outdir, "curated_graph.gfa"), "w") as fo:
        for line in fin:
            p = line.split("\t")
            if p[0] == "S" and p[1] in excluded:
                continue
            if p[0] == "L" and (p[1] in excluded or p[3] in excluded):
                continue
            fo.write(line)
    if excluded:
        print(f"[{tgt}] excluded {len(excluded)} {other_name} node(s) from curated_graph.gfa: {', '.join(excluded)}")

    with open(os.path.join(args.outdir, "bandage_colours.csv"), "w") as f:
        f.write("Name,Colour,Class,Copies\n")
        for name, v, r, cls, copies, span, flag in rows:
            if cls:
                f.write(f"{name},{COLOURS[cls]},{cls},{copies}\n")

    print(f"[{tgt}] baseline depth D={D:.2f}" + (f", {other_name} baseline={D_other:.2f}" if D_other else ""))
    if D > 0:
        summ = defaultdict(lambda: [0, 0])
        for name, v, r, cls, *_ in rows:
            summ[cls][0] += 1
            summ[cls][1] += v["len"]
        for cls in sorted(summ, key=lambda c: -summ[c][1]):
            print(f"[{tgt}]   {cls:<18} {summ[cls][0]:>4} nodes  {summ[cls][1]:>9,} bp")
        keep = {tgt, f"{tgt}_repeat", "mtpt", "pt_shared_with_mt"}
        est = sum(v["len"] * copies for name, v, r, cls, copies, *_ in rows if cls in keep)
        print(f"[{tgt}] estimated genome size (target classes, repeats x copies): {est:,} bp")
        unres = [x for x in rows if x[5] == "exceeds_insert"]
        if unres:
            print(f"[{tgt}] {len(unres)} repeat node(s) longer than the median insert ({args.insert} bp) -- "
                  "alternative paths through them cannot be resolved with these reads; report the graph, "
                  "do not force a single circle.")
    print(f"Repeat candidate table written to: {out_tsv}  (Bandage colours: bandage_colours.csv)")


if __name__ == "__main__":
    main()
