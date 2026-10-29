#!/usr/bin/env python3
"""
trace_dead_ends.py -- Stage C3: follow every dead end of the curated organelle
graph through the UNFILTERED GetOrganelle graph (slimmed_assembly_graph.fastg,
written before GetOrganelle's depth selection) and report where it goes.

Why: in plant mitogenomes assembled from short reads, dead ends frequently sit
at the boundary of a plastid-derived insert (MTPT) whose sequence is identical
to the plastome. Reads of the insert cannot be told apart from plastid reads,
the insert's edges carry plastid-level depth, and GetOrganelle's depth filter
cuts them, leaving an apparent "end". Tracing through those discarded edges
(i) identifies what each dead end runs into, and (ii) finds candidate
connections between dead ends through other-organelle sequence -- hypotheses
to test by PCR or long reads, never joins to be made automatically.
(Fuchsia magellanica: all four mitochondrial dead ends were MTPT boundaries;
two of them were linked through ~3.1 kb of plastid inverted repeat.)

Inputs
  --gfa        curated/selected graph (GetOrganelle *.selected_graph.gfa)
  --fastg      slimmed_assembly_graph.fastg from the same GetOrganelle run
  --labels     slimmed_assembly_graph.csv (gene labels per edge), optional
  --other-ref  FASTA of the other organelle (plastome for an mt graph);
               optional: with blastn in PATH, bridge edges get coordinates
  --exclude    comma-separated node names to ignore as graph members
               (e.g. nodes Stage C classified as plastid)
  --depth      baseline depth D of the curated graph (k-mer); if omitted the
               length-weighted median depth of the curated nodes is used
  --max-bridge maximum bridge length in bp (default 40000)

Output: <outdir>/dead_ends.tsv, one row per (dead end, reachable node):
  dead_end, end, first_edge, first_edge_depth_xD, first_edge_loci,
  reaches, bridge_bp, bridge_edges, bridge_depth_xD, bridge_loci,
  bridge_other_coords
and a plain-language summary on stdout.
"""
import argparse
import collections
import csv
import heapq
import os
import re
import shutil
import subprocess

HDR = re.compile(r"EDGE_(\d+)_length_(\d+)_cov_([\d.]+)(')?")
FLIP = {"+": "-", "-": "+"}


def read_fastg(path):
    info, seq, out = {}, {}, collections.defaultdict(set)
    cur = None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip().rstrip(";")
                head, _, tail = h.partition(":")
                m = HDR.match(head)
                e, L, c, rv = m.groups()
                info[e] = (int(L), float(c))
                o = "-" if rv else "+"
                cur = e if o == "+" else None
                if cur:
                    seq[e] = []
                for t in tail.split(",") if tail else []:
                    mm = HDR.match(t)
                    out[(e, o)].add((mm.group(1), "-" if mm.group(4) else "+"))
            elif cur:
                seq[cur].append(line.strip())
    return info, {k: "".join(v) for k, v in seq.items()}, out


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))


def read_gfa(path):
    nodes, links, seqs = {}, [], {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "S":
                tags = {t.split(":", 2)[0]: t.split(":", 2)[2] for t in p[3:] if t.count(":") >= 2}
                L = int(tags.get("LN", len(p[2])))
                d = float(tags["DP"]) if "DP" in tags else (
                    float(tags["KC"]) / L if "KC" in tags else (float(tags["RC"]) / L if "RC" in tags else 0.0))
                nodes[p[1]] = (L, d)
                seqs[p[1]] = p[2]
            elif p[0] == "L":
                links.append((p[1], p[2], p[3], p[4]))
    return nodes, links, seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--fastg", required=True)
    ap.add_argument("--labels")
    ap.add_argument("--other-ref")
    ap.add_argument("--exclude", default="")
    ap.add_argument("--depth", type=float)
    ap.add_argument("--max-bridge", type=int, default=40000)
    ap.add_argument("--overlap", type=int, default=127, help="k-mer overlap between graph edges")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    nodes, links, nseq = read_gfa(a.gfa)
    excl = {x for x in a.exclude.split(",") if x}
    nodes = {n: v for n, v in nodes.items() if n not in excl}
    links = [l for l in links if l[0] in nodes and l[2] in nodes]
    info, seq, out = read_fastg(a.fastg)
    labels = {}
    if a.labels and os.path.exists(a.labels):
        with open(a.labels) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                labels[r["EDGE"]] = f"{r.get('database', '')}:{r.get('loci', '')}"

    if a.depth:
        D = a.depth
    else:
        pairs = sorted((d, L) for L, d in nodes.values())
        half, acc, D = sum(L for _, L in pairs) / 2, 0, pairs[-1][0]
        for d, L in pairs:
            acc += L
            if acc >= half:
                D = d
                break

    # curated node -> component SPAdes edges (GetOrganelle joins edge ids with "_")
    comp_of = {}
    for n in nodes:
        for e in n.split("_"):
            comp_of[e] = n
    member = set(comp_of)

    # dead ends of the curated graph
    deg = collections.Counter()
    for x, xo, y, yo in links:
        deg[(x, "R" if xo == "+" else "L")] += 1
        deg[(y, "L" if yo == "+" else "R")] += 1
    dead = [(n, e) for n in nodes for e in "LR" if deg[(n, e)] == 0]

    def exits(n, end):
        """oriented fastg edge leaving curated node n at the given end (L or R),
        located by matching the node's terminal sequence to a component edge,
        and its neighbours outside the curated graph"""
        s, k = nseq.get(n, ""), 200
        exit_oe = None
        for e in n.split("_"):
            if e not in seq or not s:
                continue
            fw, bw = seq[e], rc(seq[e])
            kk = min(k, len(fw), len(s))
            if end == "R":
                if s[-kk:] == fw[-kk:]:
                    exit_oe = (e, "+")
                elif s[-kk:] == bw[-kk:]:
                    exit_oe = (e, "-")
            else:
                if s[:kk] == fw[:kk]:
                    exit_oe = (e, "-")
                elif s[:kk] == bw[:kk]:
                    exit_oe = (e, "+")
            if exit_oe:
                break
        if not exit_oe:
            return []
        return [(exit_oe, nb) for nb in out[exit_oe] if nb[0] not in member]

    rows = []
    for n, end in dead:
        firsts = exits(n, end)
        if not firsts:
            rows.append(dict(dead_end=n, end=end, first_edge="", first_edge_depth_xD="", first_edge_loci="",
                             reaches="none (no edge beyond this end in the unfiltered graph)", bridge_bp="",
                             bridge_edges="", bridge_depth_xD="", bridge_loci="", bridge_other_coords=""))
            continue
        for (src, nb) in firsts:
            fe = nb[0]
            # shortest bridge (by bp) through non-member edges to a different curated node
            pq = [(info[fe][0] - a.overlap, nb, [nb])]
            best, hit = {}, None
            while pq:
                d, cur, path = heapq.heappop(pq)
                for nxt in out[cur]:
                    if nxt[0] in member:
                        if comp_of[nxt[0]] != n:
                            hit = (d, comp_of[nxt[0]], path)
                            break
                        continue
                    nd = d + info[nxt[0]][0] - a.overlap
                    if nd > a.max_bridge or best.get(nxt, 1e18) <= nd:
                        continue
                    best[nxt] = nd
                    heapq.heappush(pq, (nd, nxt, path + [nxt]))
                if hit:
                    break
            base = dict(dead_end=n, end=end, first_edge=fe + nb[1],
                        first_edge_depth_xD=f"{info[fe][1] / D:.2f}", first_edge_loci=labels.get(fe, ""))
            if hit:
                d, tgt, path = hit
                base.update(reaches=tgt, bridge_bp=d, bridge_edges=",".join(e + o for e, o in path),
                            bridge_depth_xD=",".join(f"{info[e][1] / D:.1f}" for e, _ in path),
                            bridge_loci=";".join(sorted({labels.get(e, "") for e, _ in path} - {""})),
                            bridge_other_coords="")
            else:
                base.update(reaches=f"no curated node within {a.max_bridge} bp", bridge_bp="", bridge_edges="",
                            bridge_depth_xD="", bridge_loci="", bridge_other_coords="")
            rows.append(base)

    # coordinates of bridge edges on the other organelle (optional, BLASTn)
    if a.other_ref and shutil.which("blastn") and any(r["bridge_edges"] for r in rows):
        fa = os.path.join(a.outdir, "bridge_edges.fasta")
        want = {e[:-1] for r in rows for e in r["bridge_edges"].split(",") if e}
        with open(fa, "w") as f:
            for e in want:
                f.write(f">{e}\n{seq.get(e, '')}\n")
        db = os.path.join(a.outdir, "other_ref_db")
        subprocess.run(["makeblastdb", "-in", a.other_ref, "-dbtype", "nucl", "-out", db],
                       check=True, stdout=subprocess.DEVNULL)
        tsv = os.path.join(a.outdir, "bridge_edges_vs_other.tsv")
        subprocess.run(["blastn", "-query", fa, "-db", db, "-outfmt", "6 qseqid sstart send pident length qlen",
                        "-evalue", "1e-20", "-perc_identity", "90", "-num_threads", str(a.threads),
                        "-out", tsv], check=True)
        coords = collections.defaultdict(list)
        with open(tsv) as f:
            for line in f:
                q, s1, s2, pid, ln, ql = line.split("\t")
                if int(ln) >= 0.8 * int(ql):
                    coords[q].append(f"{min(int(s1), int(s2))}-{max(int(s1), int(s2))}({float(pid):.1f}%)")
        for r in rows:
            if r["bridge_edges"]:
                r["bridge_other_coords"] = " | ".join(
                    f"{e[:-1]}:" + "/".join(coords.get(e[:-1], ["nd"])) for e in r["bridge_edges"].split(","))

    cols = ["dead_end", "end", "first_edge", "first_edge_depth_xD", "first_edge_loci", "reaches", "bridge_bp",
            "bridge_edges", "bridge_depth_xD", "bridge_loci", "bridge_other_coords"]
    with open(os.path.join(a.outdir, "dead_ends.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"[dead ends] {len(dead)} dead end(s) in the curated graph (baseline D={D:.1f})")
    for r in rows:
        if r["first_edge"]:
            print(f"[dead ends]   {r['dead_end']} ({r['end']}) -> edge {r['first_edge']} at "
                  f"{r['first_edge_depth_xD']} x D [{r['first_edge_loci'] or 'unlabelled'}]; "
                  + (f"candidate bridge to {r['reaches']} via {r['bridge_bp']} bp"
                     if r["bridge_bp"] != "" else r["reaches"]))
        else:
            print(f"[dead ends]   {r['dead_end']} ({r['end']}): {r['reaches']}")
    print("[dead ends] Bridges are hypotheses (other-organelle-identical sequence cannot be assigned with short "
          f"reads); validate by PCR/long reads. Table: {os.path.join(a.outdir, 'dead_ends.tsv')}")


if __name__ == "__main__":
    main()
