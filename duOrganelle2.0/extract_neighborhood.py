#!/usr/bin/env python3
"""
extract_neighborhood.py -- pull the neighbourhood of selected nodes out of a
large SPAdes assembly graph (GFA), so it can be inspected or shared without
moving the whole ~1 GB graph.

Typical use: after DuOrganelle's graph route, check whether the dead ends of
the mitochondrial graph connect, in the GLOBAL SPAdes graph, to nodes that
GetOrganelle's slimming step discarded (e.g. nodes at mitochondrial depth
without mitochondrial genes).

  python3 extract_neighborhood.py \
     --gfa out_duo_v2f/spades/assembly_graph_with_scaffolds.gfa \
     --nodes 7564413,7570939,7576910,7523053 --distance 3 \
     --out tips_neighbourhood

Node names are SPAdes edge numbers (EDGE_7564413 in the fastg = segment
"7564413" in the GFA). Composite GetOrganelle names such as
7523053_7486118_2987897_7593259 can be given as the full name; every number
in it is used as a start node.

Writes <out>.gfa (segments + links within the neighbourhood) and <out>.tsv
(node, length, depth KC/len or DP, distance from the nearest start node,
start node). Depth units match the GetOrganelle graph (k-mer coverage).
"""
import argparse
from collections import defaultdict, deque


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--nodes", required=True, help="comma-separated node names")
    ap.add_argument("--distance", type=int, default=3)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    starts = set()
    for n in a.nodes.split(","):
        starts.update(x for x in n.strip().split("_") if x)

    adj = defaultdict(set)
    with open(a.gfa) as f:                      # pass 1: links only
        for line in f:
            if line.startswith("L\t"):
                p = line.split("\t", 6)
                if p[1] != p[3]:
                    adj[p[1]].add(p[3])
                    adj[p[3]].add(p[1])

    dist, origin = {}, {}
    q = deque()
    for s in starts:
        dist[s], origin[s] = 0, s
        q.append(s)
    while q:
        x = q.popleft()
        if dist[x] >= a.distance:
            continue
        for y in adj.get(x, ()):
            if y not in dist:
                dist[y], origin[y] = dist[x] + 1, origin[x]
                q.append(y)
    keep = set(dist)

    rows = []
    with open(a.gfa) as f, open(a.out + ".gfa", "w") as og:   # pass 2
        og.write("H\tVN:Z:1.0\n")
        for line in f:
            if line.startswith("S\t"):
                p = line.rstrip("\n").split("\t")
                if p[1] not in keep:
                    continue
                og.write(line if line.endswith("\n") else line + "\n")
                tags = {t.split(":", 2)[0]: t.split(":", 2)[2] for t in p[3:] if t.count(":") >= 2}
                ln = int(tags.get("LN", len(p[2]) if p[2] != "*" else 0))
                if "DP" in tags:
                    dp = float(tags["DP"])
                elif "KC" in tags and ln:
                    dp = float(tags["KC"]) / ln
                elif "RC" in tags and ln:
                    dp = float(tags["RC"]) / ln
                else:
                    dp = float("nan")
                rows.append((p[1], ln, dp, dist[p[1]], origin[p[1]]))
            elif line.startswith("L\t"):
                p = line.split("\t", 6)
                if p[1] in keep and p[3] in keep:
                    og.write(line if line.endswith("\n") else line + "\n")

    with open(a.out + ".tsv", "w") as ot:
        ot.write("node\tlength\tdepth\tdistance\tstart_node\n")
        for r in sorted(rows, key=lambda r: (r[4], r[3], -r[1])):
            ot.write(f"{r[0]}\t{r[1]}\t{r[2]:.2f}\t{r[3]}\t{r[4]}\n")
    missing = starts - {r[0] for r in rows}
    print(f"{len(rows)} nodes within distance {a.distance} of {len(starts)} start node(s) -> {a.out}.gfa / {a.out}.tsv")
    if missing:
        print("WARNING: start node(s) not found in this GFA:", ",".join(sorted(missing)))


if __name__ == "__main__":
    main()
