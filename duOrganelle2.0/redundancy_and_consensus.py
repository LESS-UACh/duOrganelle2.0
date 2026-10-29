#!/usr/bin/env python3
"""
redundancy_and_consensus.py

Given a fresh set of NOVOPlasty contigs and the segment set accumulated from
previous rounds, this script:
  1. Runs all-vs-all BLASTn to detect redundant/overlapping contigs.
  2. Reorients redundant contigs to a common strand and builds a majority-rule
     consensus (minimap2 + pileup) for the largest connected cluster.
  3. Masks all content already present in the accumulated segment set and
     extracts genuinely novel, non-redundant sequence.

This formalizes the manual workflow used to characterize the repeat-mediated
multipartite mitochondrial genome of Fuchsia magellanica.
"""
import argparse
import os
import subprocess
import sys
from collections import defaultdict, Counter


def sh(cmd):
    """Run cmd, and if it fails, SHOW why (stdout+stderr), instead of swallowing
    the real error behind a generic CalledProcessError. Bug found in real use:
    the previous version had callers append '2>/dev/null' to each command
    string, which discarded genuinely useful errors (e.g. makeblastdb's "File
    ... does not exist" when a --anti-reference/-r_mt/-r_pt path was wrong) --
    all the user ever saw was a bare CalledProcessError with no explanation.
    Capturing output here (and NOT redirecting stderr away in the command
    strings themselves) means routine tool banners stay invisible on success
    but the real message is shown the moment something actually fails.

    2nd bug found in real use: text=True makes subprocess.run() decode
    stdout/stderr as UTF-8 the moment the process exits, BEFORE we ever get to
    check the return code -- so if a tool writes even one stray non-UTF-8 byte
    to stdout/stderr (observed: samtools sort emitting binary bytes on stdout
    on some builds even with -o given), subprocess.run() itself raises
    UnicodeDecodeError and crashes the whole pipeline before this function's
    own error-reporting logic ever runs -- hiding the real problem exactly
    like the 2>/dev/null bug above did. Fixed by capturing raw bytes and
    decoding leniently (errors='replace') only when building the message, so
    a stray binary byte degrades to a replacement character instead of
    crashing the run."""
    result = subprocess.run(cmd, shell=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"\nCommand failed (exit {result.returncode}): {cmd}\n"
            f"--- stdout ---\n{result.stdout.decode('utf-8', errors='replace')}"
            f"--- stderr ---\n{result.stderr.decode('utf-8', errors='replace')}"
        )
    return result


def load_fasta(path):
    seqs = {}
    name = None
    seq = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return seqs
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            seqs[name] = "".join(seq)
    return seqs


def write_fasta(seqs, path):
    with open(path, "w") as f:
        for k, v in seqs.items():
            f.write(f">{k}\n")
            for i in range(0, len(v), 70):
                f.write(v[i:i + 70] + "\n")


def revcomp(s):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(comp)[::-1]


def run_blast(query, db_fasta, out_tsv, identity, minlen, threads):
    sh(f"makeblastdb -in {db_fasta} -dbtype nucl -out {db_fasta}.db")
    sh(
        f"blastn -query {query} -db {db_fasta}.db "
        f'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" '
        f"-evalue 1e-20 -word_size 20 -num_threads {threads} -out {out_tsv}"
    )


def parse_hits(tsv, min_identity, min_len):
    hits = []
    if not os.path.exists(tsv):
        return hits
    with open(tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 13:
                continue
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen = parts[:13]
            if qseqid == sseqid:
                continue
            if float(pident) < min_identity or int(length) < min_len:
                continue
            hits.append((qseqid, sseqid, int(length), int(sstart), int(send)))
    return hits


def find_redundant_clusters(seqs, hits, coverage_frac=0.15):
    lens = {k: len(v) for k, v in seqs.items()}
    best = defaultdict(int)
    for q, s, length, sstart, send in hits:
        key = tuple(sorted([q, s]))
        best[key] = max(best[key], length)

    adj = defaultdict(set)
    for (a, b), length in best.items():
        frac = length / min(lens[a], lens[b])
        if frac > coverage_frac:
            adj[a].add(b)
            adj[b].add(a)

    seen = set()
    clusters = []
    for node in seqs:
        if node in seen:
            continue
        stack = [node]
        comp = set()
        while stack:
            n = stack.pop()
            if n in comp:
                continue
            comp.add(n)
            stack.extend(adj[n] - comp)
        seen |= comp
        clusters.append(comp)
    return clusters


def build_consensus(seqs, cluster, hits, outdir, threads):
    """Reorient cluster members to the longest member's strand and build a
    majority-rule consensus via minimap2 + pileup."""
    members = sorted(cluster, key=lambda k: -len(seqs[k]))
    ref_name = members[0]
    ref_seq = seqs[ref_name]

    orient = {ref_name: "+"}
    for m in members[1:]:
        best_hit = None
        for q, s, length, sstart, send in hits:
            if q == m and s == ref_name:
                if best_hit is None or length > best_hit[0]:
                    best_hit = (length, sstart, send)
        if best_hit:
            _, sstart, send = best_hit
            orient[m] = "+" if sstart < send else "-"
        else:
            orient[m] = "+"  # fallback, no direct hit to reference

    reoriented = {}
    for m in members:
        s = seqs[m]
        if orient.get(m) == "-":
            s = revcomp(s)
        reoriented[m] = s

    ref_path = os.path.join(outdir, "cluster_reference.fasta")
    reads_path = os.path.join(outdir, "cluster_reads.fasta")
    write_fasta({ref_name: reoriented[ref_name]}, ref_path)
    write_fasta({m: reoriented[m] for m in members[1:]}, reads_path)

    if len(members) == 1:
        return reoriented[ref_name]

    sam = os.path.join(outdir, "aln.sam")
    bam = os.path.join(outdir, "aln.sorted.bam")
    sh(f"minimap2 -ax asm5 --secondary=no -t {threads} {ref_path} {reads_path} > {sam}")
    sh(f"samtools sort -o {bam} {sam}")
    sh(f"samtools index {bam}")
    pileup = os.path.join(outdir, "pileup.txt")
    sh(f"samtools mpileup -f {ref_path} -a {bam} > {pileup}")

    def parse_pileup_bases(bases, ref_base):
        calls = []
        i, n = 0, len(bases)
        while i < n:
            c = bases[i]
            if c == "^":
                i += 2
                continue
            if c == "$":
                i += 1
                continue
            if c in "+-":
                i += 1
                num = ""
                while i < n and bases[i].isdigit():
                    num += bases[i]
                    i += 1
                i += int(num)
                continue
            if c in ".,":
                calls.append(ref_base.upper())
                i += 1
                continue
            if c in "ACGTNacgtn":
                calls.append(c.upper())
                i += 1
                continue
            if c == "*":
                i += 1
                continue
            i += 1
        return calls

    consensus = list(reoriented[ref_name])
    with open(pileup) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            _, pos, ref_base, depth, bases, _ = parts[:6]
            pos = int(pos)
            calls = parse_pileup_bases(bases, ref_base)
            votes = calls + [ref_base.upper()]
            best_base, _ = Counter(votes).most_common(1)[0]
            if best_base in "ACGT":
                consensus[pos - 1] = best_base
    return "".join(consensus)


def merge_intervals(ivs):
    ivs = sorted(ivs)
    merged = []
    for lo, hi in ivs:
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def screen_against_other_organelle(seqs, anti_reference_fasta, outdir, identity, threads, frac_threshold):
    """Flag/remove any sequence whose length is mostly explained by near-identical
    matches to anti_reference_fasta -- the CURRENT best assembly of the OTHER
    organelle track (e.g. the plastid accumulation, when screening the mitochondrial
    track, or vice versa) -- BEFORE it is treated as novel content for THIS track.
    This is the DuOrganelle generalization of repeatMito's original chloroplast-only
    screen: same mechanism, but symmetric and applied bidirectionally between two
    organelle tracks that are being assembled together rather than against one fixed,
    pre-existing reference genome.

    Returns (clean_seqs, flagged_seqs, coverage_report: dict[name] -> pct_other_organelle_covered).
    """
    if not anti_reference_fasta or not seqs:
        return seqs, {}, {}

    tmp_fasta = os.path.join(outdir, "_organelle_screen_query.fasta")
    write_fasta(seqs, tmp_fasta)
    tsv = os.path.join(outdir, "vs_other_organelle.tsv")
    run_blast(tmp_fasta, anti_reference_fasta, tsv, identity, 30, threads)

    covered_q = defaultdict(list)
    with open(tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 13:
                continue
            q, s, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen = parts[:13]
            if float(pident) < identity:
                continue
            lo, hi = sorted([int(qstart), int(qend)])
            covered_q[q].append((lo, hi))

    clean, flagged, report = {}, {}, {}
    for name, seq in seqs.items():
        merged = merge_intervals(covered_q.get(name, []))
        covered_bp = sum(hi - lo + 1 for lo, hi in merged)
        pct = 100.0 * covered_bp / len(seq) if len(seq) else 0.0
        report[name] = pct
        if pct >= frac_threshold * 100:
            flagged[name] = seq
        else:
            clean[name] = seq
    return clean, flagged, report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--accumulated", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--identity", type=float, default=90)
    ap.add_argument("--minlen", type=int, default=200)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--anti-reference", "--chloroplast", dest="anti_reference",
                     help="FASTA of the CURRENT best assembly of the OTHER organelle track (e.g. the "
                          "plastid accumulation so far, when this invocation is processing the "
                          "mitochondrial track, or vice versa -- updated dynamically round to round by "
                          "the DuOrganelle orchestrator, NOT a fixed pre-existing reference). When given, "
                          "every per-cluster consensus is screened against it BEFORE being treated as "
                          "novel content for THIS track; any consensus with >= --other-frac-threshold of "
                          "its length explained by near-identical matches to the other organelle is "
                          "written to flagged_other_organelle.fasta instead of novel_segment.fasta, and "
                          "is excluded from this track's accumulated segment set / next seed. "
                          "'--chloroplast' is kept as a backward-compatible alias for classic single-"
                          "organelle repeatMito usage (a fixed, pre-assembled chloroplast genome).")
    ap.add_argument("--other-frac-threshold", "--plastid-frac-threshold", dest="other_frac_threshold",
                     type=float, default=0.5,
                     help="Fraction (0-1) of a consensus sequence's length that must be covered by "
                          "near-identical BLASTn hits to --anti-reference to flag it as cross-organelle "
                          "contamination (default 0.5). Genuine MTPT-bearing segments normally show small "
                          "fractional coverage (a few percent to ~20%%), not a majority.")
    ap.add_argument("--raw-other-frac-threshold", "--raw-plastid-frac-threshold", dest="raw_other_frac_threshold",
                     type=float, default=0.3,
                     help="Same idea as --other-frac-threshold, but applied to each RAW NOVOPlasty contig "
                          "BEFORE all-vs-all clustering/consensus-building (default 0.3, stricter than the "
                          "post-consensus check). This matters because single-linkage clustering can "
                          "transitively merge a mostly-other-organelle contig into the same cluster as "
                          "genuinely-this-track contigs if they share even a short region of homology "
                          "(e.g. a real MTPT) -- the resulting consensus would then be a chimera that the "
                          "post-consensus screen alone might not flag if the contaminating fraction of the "
                          "FINAL merged consensus dilutes below --other-frac-threshold.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    seqs = load_fasta(args.contigs)
    if not seqs:
        write_fasta({}, os.path.join(args.outdir, "novel_segment.fasta"))
        return

    # 0) cross-organelle pre-screen on RAW contigs, before they can be clustered/merged
    # with other contigs via shared homology (e.g. a real MTPT bridging a mostly-other-organelle
    # contig to genuinely-this-track ones).
    if args.anti_reference:
        seqs, flagged_raw, raw_report = screen_against_other_organelle(
            seqs, args.anti_reference, args.outdir, args.identity, args.threads,
            args.raw_other_frac_threshold,
        )
        raw_report_path = os.path.join(args.outdir, "other_organelle_screen_report_RAW_contigs.tsv")
        with open(raw_report_path, "w") as rf:
            rf.write("raw_contig\tpct_length_matching_other_organelle\tverdict\n")
            for name, pct in sorted(raw_report.items(), key=lambda x: -x[1]):
                verdict = "FLAGGED_AS_OTHER_ORGANELLE_excluded_before_clustering" if name in flagged_raw else "ok"
                rf.write(f"{name}\t{pct:.1f}\t{verdict}\n")
        if flagged_raw:
            write_fasta(flagged_raw, os.path.join(args.outdir, "flagged_other_organelle_raw_contigs.fasta"))
            print(f"WARNING: {len(flagged_raw)} RAW contig(s) excluded BEFORE clustering "
                  f"(>= {args.raw_other_frac_threshold*100:.0f}% matching the other organelle's current "
                  f"assembly): {', '.join(flagged_raw)}. This prevents them from dragging the other "
                  f"organelle's content into a cluster/consensus with genuinely-this-track contigs via "
                  f"shared homology (e.g. a real MTPT). See {raw_report_path}.")
        if not seqs:
            write_fasta({}, os.path.join(args.outdir, "novel_segment.fasta"))
            print("All raw contigs this round were flagged as cross-organelle contamination; nothing left to cluster.")
            return

    # 1) all-vs-all redundancy within this round's contigs (post cross-organelle pre-screen)
    clean_contigs_fasta = os.path.join(args.outdir, "_raw_contigs_post_cross_organelle_screen.fasta")
    write_fasta(seqs, clean_contigs_fasta)
    aa_tsv = os.path.join(args.outdir, "all_vs_all.tsv")
    run_blast(clean_contigs_fasta, clean_contigs_fasta, aa_tsv, args.identity, args.minlen, args.threads)
    hits = parse_hits(aa_tsv, args.identity, args.minlen)
    clusters = find_redundant_clusters(seqs, hits)

    consensus_seqs = {}
    for i, cluster in enumerate(clusters):
        cdir = os.path.join(args.outdir, f"cluster_{i}")
        os.makedirs(cdir, exist_ok=True)
        if len(cluster) == 1:
            name = next(iter(cluster))
            consensus_seqs[f"cluster{i}_{name}"] = seqs[name]
        else:
            cons = build_consensus(seqs, cluster, hits, cdir, args.threads)
            consensus_seqs[f"cluster{i}_consensus_n{len(cluster)}"] = cons

    consolidated_round = os.path.join(args.outdir, "round_consolidated.fasta")
    write_fasta(consensus_seqs, consolidated_round)

    # 1b) cross-organelle safeguard: screen every cluster consensus against the other
    # organelle's current assembly BEFORE it can be folded into this track's accumulated set.
    if args.anti_reference:
        consensus_seqs, flagged_seqs, other_report = screen_against_other_organelle(
            consensus_seqs, args.anti_reference, args.outdir, args.identity, args.threads,
            args.other_frac_threshold,
        )
        report_path = os.path.join(args.outdir, "other_organelle_screen_report.tsv")
        with open(report_path, "w") as rf:
            rf.write("sequence\tpct_length_matching_other_organelle\tverdict\n")
            for name, pct in sorted(other_report.items(), key=lambda x: -x[1]):
                verdict = "FLAGGED_AS_OTHER_ORGANELLE" if name in flagged_seqs else "ok"
                rf.write(f"{name}\t{pct:.1f}\t{verdict}\n")
        if flagged_seqs:
            write_fasta(flagged_seqs, os.path.join(args.outdir, "flagged_other_organelle_contamination.fasta"))
            print(f"WARNING: {len(flagged_seqs)} cluster consensus sequence(s) flagged as cross-organelle "
                  f"contamination (>= {args.other_frac_threshold*100:.0f}% of length matches the other "
                  f"organelle's current assembly at >= {args.identity}% identity) and EXCLUDED from this "
                  f"track's segment set. See {report_path} and flagged_other_organelle_contamination.fasta.")
            for name, pct in other_report.items():
                if name in flagged_seqs:
                    print(f"    - {name}: {pct:.1f}% other-organelle-matching -> excluded")
    else:
        print("NOTE: no --anti-reference given to redundancy_and_consensus.py; skipping the cross-organelle "
              "safeguard for this round (expected in round 1 of a track before the other track has "
              "produced anything to screen against; should not remain true from round 2 onward).", file=sys.stderr)

    # 2) mask content already present in the accumulated segment set -> novel segment(s)
    if os.path.exists(args.accumulated) and os.path.getsize(args.accumulated) > 0:
        vs_acc_tsv = os.path.join(args.outdir, "vs_accumulated.tsv")
        run_blast(consolidated_round, args.accumulated, vs_acc_tsv, 90, 100, args.threads)
        known_hits = parse_hits(vs_acc_tsv, 90, 100)
        covered = defaultdict(list)
        for q, s, length, sstart, send in known_hits:
            pass
        # need query-side coordinates, not subject; rerun parse with qstart/qend
        covered_q = defaultdict(list)
        with open(vs_acc_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 13:
                    continue
                q, s, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen = parts[:13]
                if float(pident) < 90 or int(length) < 100:
                    continue
                lo, hi = sorted([int(qstart), int(qend)])
                covered_q[q].append((lo, hi))

        def merge(ivs):
            ivs = sorted(ivs)
            merged = []
            for lo, hi in ivs:
                if merged and lo <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
                else:
                    merged.append((lo, hi))
            return merged

        novel_parts = {}
        for name, seq in consensus_seqs.items():
            L = len(seq)
            merged = merge(covered_q.get(name, []))
            pos = 1
            gaps = []
            for lo, hi in merged:
                if lo > pos:
                    gaps.append((pos, lo - 1))
                pos = max(pos, hi + 1)
            if pos <= L:
                gaps.append((pos, L))
            for gi, (lo, hi) in enumerate(gaps):
                if hi - lo + 1 < 200:
                    continue
                novel_parts[f"{name}_novel_{lo}_{hi}"] = seq[lo - 1:hi]
    else:
        novel_parts = consensus_seqs

    write_fasta(novel_parts, os.path.join(args.outdir, "novel_segment.fasta"))
    print(f"Round complete: {len(clusters)} redundancy cluster(s), "
          f"{sum(len(v) for v in novel_parts.values())} bp novel sequence extracted.")


if __name__ == "__main__":
    main()
