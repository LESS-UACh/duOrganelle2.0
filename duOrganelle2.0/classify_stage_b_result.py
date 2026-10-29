#!/usr/bin/env python3
"""
classify_stage_b_result.py -- Stage B (GetOrganelle) plain-language status flag.

MOTIVATION (added after the Oenothera villaricae validation run): GetOrganelle's
own get_org.log.txt already says exactly how much to trust its Stage B result --
a clean, closed result never prints "confirm the incomplete result in Bandage",
while a fragmented/ambiguous one always does, alongside specific warnings
(self-loop nodes, conflicting connections, "Multiple isolated components",
"contamination or polymorphism"). Before this script, that signal was buried in
round_logs/getorganelle.log, which nobody reads line-by-line; report/run.log's
own last status line was just "Both tracks converged after round N", giving no
indication that, e.g., Stage B produced 15 disconnected mt scaffolds or that a
pt result was flagged for possible contamination. This script parses that
signal once per organelle and prints one unambiguous line, meant to be
tee'd into report/run.log immediately after Stage B, so the run's headline
summary can no longer look more finished than the assembly actually is.

REVISED (after the first re-run on Fuchsia magellanica with this script wired
in -- see duorganelle_jjnunez_final): the initial version only recognized
"Result status of X: N scaffold(s)". A genuinely successful GetOrganelle run
prints "Result status of X: circular genome" instead (no digit), which fell
through to the catch-all UNKNOWN branch -- a clean, fully circularized plastid
result was reported as unclassifiable. That run also surfaced a second,
legitimate GetOrganelle message this script didn't know about: "More than one
circular genome structure produced" / "different in SSC direction (two
flip-flop configurations)", which for a plastid genome with an inverted repeat
is normally a benign orientation ambiguity, not a sign of a bad assembly. Both
are handled explicitly now rather than defaulting to UNKNOWN or being silently
absorbed into the failure branch.

This is a classifier over GetOrganelle's own stated confidence, NOT an
independent quality assessment -- it does not decide whether a fragmented
result reflects genuine multipartite biology (see Fischer et al. 2022 for
Oenothera) versus a real data limitation (e.g. short-insert-only reads unable
to span a repeat) versus true contamination. It only makes sure that
distinction is surfaced instead of silently defaulting to "converged" =
"done".
"""
import argparse
import re
import sys

# Patterns that force an INCOMPLETE / FRAGMENTED verdict.
WARNING_PATTERNS = [
    ("Multiple isolated", "multiple isolated graph components (GetOrganelle's own 'Broken or contamination?' flag)"),
    ("contamination or polymorphism", "contamination/polymorphism warning"),
    ("insufficient graph", "graph disentangling fell back to insufficient-graph/linear mode (could not resolve a single path)"),
    ("Self-loop contig detected", "self-loop node(s) detected (ambiguous repeat junction)"),
    ("conflicting connections existed", "conflicting graph connections at a junction (ambiguous overlap)"),
]

# Patterns that are worth surfacing alongside EITHER verdict but do not by
# themselves mean the result is bad -- both are routine GetOrganelle output
# for repeat-rich or inverted-repeat-bearing genomes.
NOTE_PATTERNS = [
    ("More than one circular genome structure produced",
     "GetOrganelle produced more than one circular path variant (its own note: \"different in SSC direction "
     "(two flip-flop configurations)\") -- usually a benign plastid SSC-orientation ambiguity, but visually "
     "confirm in Bandage that the variants represent the same genome before picking one"),
    ("Multiple repeat patterns appeared in your data",
     "GetOrganelle flagged multiple repeat patterns in the data (several alternative PATH files written) -- "
     "inspect which path is best supported before treating any single one as final"),
]

INCOMPLETE_MARKER = "confirm the incomplete result in Bandage"


def classify(log_text: str, organelle: str) -> str:
    start_marker = f"Extracting {organelle} from the assemblies ..."
    end_marker = f"Extracting {organelle} from the assemblies finished."
    start = log_text.find(start_marker)
    end = log_text.find(end_marker)
    if start == -1:
        # get_organelle_from_assembly.py (graph route) writes a single-organelle
        # log without the "Extracting X from the assemblies" markers; in that
        # case the whole log belongs to this organelle.
        if f"Result status of {organelle}" in log_text or INCOMPLETE_MARKER in log_text:
            start, end = 0, -1
        else:
            return f"[{organelle}] Stage B verdict: NO BLOCK FOUND in log (organelle not requested, or GetOrganelle failed before reaching it) -- inspect the GetOrganelle log manually."
    block = log_text[start: end + len(end_marker) if end != -1 else len(log_text)]

    circular = bool(re.search(rf"Result status of {re.escape(organelle)}:\s*circular genome", block))
    m = re.search(rf"Result status of {re.escape(organelle)}:\s*(\d+)\s*scaffold", block)
    n_scaffolds = int(m.group(1)) if m else None

    incomplete = INCOMPLETE_MARKER in block
    hits = [label for pattern, label in WARNING_PATTERNS if pattern in block]
    notes = [label for pattern, label in NOTE_PATTERNS if pattern in block]
    note_str = (" NOTE: " + "; ".join(notes)) if notes else ""

    if circular and not incomplete:
        verdict = (
            "CLEAN -- GetOrganelle resolved a single circular genome. Still confirm visually in Bandage "
            "before deposition." + note_str
        )
    elif n_scaffolds == 1 and not incomplete and not hits:
        verdict = (
            "CLEAN, single scaffold, no repeat/contamination warnings -- consistent with a resolved "
            "circular/linear genome (still confirm visually in Bandage before deposition)." + note_str
        )
    elif hits or incomplete or (n_scaffolds is not None and n_scaffolds > 1):
        reason = "; ".join(hits) if hits else "GetOrganelle flagged the result incomplete"
        verdict = (
            f"INCOMPLETE / FRAGMENTED -- {reason}. "
            "Do NOT treat this as a closed genome. Could reflect genuine multipartite/repeat-mediated "
            "structure (see Fischer et al. 2022 for a documented case in Oenothera) or a data limitation "
            "(e.g. short-insert-only reads cannot span a repeat longer than the insert size) rather than "
            "a pipeline error -- inspect Stage C's repeat_candidates.tsv (class column) and the Bandage graph before "
            "drawing conclusions either way." + note_str
        )
    else:
        verdict = (
            "UNKNOWN -- no scaffold count or circular-genome status parsed and no recognized warning "
            "pattern matched; inspect getorganelle.log manually." + note_str
        )

    if circular:
        scaffold_str = "circular genome (multiple path variants)" if notes and "flip-flop" in note_str else "circular genome"
    elif n_scaffolds is not None:
        scaffold_str = f"{n_scaffolds} scaffold(s)"
    else:
        scaffold_str = "scaffold count not found"
    return f"[{organelle}] Stage B result: {scaffold_str}. Verdict: {verdict}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--log", required=True, help="path to round_logs/getorganelle.log (GetOrganelle's get_org.log.txt)")
    ap.add_argument("--organelle", required=True, choices=["embplant_mt", "embplant_pt"])
    args = ap.parse_args()

    with open(args.log, "r", errors="replace") as f:
        log_text = f.read()

    print(classify(log_text, args.organelle))


if __name__ == "__main__":
    main()
