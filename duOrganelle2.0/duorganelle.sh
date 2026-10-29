#!/bin/bash
# =============================================================================
# DuOrganelle
#
# Jointly, iteratively assembles BOTH the plastid and mitochondrial genomes of
# a plant from the SAME Illumina short-read dataset, with each organelle track
# used as the OTHER track's dynamically-updated cross-contamination screen --
# generalizing repeatMito's chloroplast-only BLASTn safeguard into a symmetric,
# bidirectional one that does not require either organelle genome to already
# be assembled beforehand.
#
# STATUS: run end-to-end on three datasets so far -- Fuchsia magellanica
# (Onagraceae), Oenothera villaricae (Onagraceae; validation run against the
# independently published, graph-resolved mitogenome of Fischer et al. 2022),
# and Punica granatum (Lythraceae; validation run against the independently
# published 382,774 bp / 7-circular-contig mitogenome, Sha et al. 2023 /
# GenBank OQ973289-OQ973295). Three things noted so far, now handled
# explicitly rather than silently:
#   (1) GetOrganelle's dual -F embplant_mt,embplant_pt mode DOES prefix
#       outputs per organelle type as assumed (embplant_mt.K85...gfa /
#       embplant_pt.K85...gfa) -- confirmed on a real run.
#   (2) Stage B can return successfully (exit 0, .gfa written) while
#       GetOrganelle's own log still flags the result as fragmented,
#       non-circularized, or possibly contaminated/polymorphic -- this was
#       previously invisible in report/run.log (which only recorded "Both
#       tracks converged"), and is now surfaced explicitly right after Stage
#       B by classify_stage_b_result.py (see STAGE B block below). A
#       fragmented Stage B result is not necessarily wrong: for a genuinely
#       multipartite/repeat-recombining mitogenome (as documented for
#       Oenothera by Fischer et al. 2022), several disconnected scaffolds can
#       be the CORRECT graph-level output. Whether a given fragmented result
#       reflects real biology or a data limitation (e.g. short-insert-only
#       reads unable to span a repeat -- no mate-pair/long-read data) is not
#       something this pipeline can determine automatically; the new status
#       line is meant to make sure that question gets asked instead of
#       silently defaulting to "converged" = "done".
#   (3) The mt track's Stage A "Genome Range" ceiling (12000-2000000, a wide
#       generic default inherited from repeatMito) is not a safety limit --
#       it is the size NOVOPlasty will keep extending toward if the read
#       graph lets it. On the Punica granatum run, mt round 1 produced a
#       contig of exactly 2,000,030 bp (the ceiling plus a small overhang):
#       NOVOPlasty hit the configured maximum rather than converging
#       naturally, against a species whose published mitogenome is 382,774
#       bp -- i.e. more than 5x the expected size, almost certainly runaway
#       extension into repetitive or nuclear-integrated (NUMT-like) content
#       rather than real single-copy mitochondrial sequence. Stage B's own
#       GetOrganelle read classification partly recovered from this (342,357
#       bp across 10 scaffolds, within ~10% of the true 382,774 bp), but the
#       Stage A contamination this produces is avoidable. New flags -g_mt /
#       -g_pt let a tighter, species-appropriate Genome Range be supplied
#       when the approximate true genome size is already known (e.g.
#       -g_mt 200000-450000 for a Punica-sized mitogenome) instead of always
#       falling back to the wide default; see make_novoplasty_config.py's
#       --genome-range, which already supported this; it just was not wired
#       through from here before.
#   (4) Organelle-ENRICHED libraries (Oct 2026, Fuchsia enriched run): with
#       most reads organellar, GetOrganelle's seed-and-extend recruitment
#       selects almost nothing (word dictionary 2 M -> 210 M in 3 rounds)
#       while the plastome, at ~100x mito depth, floods the mt graph. Stage
#       A2 (coverage_profile.py) now measures pt/mt median depth, organellar
#       read fraction and insert size, and picks the Stage B route:
#       RECRUIT (total DNA; unchanged behaviour) or GRAPH (partial plastid
#       depletion to ~2.5x mito depth, keeping MTPT-spanning pairs ->
#       SPAdes -> get_organelle_from_assembly.py). Stage C
#       (graph_depth_screen.py) now also labels graph nodes as mt /
#       repeat (n copies) / MTPT / plastid / low-depth, flags repeats longer
#       than the insert as unresolvable with short reads, and writes Bandage
#       colours. The k-mer list is now also cleaned automatically (odd,
#       <= 127, < read length).
#
# Usage (paired-end):
#   ./duorganelle.sh -m paired -1 R1.fq.gz -2 R2.fq.gz \
#       -s_mt seed_mt.fasta -s_pt seed_pt.fasta \
#       [-r_mt known_mito.fasta] [-r_pt known_plastid.fasta] \
#       -o outdir [-t threads] [-r max_rounds]
#
#   -s_mt / -s_pt : seed genes to start each track (e.g. matR/cox1/nad1 for mt;
#                   rbcL/matK/psbA for pt). Both required.
#   -r_mt / -r_pt : OPTIONAL warm-start reference for either organelle, if you
#                   already have one (e.g. a plastome from a related species,
#                   or your own prior assembly). Removes the round-1 blind
#                   spot described below for that track. Neither is required --
#                   DuOrganelle's point is to work without them -- but supply
#                   whichever you have.
#   -g_mt / -g_pt : OPTIONAL NOVOPlasty "Genome Range" override (format
#                   MIN-MAX) for either track, if you already know roughly
#                   how big the true genome is (e.g. from a related species'
#                   published assembly). Strongly recommended whenever
#                   available: the built-in defaults (mt: 12000-2000000, pt:
#                   120000-200000) are wide generic fallbacks, and NOVOPlasty
#                   will extend toward whichever bound it hits rather than
#                   treating it as a mere safety cap -- see header comment
#                   point (3) for a real case where this produced a >5x
#                   oversized mt contig.
#
# THE ROUND-1 BOOTSTRAP GAP (read this before running blind on precious data):
#   In round 1, if NEITHER -r_mt nor -r_pt is given, neither track has
#   anything yet to screen the other against -- both run with zero
#   cross-organelle protection for exactly one round, identical in principle
#   to the original repeatMito failure mode. From round 2 onward, each track
#   screens against the OTHER track's round-1+ accumulation, which is when
#   protection actually starts. This is a real, structural limitation of fully
#   blind bootstrapping, not an oversight -- it is presented as such in
#   README.md ("Limitations") rather than hidden. Passing
#   either -r_mt or -r_pt (even an imperfect/distantly-related one) removes
#   the gap for whichever direction you can supply it, and is recommended
#   whenever available.
# =============================================================================
set -euo pipefail

DUORGANELLE_VERSION="2.0.1"
CMDLINE="$0 $*"
MODE="paired"
THREADS=8
MAX_ROUNDS=5
MIN_NOVEL_BP=500
BLAST_IDENTITY=90
BLAST_MINLEN=200
OTHER_FRAC_THRESHOLD=0.5
RAW_OTHER_FRAC_THRESHOLD=0.3
REF_MT=""
REF_PT=""
GENOME_RANGE_MT=""
GENOME_RANGE_PT=""
# K-mer list passed to Stage B's GetOrganelle call (-k). Explicit rather than
# relying on GetOrganelle's own built-in default, so the range used here is
# guaranteed regardless of GetOrganelle version/organelle-type defaults.
# 127 is included because it is the k-mer that has previously produced the
# most-resolved single-connected-graph mitochondrial assembly for this
# project (see manuscript Results 3.2); 135 is included speculatively for
# read sets long enough to support it. Every value must be ODD and STRICTLY
# LESS than your read length (SPAdes requirement) -- check read length first
# (e.g. `zcat R1.fq.gz | head -400 | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c`)
# and drop any values >= read length via -k before running on new data.
KMER_LIST="21,45,65,85,105,115,127,135"
# NOTE (Oct 2026): default SPAdes builds reject k > 127 (MAX_K), so 135 is
# dropped automatically below, together with any even value or any value
# >= the detected read length. The cleaned list is logged.

# Stage A2 / Stage B route options (see coverage_profile.py):
ROUTE="auto"            # auto | recruit | graph
ENRICH_THRESHOLD=0.5    # organellar read fraction at/above which the library counts as enriched
PT_TARGET_RATIO=2.5     # plastid:mito depth kept after partial plastid depletion (graph route)
SPADES_MEM=120          # GB, SPAdes -m (graph route)
PROFILE_TIMEOUT=3h      # hard limit on the Stage A2 BBMap pass

usage() {
  echo "Usage: $0 -m paired -1 R1.fq.gz -2 R2.fq.gz -s_mt seed_mt.fasta -s_pt seed_pt.fasta -o outdir"
  echo "          [-r_mt known_mito.fasta] [-r_pt known_plastid.fasta] [-t threads] [-r max_rounds]"
  echo "          [-p other_frac_threshold] [-q raw_other_frac_threshold] [-k kmer_list]"
  echo "  -k kmer_list : comma-separated k-mer values for Stage B's GetOrganelle run"
  echo "                 (default: $KMER_LIST). Each value must be odd and < read length."
  echo "  -g_mt / -g_pt MIN-MAX : override NOVOPlasty's Genome Range for that track (see header"
  echo "                 comment point 3) -- strongly recommended when the approximate true genome"
  echo "                 size is already known."
  echo "  --route auto|recruit|graph : Stage B route (default auto: decided from the organellar read"
  echo "                 fraction measured in Stage A2). recruit = GetOrganelle from reads (total DNA);"
  echo "                 graph = partial plastid depletion + SPAdes + GetOrganelle from assembly (enriched)."
  echo "  --enrich-threshold F : organellar fraction for auto -> graph (default $ENRICH_THRESHOLD)"
  echo "  --pt-ratio F   : plastid:mito depth to keep after depletion (default $PT_TARGET_RATIO)"
  echo "  --spades-mem G : SPAdes memory limit in GB, graph route (default $SPADES_MEM)"
  echo "  See the header comment for the round-1 bootstrap caveat before running on irreplaceable data."
  exit 1
}

# Using long-form flags via getopts is awkward in bash; parse manually for the
# two-organelle-specific ones and fall back to getopts for the rest.
ARGS=()
while [ $# -gt 0 ]; do
  case "$1" in
    -v|--version) echo "DuOrganelle $DUORGANELLE_VERSION"; exit 0 ;;
    -s_mt) SEED_MT="$2"; shift 2 ;;
    -s_pt) SEED_PT="$2"; shift 2 ;;
    -r_mt) REF_MT="$2"; shift 2 ;;
    -r_pt) REF_PT="$2"; shift 2 ;;
    -g_mt) GENOME_RANGE_MT="$2"; shift 2 ;;
    -g_pt) GENOME_RANGE_PT="$2"; shift 2 ;;
    --route) ROUTE="$2"; shift 2 ;;
    --enrich-threshold) ENRICH_THRESHOLD="$2"; shift 2 ;;
    --pt-ratio) PT_TARGET_RATIO="$2"; shift 2 ;;
    --spades-mem) SPADES_MEM="$2"; shift 2 ;;
    *) ARGS+=("$1"); shift ;;
  esac
done
set -- "${ARGS[@]}"

while getopts "m:1:2:U:o:t:r:p:q:k:h" opt; do
  case $opt in
    m) MODE="$OPTARG" ;;
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    U) SEREADS="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    r) MAX_ROUNDS="$OPTARG" ;;
    p) OTHER_FRAC_THRESHOLD="$OPTARG" ;;
    q) RAW_OTHER_FRAC_THRESHOLD="$OPTARG" ;;
    k) KMER_LIST="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[ -z "${SEED_MT:-}" ] && { echo "ERROR: -s_mt (mitochondrial seed) required"; usage; }
[ -z "${SEED_PT:-}" ] && { echo "ERROR: -s_pt (plastid seed) required"; usage; }
[ -z "${OUTDIR:-}" ] && usage
if [ "$MODE" = "paired" ]; then
  { [ -z "${R1:-}" ] || [ -z "${R2:-}" ]; } && { echo "ERROR: paired mode requires -1 and -2"; usage; }
else
  [ -z "${SEREADS:-}" ] && { echo "ERROR: single mode requires -U"; usage; }
fi

# Check every file argument up front, so a missing/misspelled path is caught
# in under a second, before anything expensive starts (this was found the
# hard way in real use -- see CHANGELOG.md).
missing=0
check_file() {  # usage: check_file <flag-name> <path>
  if [ -n "$2" ] && [ ! -s "$2" ]; then
    echo "ERROR: $1 file not found or empty: $2" >&2
    missing=1
  fi
}
check_file "-s_mt" "$SEED_MT"
check_file "-s_pt" "$SEED_PT"
check_file "-r_mt" "$REF_MT"
check_file "-r_pt" "$REF_PT"
if [ "$MODE" = "paired" ]; then
  check_file "-1" "$R1"
  check_file "-2" "$R2"
else
  check_file "-U" "$SEREADS"
fi
if [ "$missing" -eq 1 ]; then
  echo "One or more input files above were not found relative to the current directory ($(pwd))." >&2
  echo "Double-check the paths (or your cwd) before rerunning -- nothing has been started yet." >&2
  exit 1
fi

case "$ROUTE" in auto|recruit|graph) ;; *) echo "ERROR: --route must be auto, recruit or graph"; usage ;; esac

# Check external tools up front (found the hard way: a run with BBTools
# missing silently skipped profiling and depletion and went on to assemble
# all reads). Profiling/depletion tools are mandatory unless --route recruit.
NEED_TOOLS=(NOVOPlasty4.3.5.pl get_organelle_from_reads.py blastn makeblastdb minimap2 samtools)
if [ "$ROUTE" != "recruit" ]; then
  NEED_TOOLS+=(bbduk.sh reformat.sh spades.py get_organelle_from_assembly.py)
fi
MISSING_TOOLS=()
for t in "${NEED_TOOLS[@]}"; do command -v "$t" >/dev/null 2>&1 || MISSING_TOOLS+=("$t"); done
if [ "${#MISSING_TOOLS[@]}" -gt 0 ]; then
  echo "ERROR: required tool(s) not in PATH: ${MISSING_TOOLS[*]}" >&2
  echo "       (BBTools: conda install -c bioconda bbmap). Nothing has been started yet." >&2
  exit 1
fi

mkdir -p "$OUTDIR"/{round_logs,novoplasty_runs/mt,novoplasty_runs/pt,consensus/mt,consensus/pt,bandage,report,profile}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Read length (max of the first 1000 reads) and a cleaned k-mer list: odd,
# <= 127 (SPAdes MAX_K) and < read length. The (zcat || true) subshell keeps
# head's early close from tripping `set -o pipefail`.
FIRST_READS="${R1:-${SEREADS:-}}"
READ_LEN=$( (zcat -f "$FIRST_READS" 2>/dev/null || true) | head -n 4000 | awk 'NR%4==2{print length($0)}' | sort -n | tail -1 )
READ_LEN=${READ_LEN:-0}
KMER_CLEAN=$(echo "$KMER_LIST" | tr ',' '\n' | awk -v rl="$READ_LEN" '$1%2==1 && $1<=127 && (rl==0 || $1<rl)' | paste -sd, -)
if [ -z "$KMER_CLEAN" ]; then KMER_CLEAN="21,33,55"; fi
if [ "$KMER_CLEAN" != "$KMER_LIST" ]; then
  echo "NOTE: k-mer list adjusted for read length ${READ_LEN} bp and SPAdes MAX_K 127: $KMER_LIST -> $KMER_CLEAN" >&2
fi
KMER_LIST="$KMER_CLEAN"

if [ -z "$REF_MT" ] && [ -z "$REF_PT" ]; then
  echo "WARNING: neither -r_mt nor -r_pt given. Round 1 of BOTH tracks will run with zero" >&2
  echo "         cross-organelle protection (see header comment: THE ROUND-1 BOOTSTRAP GAP)." >&2
  echo "         Protection begins in round 2, once each track has something from the other" >&2
  echo "         to screen against. Strongly consider supplying at least one of -r_mt/-r_pt" >&2
  echo "         if any reference (even a distantly related one) is available." >&2
fi

echo "=== DuOrganelle pipeline started: $(date) (mode=$MODE) ===" | tee "$OUTDIR/report/run.log"
echo "Read length: ${READ_LEN} bp; Stage B k-mers: $KMER_LIST; route: $ROUTE" | tee -a "$OUTDIR/report/run.log"

# Provenance: exact command and tool versions (report/versions.txt), for Methods.
first_ver() { "$@" 2>&1 | grep -m1 -oE '[vV]?[0-9]+\.[0-9]+(\.[0-9]+)*' || echo "n/a"; }
{
  echo "DuOrganelle      $DUORGANELLE_VERSION"
  echo "Command          $CMDLINE"
  echo "Working dir      $(pwd)"
  echo "Host             $(hostname 2>/dev/null || echo unknown)"
  echo "NOVOPlasty       4.3.5 (NOVOPlasty4.3.5.pl)"
  echo "GetOrganelle     $(first_ver get_organelle_from_reads.py --version)"
  echo "SPAdes           $(first_ver spades.py --version)"
  echo "BLAST+           $(first_ver blastn -version)"
  echo "minimap2         $(first_ver minimap2 --version)"
  echo "samtools         $(first_ver samtools --version)"
  echo "BBTools          $( (command -v bbversion.sh >/dev/null && first_ver bbversion.sh) || first_ver bbduk.sh --version)"
  echo "Python           $(first_ver python3 --version)"
} > "$OUTDIR/report/versions.txt" 2>/dev/null || true
tee -a "$OUTDIR/report/run.log" < "$OUTDIR/report/versions.txt"

ACC_MT="$OUTDIR/consensus/mt/accumulated_nonredundant.fasta"; : > "$ACC_MT"
ACC_PT="$OUTDIR/consensus/pt/accumulated_nonredundant.fasta"; : > "$ACC_PT"
CURRENT_SEED_MT="$SEED_MT"
CURRENT_SEED_PT="$SEED_PT"
MT_DONE=0
PT_DONE=0

# The anti-reference each track screens against: the user-supplied warm-start
# reference if given, falling back to the OTHER track's live accumulation
# (which starts empty in round 1 -- see bootstrap-gap warning above).
anti_ref_for_mt() { if [ -n "$REF_PT" ]; then echo "$REF_PT"; elif [ -s "$ACC_PT" ]; then echo "$ACC_PT"; else echo ""; fi; }
anti_ref_for_pt() { if [ -n "$REF_MT" ]; then echo "$REF_MT"; elif [ -s "$ACC_MT" ]; then echo "$ACC_MT"; else echo ""; fi; }

for round in $(seq 1 "$MAX_ROUNDS"); do
  # v2.0.1: the plastid track runs FIRST in every round. At ~10x mitochondrial
  # depth the plastome assembles easily from rbcL, and the mt track then starts
  # round 1 with that plastid sequence as its anti-reference (NOVOPlasty's
  # "Chloroplast sequence" field + the BLASTn screen). In v2.0.0 the mt track ran
  # first, extended into plastid sequence without any filter, and its round-1
  # output then caused every plastid contig of the same round to be rejected.
  if [ "$PT_DONE" -eq 0 ]; then
    echo "--- [pt] Round $round: NOVOPlasty seed=$CURRENT_SEED_PT ---" | tee -a "$OUTDIR/report/run.log"
    RDIR="$OUTDIR/novoplasty_runs/pt/round_${round}"; mkdir -p "$RDIR"
    ANTI_PT="$(anti_ref_for_pt)"
    # NOTE: --anti-reference is intentionally NOT passed to make_novoplasty_config.py
    # here -- NOVOPlasty's Type=chloro config has no field for it (verified against
    # NOVOPlasty's own config.txt; see make_novoplasty_config.py's docstring). The
    # pt track's cross-organelle protection comes entirely from the
    # redundancy_and_consensus.py call below, not from NOVOPlasty itself.
    if [ "$MODE" = "paired" ]; then
      python3 "$SCRIPT_DIR/make_novoplasty_config.py" \
        --seed "$CURRENT_SEED_PT" --target pt --mode paired --r1 "$R1" --r2 "$R2" \
        ${GENOME_RANGE_PT:+--genome-range "$GENOME_RANGE_PT"} \
        --out "$RDIR/config.txt" --project "pt_round${round}"
    else
      python3 "$SCRIPT_DIR/make_novoplasty_config.py" \
        --seed "$CURRENT_SEED_PT" --target pt --mode single --reads "$SEREADS" \
        ${GENOME_RANGE_PT:+--genome-range "$GENOME_RANGE_PT"} \
        --out "$RDIR/config.txt" --project "pt_round${round}"
    fi
    ( cd "$RDIR" && NOVOPlasty4.3.5.pl -c config.txt ) \
      > "$OUTDIR/round_logs/novoplasty_pt_round${round}.log" 2>&1 || {
        echo "WARNING: [pt] NOVOPlasty round $round failed; freezing pt track." | tee -a "$OUTDIR/report/run.log"
        PT_DONE=1
      }
    if [ "$PT_DONE" -eq 0 ]; then
      CONTIGS=$(find "$RDIR" -iname "Contigs_*.fasta" | head -1)
      if [ -z "$CONTIGS" ]; then
        echo "[pt] No contigs in round $round; freezing pt track." | tee -a "$OUTDIR/report/run.log"
        PT_DONE=1
      else
        python3 "$SCRIPT_DIR/redundancy_and_consensus.py" \
          --contigs "$CONTIGS" --accumulated "$ACC_PT" \
          --outdir "$OUTDIR/consensus/pt/round_${round}" \
          --identity "$BLAST_IDENTITY" --minlen "$BLAST_MINLEN" --threads "$THREADS" \
          ${ANTI_PT:+--anti-reference "$ANTI_PT" --other-frac-threshold "$OTHER_FRAC_THRESHOLD" --raw-other-frac-threshold "$RAW_OTHER_FRAC_THRESHOLD"}
        NEW_SEG="$OUTDIR/consensus/pt/round_${round}/novel_segment.fasta"
        NEW_BP=$(python3 -c "
s=''
try:
    with open('$NEW_SEG') as f:
        for l in f:
            if not l.startswith('>'): s+=l.strip()
except FileNotFoundError: pass
print(len(s))")
        echo "[pt] Round $round: ${NEW_BP} bp novel." | tee -a "$OUTDIR/report/run.log"
        cat "$NEW_SEG" >> "$ACC_PT" 2>/dev/null || true
        if [ "$NEW_BP" -lt "$MIN_NOVEL_BP" ]; then
          echo "[pt] Novel yield below threshold; pt track converged after round $round." | tee -a "$OUTDIR/report/run.log"
          PT_DONE=1
        else
          CURRENT_SEED_PT="$NEW_SEG"
        fi
      fi
    fi
  fi

  if [ "$MT_DONE" -eq 0 ]; then
    echo "--- [mt] Round $round: NOVOPlasty seed=$CURRENT_SEED_MT ---" | tee -a "$OUTDIR/report/run.log"
    RDIR="$OUTDIR/novoplasty_runs/mt/round_${round}"; mkdir -p "$RDIR"
    ANTI_MT="$(anti_ref_for_mt)"
    if [ "$MODE" = "paired" ]; then
      python3 "$SCRIPT_DIR/make_novoplasty_config.py" \
        --seed "$CURRENT_SEED_MT" --target mt --mode paired --r1 "$R1" --r2 "$R2" \
        ${ANTI_MT:+--anti-reference "$ANTI_MT"} \
        ${GENOME_RANGE_MT:+--genome-range "$GENOME_RANGE_MT"} \
        --out "$RDIR/config.txt" --project "mt_round${round}"
    else
      python3 "$SCRIPT_DIR/make_novoplasty_config.py" \
        --seed "$CURRENT_SEED_MT" --target mt --mode single --reads "$SEREADS" \
        ${ANTI_MT:+--anti-reference "$ANTI_MT"} \
        ${GENOME_RANGE_MT:+--genome-range "$GENOME_RANGE_MT"} \
        --out "$RDIR/config.txt" --project "mt_round${round}"
    fi
    ( cd "$RDIR" && NOVOPlasty4.3.5.pl -c config.txt ) \
      > "$OUTDIR/round_logs/novoplasty_mt_round${round}.log" 2>&1 || {
        echo "WARNING: [mt] NOVOPlasty round $round failed; freezing mt track." | tee -a "$OUTDIR/report/run.log"
        MT_DONE=1
      }
    if [ "$MT_DONE" -eq 0 ]; then
      CONTIGS=$(find "$RDIR" -iname "Contigs_*.fasta" | head -1)
      if [ -z "$CONTIGS" ]; then
        echo "[mt] No contigs in round $round; freezing mt track." | tee -a "$OUTDIR/report/run.log"
        MT_DONE=1
      else
        python3 "$SCRIPT_DIR/redundancy_and_consensus.py" \
          --contigs "$CONTIGS" --accumulated "$ACC_MT" \
          --outdir "$OUTDIR/consensus/mt/round_${round}" \
          --identity "$BLAST_IDENTITY" --minlen "$BLAST_MINLEN" --threads "$THREADS" \
          ${ANTI_MT:+--anti-reference "$ANTI_MT" --other-frac-threshold "$OTHER_FRAC_THRESHOLD" --raw-other-frac-threshold "$RAW_OTHER_FRAC_THRESHOLD"}
        NEW_SEG="$OUTDIR/consensus/mt/round_${round}/novel_segment.fasta"
        NEW_BP=$(python3 -c "
s=''
try:
    with open('$NEW_SEG') as f:
        for l in f:
            if not l.startswith('>'): s+=l.strip()
except FileNotFoundError: pass
print(len(s))")
        echo "[mt] Round $round: ${NEW_BP} bp novel." | tee -a "$OUTDIR/report/run.log"
        cat "$NEW_SEG" >> "$ACC_MT" 2>/dev/null || true
        if [ "$NEW_BP" -lt "$MIN_NOVEL_BP" ]; then
          echo "[mt] Novel yield below threshold; mt track converged after round $round." | tee -a "$OUTDIR/report/run.log"
          MT_DONE=1
        else
          CURRENT_SEED_MT="$NEW_SEG"
        fi
      fi
    fi
  fi

  if [ "$MT_DONE" -eq 1 ] && [ "$PT_DONE" -eq 1 ]; then
    echo "Both tracks converged after round $round." | tee -a "$OUTDIR/report/run.log"
    break
  fi
done

echo "mt accumulated: $ACC_MT ($(grep -c '^>' "$ACC_MT" 2>/dev/null || true) segments)" | tee -a "$OUTDIR/report/run.log"
echo "pt accumulated: $ACC_PT ($(grep -c '^>' "$ACC_PT" 2>/dev/null || true) segments)" | tee -a "$OUTDIR/report/run.log"

# -----------------------------------------------------------------------
# STAGE A2: coverage profile and Stage B route decision (coverage_profile.py).
# One minimap2 pass of the raw reads against plastid + mitochondrial references,
# giving median depth per compartment, organellar read fraction and insert
# size. References: pt = -r_pt > Stage A pt accumulation > -s_pt;
# mt = -r_mt > Stage A mt accumulation > -s_mt (v2.0.1; the median is robust to
# a minority of runaway extension, and -g_mt limits it). The plastid reference also drives the graph-route depletion.
# -----------------------------------------------------------------------
echo "--- Stage A2: coverage profile and Stage B route ---" | tee -a "$OUTDIR/report/run.log"
PROFDIR="$OUTDIR/profile"
STAGEB_FORCED_GRAPH=0; [ "$ROUTE" = "graph" ] && STAGEB_FORCED_GRAPH=1
if [ -n "$REF_PT" ]; then PT_REF_ALL="$REF_PT"; elif [ -s "$ACC_PT" ]; then PT_REF_ALL="$ACC_PT"; else PT_REF_ALL="$SEED_PT"; fi
# v2.0.1: without -r_mt, use the Stage A mt accumulation (now screened against
# the plastid track from round 1) rather than the seed genes alone, so the
# median mt depth is taken over hundreds of kb instead of one ~2-kb gene.
if [ -n "$REF_MT" ]; then MT_REF_PROF="$REF_MT"; elif [ -s "$ACC_MT" ]; then MT_REF_PROF="$ACC_MT"; else MT_REF_PROF="$SEED_MT"; fi
{ awk '/^>/{sub(/^>/,">pt__")}1' "$PT_REF_ALL"; awk '/^>/{sub(/^>/,">mt__")}1' "$MT_REF_PROF"; } > "$PROFDIR/ref_org.fa"

PROFILE_OK=0
# minimap2 (-ax sr) streamed into sam_profile.py, which writes BBMap-style
# covstats/ihist/log files for coverage_profile.py. Replaced BBMap here
# because BBMap 40.05 stalled indefinitely in index construction for this
# small reference on svrlinux01 (>19 h, and again >10 min, at "Indexing
# threads started"). Time-limited; with --route graph a failed profile aborts.
if [ "$MODE" = "paired" ]; then MM_IN=("$R1" "$R2"); else MM_IN=("$SEREADS"); fi
if timeout "$PROFILE_TIMEOUT" minimap2 -ax sr -t "$THREADS" --secondary=no "$PROFDIR/ref_org.fa" "${MM_IN[@]}" \
     2> "$PROFDIR/minimap2.log" \
   | python3 "$SCRIPT_DIR/sam_profile.py" --ref "$PROFDIR/ref_org.fa" \
       --covstats "$PROFDIR/covstats.txt" --ihist "$PROFDIR/ihist.txt" --log "$PROFDIR/map.log" \
       2>> "$PROFDIR/minimap2.log"; then
  tail -n 1 "$PROFDIR/minimap2.log" | tee -a "$OUTDIR/report/run.log"
  python3 "$SCRIPT_DIR/coverage_profile.py" --covstats "$PROFDIR/covstats.txt" \
    --ihist "$PROFDIR/ihist.txt" --bbmap-log "$PROFDIR/map.log" --sample-fraction 1 \
    --target-ratio "$PT_TARGET_RATIO" --enrich-threshold "$ENRICH_THRESHOLD" --route "$ROUTE" \
    --out-env "$PROFDIR/profile.env" | tee -a "$OUTDIR/report/run.log"
  PROFILE_OK=1
else
  echo "WARNING: minimap2 profiling failed or exceeded $PROFILE_TIMEOUT (see $PROFDIR/minimap2.log)." | tee -a "$OUTDIR/report/run.log"
  if [ "$STAGEB_FORCED_GRAPH" = "1" ]; then
    echo "ERROR: --route graph needs the profile for depletion; aborting instead of assembling all reads." | tee -a "$OUTDIR/report/run.log"
    exit 1
  fi
fi
if [ "$PROFILE_OK" -eq 0 ]; then
  FALLBACK_ROUTE="$ROUTE"; [ "$ROUTE" = "auto" ] && FALLBACK_ROUTE="recruit"
  printf 'COV_PT=0\nCOV_MT=0\nORG_FRAC=0\nINSERT_MEDIAN=0\nSAMPLERATE=1\nSTAGEB_ROUTE=%s\n' "$FALLBACK_ROUTE" > "$PROFDIR/profile.env"
  echo "         Stage B route: $FALLBACK_ROUTE (no depletion; repeat spannability unknown)." | tee -a "$OUTDIR/report/run.log"
fi
# shellcheck disable=SC1090
source "$PROFDIR/profile.env"
if [ "$STAGEB_FORCED_GRAPH" = "1" ] && awk -v c="$COV_MT" 'BEGIN{exit !(c<=0)}'; then
  echo "ERROR: --route graph but no mitochondrial depth was measured (COV_MT=0), so the plastid" | tee -a "$OUTDIR/report/run.log"
  echo "       depletion rate cannot be set. Check $PROFDIR/covstats.txt and the -r_mt/-s_mt files." | tee -a "$OUTDIR/report/run.log"
  exit 1
fi

# -----------------------------------------------------------------------
# STAGE B: single joint GetOrganelle run, both organelles at once, informed
# by BOTH accumulated seed sets -- lets GetOrganelle's own gene-database-based
# read classification (Jin et al. 2020) decide with more context than either
# single-organelle mode or a static reference could give it.
# -----------------------------------------------------------------------
echo "--- Stage B: organelle assembly (route: $STAGEB_ROUTE) ---" | tee -a "$OUTDIR/report/run.log"

# If a track's Stage A accumulation ended up EMPTY (e.g. its seed never
# extended, or every raw contig got screened out as cross-organelle
# contamination), passing that empty file straight to
# get_organelle_from_reads.py's -s produced an unusable seed list in earlier
# testing. Guard against it explicitly and fall back to that track's
# ORIGINAL seed so Stage B still gets something usable to seed from.
SEED_ARG_MT="$ACC_MT"
# With -r 0 (Stage A skipped on purpose) empty accumulations are expected;
# fall back to the seeds silently.
if [ "$MAX_ROUNDS" -eq 0 ]; then
  [ -s "$ACC_MT" ] || SEED_ARG_MT="$SEED_MT"
elif [ ! -s "$ACC_MT" ]; then
  echo "WARNING: mt accumulated set is EMPTY after Stage A (0 bp across all rounds)." | tee -a "$OUTDIR/report/run.log"
  echo "         Falling back to the original -s_mt seed ($SEED_MT) for Stage B instead." | tee -a "$OUTDIR/report/run.log"
  echo "         This usually means NOVOPlasty's mt track never produced usable contigs, or every" | tee -a "$OUTDIR/report/run.log"
  echo "         raw contig was screened out as cross-organelle contamination -- inspect" | tee -a "$OUTDIR/report/run.log"
  echo "         $OUTDIR/novoplasty_runs/mt/round_*/ and $OUTDIR/consensus/mt/round_*/*.tsv before trusting Stage B's result." | tee -a "$OUTDIR/report/run.log"
  SEED_ARG_MT="$SEED_MT"
fi
SEED_ARG_PT="$ACC_PT"
if [ "$MAX_ROUNDS" -eq 0 ]; then
  [ -s "$ACC_PT" ] || SEED_ARG_PT="$SEED_PT"
elif [ ! -s "$ACC_PT" ]; then
  echo "WARNING: pt accumulated set is EMPTY after Stage A (0 bp across all rounds)." | tee -a "$OUTDIR/report/run.log"
  echo "         Falling back to the original -s_pt seed ($SEED_PT) for Stage B instead." | tee -a "$OUTDIR/report/run.log"
  SEED_ARG_PT="$SEED_PT"
fi

if [ "$STAGEB_ROUTE" = "recruit" ]; then
  # ---- RECRUIT route (total-DNA library): GetOrganelle seed-and-extend ----
  echo "[Stage B] route RECRUIT: get_organelle_from_reads.py (-F embplant_mt,embplant_pt)" | tee -a "$OUTDIR/report/run.log"
  if [ "$MODE" = "paired" ]; then
    get_organelle_from_reads.py \
      -1 "$R1" -2 "$R2" \
      -s "$SEED_ARG_MT,$SEED_ARG_PT" \
      -F embplant_mt,embplant_pt -R 30 -t "$THREADS" -k "$KMER_LIST" \
      -o "$OUTDIR/getorganelle" --overwrite \
      > "$OUTDIR/round_logs/getorganelle.log" 2>&1 || true
  else
    get_organelle_from_reads.py \
      -u "$SEREADS" \
      -s "$SEED_ARG_MT,$SEED_ARG_PT" \
      -F embplant_mt,embplant_pt -R 30 -t "$THREADS" -k "$KMER_LIST" \
      -o "$OUTDIR/getorganelle" --overwrite \
      > "$OUTDIR/round_logs/getorganelle.log" 2>&1 || true
  fi

  # Check the output directory actually exists BEFORE globbing it -- `find` on a
  # missing directory produces an uninformative error ("No such file or
  # directory") instead of pointing at GetOrganelle's own log. Show the log's
  # tail directly so the real reason is visible immediately.
  if [ ! -d "$OUTDIR/getorganelle" ]; then
    echo "ERROR: GetOrganelle did not create $OUTDIR/getorganelle -- it failed before producing any output." | tee -a "$OUTDIR/report/run.log"
    echo "Last 40 lines of $OUTDIR/round_logs/getorganelle.log:" | tee -a "$OUTDIR/report/run.log"
    tail -n 40 "$OUTDIR/round_logs/getorganelle.log" | tee -a "$OUTDIR/report/run.log"
    exit 1
  fi

  # Confirmed on the O. villaricae run: dual -F embplant_mt,embplant_pt mode
  # does prefix outputs per organelle type, analogous to single-organelle
  # mode's "embplant_mt.K*...gfa".
  GFA_MT=$(find "$OUTDIR/getorganelle" -iname "*embplant_mt*.gfa" 2>/dev/null | head -1 || true)
  GFA_PT=$(find "$OUTDIR/getorganelle" -iname "*embplant_pt*.gfa" 2>/dev/null | head -1 || true)
  LOG_MT="$OUTDIR/round_logs/getorganelle.log"
  LOG_PT="$OUTDIR/round_logs/getorganelle.log"

else
  # ---- GRAPH route (organelle-enriched library) ----
  # 1) Partial plastid depletion (BBDuk + reformat.sh): pairs whose BOTH mates
  #    are >= 90 % covered by plastome 31-mers are classified plastid
  #    (rieb=f); only SAMPLERATE of them is kept, bringing plastid depth to
  #    ~PT_TARGET_RATIO x mito depth. All other pairs -- including those that
  #    span MTPT/mito junctions -- are kept, so MTPTs identical to the
  #    plastome are not turned into gaps (which full removal would do).
  # 2) SPAdes on the whole (depleted) read set.
  # 3) get_organelle_from_assembly.py, once per organelle, on the SPAdes graph.
  echo "[Stage B] route GRAPH: partial plastid depletion + SPAdes + get_organelle_from_assembly.py" | tee -a "$OUTDIR/report/run.log"
  DEPLDIR="$OUTDIR/depletion"; mkdir -p "$DEPLDIR"
  ASM_R1="${R1:-}"; ASM_R2="${R2:-}"; ASM_U="${SEREADS:-}"
  DO_DEPLETE=$(awk -v s="$SAMPLERATE" 'BEGIN{print (s>0 && s<1) ? 1 : 0}')
  if [ "$DO_DEPLETE" -eq 1 ] && [ "$PT_REF_ALL" = "$SEED_PT" ]; then
    echo "WARNING: no plastome available for depletion (only -s_pt seed genes); skipping depletion." | tee -a "$OUTDIR/report/run.log"
    DO_DEPLETE=0
  fi
  if [ "$DO_DEPLETE" -eq 1 ]; then
    echo "[Stage B] depleting plastid pairs to samplerate=$SAMPLERATE (reference: $PT_REF_ALL)" | tee -a "$OUTDIR/report/run.log"
    if [ "$MODE" = "paired" ]; then
      timeout "$PROFILE_TIMEOUT" bbduk.sh in1="$R1" in2="$R2" ref="$PT_REF_ALL" k=31 mcf=0.9 rieb=f threads="$THREADS" \
        outu1="$DEPLDIR/nopt_1.fq.gz" outu2="$DEPLDIR/nopt_2.fq.gz" \
        outm1="$DEPLDIR/pt_1.fq.gz" outm2="$DEPLDIR/pt_2.fq.gz" > "$DEPLDIR/bbduk.log" 2>&1 || {
          echo "ERROR: BBDuk depletion failed or exceeded $PROFILE_TIMEOUT -- see $DEPLDIR/bbduk.log" | tee -a "$OUTDIR/report/run.log"
          exit 1; }
      reformat.sh in1="$DEPLDIR/pt_1.fq.gz" in2="$DEPLDIR/pt_2.fq.gz" \
        out1="$DEPLDIR/ptS_1.fq.gz" out2="$DEPLDIR/ptS_2.fq.gz" \
        samplerate="$SAMPLERATE" sampleseed=13 > "$DEPLDIR/reformat.log" 2>&1
      cat "$DEPLDIR/nopt_1.fq.gz" "$DEPLDIR/ptS_1.fq.gz" > "$DEPLDIR/mix_1.fq.gz"
      cat "$DEPLDIR/nopt_2.fq.gz" "$DEPLDIR/ptS_2.fq.gz" > "$DEPLDIR/mix_2.fq.gz"
      ASM_R1="$DEPLDIR/mix_1.fq.gz"; ASM_R2="$DEPLDIR/mix_2.fq.gz"
    else
      bbduk.sh in="$SEREADS" ref="$PT_REF_ALL" k=31 mcf=0.9 threads="$THREADS" \
        out="$DEPLDIR/nopt.fq.gz" outm="$DEPLDIR/pt.fq.gz" > "$DEPLDIR/bbduk.log" 2>&1
      reformat.sh in="$DEPLDIR/pt.fq.gz" out="$DEPLDIR/ptS.fq.gz" \
        samplerate="$SAMPLERATE" sampleseed=13 > "$DEPLDIR/reformat.log" 2>&1
      cat "$DEPLDIR/nopt.fq.gz" "$DEPLDIR/ptS.fq.gz" > "$DEPLDIR/mix.fq.gz"
      ASM_U="$DEPLDIR/mix.fq.gz"
    fi
    grep -E "Contaminants|Result" "$DEPLDIR/bbduk.log" | sed 's/^/[bbduk] /' | tee -a "$OUTDIR/report/run.log" || true
  else
    echo "[Stage B] no depletion (samplerate=$SAMPLERATE); assembling all reads." | tee -a "$OUTDIR/report/run.log"
  fi

  echo "[Stage B] SPAdes (k=$KMER_LIST, -m $SPADES_MEM)" | tee -a "$OUTDIR/report/run.log"
  if [ "$MODE" = "paired" ]; then
    SPADES_IN=(-1 "$ASM_R1" -2 "$ASM_R2")
  else
    SPADES_IN=(-s "$ASM_U")
  fi
  FASTG="$OUTDIR/spades/assembly_graph.fastg"
  # Resume: if a finished SPAdes graph already exists in this output folder,
  # reuse it (SPAdes is the expensive step; a failure after it should not
  # force a multi-hour re-assembly).
  if [ -s "$FASTG" ] && grep -q "SPAdes pipeline finished" "$OUTDIR/round_logs/spades.log" 2>/dev/null; then
    echo "[Stage B] reusing existing SPAdes graph: $FASTG" | tee -a "$OUTDIR/report/run.log"
  else
  spades.py "${SPADES_IN[@]}" -k "$KMER_LIST" -t "$THREADS" -m "$SPADES_MEM" \
    -o "$OUTDIR/spades" > "$OUTDIR/round_logs/spades.log" 2>&1 || {
      echo "ERROR: SPAdes failed. Last 40 lines of $OUTDIR/round_logs/spades.log:" | tee -a "$OUTDIR/report/run.log"
      tail -n 40 "$OUTDIR/round_logs/spades.log" | tee -a "$OUTDIR/report/run.log"
      exit 1
    }
  fi
  [ -s "$FASTG" ] || { echo "ERROR: $FASTG missing or empty." | tee -a "$OUTDIR/report/run.log"; exit 1; }

  # get_organelle_from_assembly.py creates -o with os.mkdir (not makedirs),
  # so the parent folder must already exist.
  mkdir -p "$OUTDIR/getorganelle"
  for ORG in mt pt; do
    GODIR="$OUTDIR/getorganelle/$ORG"
    get_organelle_from_assembly.py -g "$FASTG" -F "embplant_$ORG" -t "$THREADS" \
      -o "$GODIR" --overwrite > "$OUTDIR/round_logs/getorganelle_from_assembly_$ORG.log" 2>&1 || true
  done
  GFA_MT=$(find "$OUTDIR/getorganelle/mt" -iname "*selected_graph*.gfa" 2>/dev/null | head -1 || true)
  [ -z "$GFA_MT" ] && GFA_MT=$(find "$OUTDIR/getorganelle/mt" -iname "*.gfa" 2>/dev/null | head -1 || true)
  GFA_PT=$(find "$OUTDIR/getorganelle/pt" -iname "*selected_graph*.gfa" 2>/dev/null | head -1 || true)
  [ -z "$GFA_PT" ] && GFA_PT=$(find "$OUTDIR/getorganelle/pt" -iname "*.gfa" 2>/dev/null | head -1 || true)
  LOG_MT="$OUTDIR/getorganelle/mt/get_org.log.txt"
  [ -s "$LOG_MT" ] || LOG_MT="$OUTDIR/round_logs/getorganelle_from_assembly_mt.log"
  LOG_PT="$OUTDIR/getorganelle/pt/get_org.log.txt"
  [ -s "$LOG_PT" ] || LOG_PT="$OUTDIR/round_logs/getorganelle_from_assembly_pt.log"
fi

if [ -z "$GFA_MT" ] && [ -z "$GFA_PT" ]; then
  echo "ERROR: no .gfa found for either organelle." | tee -a "$OUTDIR/report/run.log"
  for lg in "$OUTDIR"/round_logs/getorganelle*.log; do
    [ -s "$lg" ] || continue
    echo "--- last 20 lines of $lg ---" | tee -a "$OUTDIR/report/run.log"
    tail -n 20 "$lg" | tee -a "$OUTDIR/report/run.log"
  done
  echo "Fix the cause and rerun the SAME command with the same -o: a finished SPAdes graph is reused." | tee -a "$OUTDIR/report/run.log"
  exit 1
fi
[ -z "$GFA_MT" ] && echo "WARNING: no mt .gfa found (naming assumption may be wrong -- check getorganelle/ listing)." | tee -a "$OUTDIR/report/run.log"
[ -z "$GFA_PT" ] && echo "WARNING: no pt .gfa found (naming assumption may be wrong -- check getorganelle/ listing)." | tee -a "$OUTDIR/report/run.log"

# Plain-language Stage B status, one line per organelle, straight into
# report/run.log -- see classify_stage_b_result.py's docstring for why this
# was added. A .gfa existing (checked above) only means Stage B produced
# SOME output; it says nothing about whether that output is a validated
# single genome or a fragmented/flagged one, which is what this line is for.
if [ -n "$GFA_MT" ]; then
  python3 "$SCRIPT_DIR/classify_stage_b_result.py" --log "$LOG_MT" --organelle embplant_mt | tee -a "$OUTDIR/report/run.log"
fi
if [ -n "$GFA_PT" ]; then
  python3 "$SCRIPT_DIR/classify_stage_b_result.py" --log "$LOG_PT" --organelle embplant_pt | tee -a "$OUTDIR/report/run.log"
fi

# -----------------------------------------------------------------------
# STAGE C: graph-node screening AND classification for EACH organelle
# (graph_depth_screen.py, merged with the former Stage C2 classifier):
# depth vs a robust baseline, cross-organelle homology, and per-node class
# (mt / mt_repeat x n / mtpt / plastid / low_depth for the mt graph;
# pt / pt_repeat / pt_shared_with_mt / mitochondrial / low_depth for the pt
# graph), repeat spannability against read length and median insert,
# Bandage colours, and an estimated genome size.
# Anti-references: the most complete OTHER-organelle sequence available --
# mt graph: -r_pt > Stage A pt accumulation; pt graph: -r_mt > mt accumulation.
# -----------------------------------------------------------------------
if [ -n "$REF_PT" ]; then ANTI_C_MT="$REF_PT"; else ANTI_C_MT="$ACC_PT"; fi
if [ -n "$REF_MT" ]; then ANTI_C_PT="$REF_MT"; else ANTI_C_PT="$ACC_MT"; fi
if [ -n "$GFA_MT" ]; then
  echo "--- Stage C [mt]: graph node screening + classification (anti-reference = $ANTI_C_MT) ---" | tee -a "$OUTDIR/report/run.log"
  python3 "$SCRIPT_DIR/graph_depth_screen.py" --organelle mt \
    --gfa "$GFA_MT" --accumulated "$ACC_MT" --outdir "$OUTDIR/bandage/mt" --threads "$THREADS" \
    --insert "${INSERT_MEDIAN:-0}" --read-len "$READ_LEN" \
    ${ANTI_C_MT:+--anti-reference "$ANTI_C_MT" --other-frac-threshold "$OTHER_FRAC_THRESHOLD"} \
    | tee -a "$OUTDIR/report/run.log"
fi
if [ -n "$GFA_PT" ]; then
  echo "--- Stage C [pt]: graph node screening + classification (anti-reference = $ANTI_C_PT) ---" | tee -a "$OUTDIR/report/run.log"
  python3 "$SCRIPT_DIR/graph_depth_screen.py" --organelle pt \
    --gfa "$GFA_PT" --accumulated "$ACC_PT" --outdir "$OUTDIR/bandage/pt" --threads "$THREADS" \
    --insert "${INSERT_MEDIAN:-0}" --read-len "$READ_LEN" \
    ${ANTI_C_PT:+--anti-reference "$ANTI_C_PT" --other-frac-threshold "$OTHER_FRAC_THRESHOLD"} \
    | tee -a "$OUTDIR/report/run.log"
fi

# -----------------------------------------------------------------------
# STAGE C3: dead-end tracing (trace_dead_ends.py). Each dead end of the
# curated mt graph is followed through GetOrganelle's unfiltered graph
# (slimmed_assembly_graph.fastg) to report what it runs into (typically a
# plastid-identical MTPT cut by the depth filter) and candidate connections
# between dead ends through other-organelle sequence. Hypotheses only.
# -----------------------------------------------------------------------
SLIM_MT=$(find "$OUTDIR/getorganelle" -path "*mt*" -name "slimmed_assembly_graph.fastg" 2>/dev/null | head -1 || true)
[ -z "$SLIM_MT" ] && SLIM_MT=$(find "$OUTDIR/getorganelle" -name "*.fastg" -path "*embplant_mt*" 2>/dev/null | head -1 || true)
if [ -n "$GFA_MT" ] && [ -n "$SLIM_MT" ] && [ -s "$OUTDIR/bandage/mt/curated_graph.gfa" ]; then
  echo "--- Stage C3 [mt]: dead-end tracing through the unfiltered graph ---" | tee -a "$OUTDIR/report/run.log"
  python3 "$SCRIPT_DIR/trace_dead_ends.py" --gfa "$OUTDIR/bandage/mt/curated_graph.gfa" \
    --fastg "$SLIM_MT" --labels "${SLIM_MT%.fastg}.csv" \
    ${ANTI_C_MT:+--other-ref "$ANTI_C_MT"} --threads "$THREADS" \
    --outdir "$OUTDIR/bandage/mt" | tee -a "$OUTDIR/report/run.log" || \
    echo "WARNING: Stage C3 dead-end tracing failed." | tee -a "$OUTDIR/report/run.log"
fi

echo "=== DuOrganelle pipeline finished: $(date) ===" | tee -a "$OUTDIR/report/run.log"
echo "mt repeat/contamination table: $OUTDIR/bandage/mt/repeat_candidates.tsv (if produced)"
echo "pt repeat/contamination table: $OUTDIR/bandage/pt/repeat_candidates.tsv (if produced)"
echo "RECOMMENDED FINAL STEP (not automated here): run one last direct BLASTn of the"
echo "final mt assembly against the final pt assembly, bidirectionally, as an end-to-end"
echo "sanity check independent of the per-round screens above."
