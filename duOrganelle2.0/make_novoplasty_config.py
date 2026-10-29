#!/usr/bin/env python3
"""Generate a NOVOPlasty config.txt for a given seed and ORGANELLE TARGET
(mitochondrial or plastid), supporting both paired-end and single-end reads.

Generalized from repeatMito's original make_novoplasty_config.py (which only
ever generated mito_plant/mito configs) to also generate Type=chloro configs,
as part of the DuOrganelle dual-organelle architecture.

IMPORTANT, verified directly against NOVOPlasty's own config.txt template
(https://github.com/ndierckx/NOVOPlasty/blob/master/config.txt) before writing
this, rather than assumed:
  - Type is exactly one of: chloro | mito | mito_plant  (there is no "pt" alias).
  - The "Chloroplast sequence" field is documented as "(only for mito_plant
    option)" -- it has NO meaning/effect for Type=chloro, and NOVOPlasty ships
    with no symmetric "Mitochondrial sequence" exclusion field for that
    direction. This is precisely the structural gap DuOrganelle's own BLASTn
    cross-screening (in redundancy_and_consensus.py) is meant to fill for the
    plastid track, since NOVOPlasty itself provides no native protection there.
  - NOVOPlasty's own documented default "Genome Range" differs by organelle:
    12000-20000 for mito, 120000-200000 for chloro. We default accordingly
    per --target rather than sharing one default across both tracks.

Paired-end usage (mitochondrial track):
    make_novoplasty_config.py --seed seed_mt.fasta --target mt --mode paired \
        --r1 R1.fq.gz --r2 R2.fq.gz --anti-reference current_pt_accumulated.fasta \
        --out config_mt.txt --project mt_round1

Paired-end usage (plastid track):
    make_novoplasty_config.py --seed seed_pt.fasta --target pt --mode paired \
        --r1 R1.fq.gz --r2 R2.fq.gz \
        --out config_pt.txt --project pt_round1
    (--anti-reference is accepted but silently ignored for --target pt, with a
    printed note -- see WARNING below. NOVOPlasty has no field for it there.)
"""
import argparse
import os
import sys

TEMPLATE_PE = """Project:
-----------------------
Project name          = {project}
Type                   = {novo_type}
Genome Range           = {genome_range}
K-mer                  = 39
Max memory             = 20
Extended log           = 0
Save assembled reads   = no
Seed Input             = {seed}
Extend seed directly   = {extend_directly}
Reference sequence     =
Variance detection     =
{chloroplast_line}
Dataset 1:
-----------------------
Read Length            = 150
Insert size             = 300
Platform                = illumina
Single/Paired           = PE
Combined reads          =
Forward reads           = {r1}
Reverse reads           = {r2}
Store Hash              =

Optional:
-----------------------
Insert size auto        = yes
Use Quality Scores      = no
"""

TEMPLATE_SE = """Project:
-----------------------
Project name          = {project}
Type                   = {novo_type}
Genome Range           = {genome_range}
K-mer                  = 39
Max memory             = 20
Extended log           = 0
Save assembled reads   = no
Seed Input             = {seed}
Extend seed directly   = {extend_directly}
Reference sequence     =
Variance detection     =
{chloroplast_line}
Dataset 1:
-----------------------
Read Length            = 150
Insert size             =
Platform                = illumina
Single/Paired           = SE
Combined reads          = {reads}
Forward reads           =
Reverse reads           =
Store Hash              =

Optional:
-----------------------
Insert size auto        = yes
Use Quality Scores      = no
"""

DEFAULT_GENOME_RANGE = {
    "mt": "12000-2000000",   # repeatMito's historical wide default for plant mito
    "pt": "120000-200000",   # NOVOPlasty's own documented chloro default
}


def seed_length(seed_path):
    seed_len = 0
    seq = []
    with open(seed_path) as f:
        for line in f:
            if line.startswith(">"):
                if seq:
                    seed_len = max(seed_len, len("".join(seq)))
                seq = []
            else:
                seq.append(line.strip())
    if seq:
        seed_len = max(seed_len, len("".join(seq)))
    return seed_len


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", required=True)
    ap.add_argument("--target", choices=["mt", "pt"], required=True,
                     help="Which organelle this config assembles: 'mt' (mitochondrial; NOVOPlasty "
                          "Type=mito_plant, or Type=mito with --animal) or 'pt' (plastid; Type=chloro).")
    ap.add_argument("--animal", action="store_true",
                     help="With --target mt, use NOVOPlasty Type=mito instead of mito_plant (animal "
                          "mitochondrial genome, no plastid to worry about). Not valid with --target pt.")
    ap.add_argument("--mode", choices=["paired", "single"], default="paired")
    ap.add_argument("--r1", help="Forward reads (paired mode)")
    ap.add_argument("--r2", help="Reverse reads (paired mode)")
    ap.add_argument("--reads", help="Combined/unpaired reads file (single mode)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--project", default="assembly")
    ap.add_argument("--genome-range", default=None,
                     help="Overrides the per-target default (mt: 12000-2000000, pt: 120000-200000).")
    ap.add_argument("--anti-reference", "--chloroplast", dest="anti_reference",
                     help="FASTA of the CURRENT best assembly of the OTHER organelle (dynamically "
                          "updated round to round in the DuOrganelle loop -- not a fixed external "
                          "reference). Only meaningful with --target mt: it populates NOVOPlasty's "
                          "native 'Chloroplast sequence' field (Type=mito_plant only), which makes "
                          "NOVOPlasty subtract plastid-matching reads before/while extending the "
                          "mitochondrial seed. Silently ignored (with a printed note) for --target pt, "
                          "since NOVOPlasty has no symmetric field to exclude mitochondrial reads while "
                          "assembling the plastid genome -- for that direction, protection comes only "
                          "from the BLASTn cross-screens in redundancy_and_consensus.py, applied "
                          "downstream of NOVOPlasty rather than inside it.")
    args = ap.parse_args()

    if args.animal and args.target == "pt":
        ap.error("--animal is not valid with --target pt (no plastid genome in animals)")

    if args.mode == "paired":
        if not args.r1 or not args.r2:
            ap.error("--mode paired requires both --r1 and --r2")
    else:
        if not args.reads:
            ap.error("--mode single requires --reads")

    seed_len = seed_length(args.seed)
    extend_directly = "yes" if seed_len > 5000 else "no"

    if args.target == "mt":
        novo_type = "mito" if args.animal else "mito_plant"
    else:
        novo_type = "chloro"

    genome_range = args.genome_range or DEFAULT_GENOME_RANGE[args.target]

    anti_ref_path = os.path.abspath(args.anti_reference) if args.anti_reference else ""
    if args.target == "pt" and anti_ref_path:
        print("NOTE: --anti-reference given with --target pt, but NOVOPlasty's config format has no "
              "field to exclude mitochondrial reads while assembling the plastid genome (verified "
              "against NOVOPlasty's own config.txt: 'Chloroplast sequence' is documented as valid "
              "'only for mito_plant option'). The value is NOT written into this config. Protection "
              "for the plastid track comes entirely from the downstream BLASTn cross-screen against "
              "the current mitochondrial accumulation in redundancy_and_consensus.py.",
              file=sys.stderr)
        chloroplast_line = ""
    elif args.target == "mt":
        if not anti_ref_path and not args.animal:
            print("WARNING: no --anti-reference given for a --target mt config. NOVOPlasty's "
                  "'Chloroplast sequence' field will be left blank, so plastid-matching reads will NOT "
                  "be filtered before/during mitochondrial seed extension for this round. In the "
                  "DuOrganelle loop this is expected/unavoidable in round 1 (before the plastid track "
                  "has produced anything to screen against) but should stop being true from round 2 "
                  "onward once the plastid accumulation is non-empty.", file=sys.stderr)
        chloroplast_line = f"Chloroplast sequence   = {anti_ref_path}"
    else:
        chloroplast_line = ""

    common = dict(
        project=args.project, seed=os.path.abspath(args.seed),
        extend_directly=extend_directly, novo_type=novo_type,
        genome_range=genome_range, chloroplast_line=chloroplast_line,
    )

    if args.mode == "paired":
        cfg = TEMPLATE_PE.format(r1=os.path.abspath(args.r1), r2=os.path.abspath(args.r2), **common)
    else:
        cfg = TEMPLATE_SE.format(reads=os.path.abspath(args.reads), **common)

    with open(args.out, "w") as f:
        f.write(cfg)

    print(f"Wrote {args.out} (target={args.target}, mode={args.mode}, seed_len={seed_len}bp, "
          f"Type={novo_type}, Genome Range={genome_range}, "
          f"anti-reference={'set: ' + anti_ref_path if (anti_ref_path and args.target=='mt') else 'N/A'})")


if __name__ == "__main__":
    main()
