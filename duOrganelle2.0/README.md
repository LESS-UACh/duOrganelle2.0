# DuOrganelle

**Joint, reference-free assembly of the plastid and mitochondrial genomes of a plant from one short-read library.**

DuOrganelle assembles both organellar genomes from the same Illumina/DNBSEQ paired-end (or single-end) library. Each organelle is used as a contamination screen for the other. The screen is bidirectional and is updated every round, so neither genome has to be assembled beforehand.

It was developed for, and first used in:

> Soto-Ivanoff, A.; Angulo, F.E.; Aranda-Díaz, M.; Nuñez, J.J. Organellar Genomes of the Hummingbird Fuchsia, *Fuchsia magellanica* (Onagraceae). *Int. J. Mol. Sci.* **2026** (submitted).

## How it works

| Stage | What it does | Tools |
|---|---|---|
| **A1. Seed-guided assembly** | NOVOPlasty extends a plastid seed (e.g. *rbcL*) and then a mitochondrial seed (e.g. *matR*), over successive rounds; the plastid track runs first so that the mitochondrial track is filtered against it. After each round, the contigs of each organelle are screened by BLASTn against the current assembly of the other. A contig is discarded if ≥ 30 % of its length matches before clustering, or ≥ 50 % after redundant contigs are collapsed. The new sequence becomes the next round's seed. Iteration stops when a round adds < 500 bp. | NOVOPlasty 4.3.5, BLAST+, `make_novoplasty_config.py`, `redundancy_and_consensus.py` |
| **A2. Library profiling** | Maps all reads to the plastid sequence and the mitochondrial seed/reference. Reports median plastid and mitochondrial depth, organellar read fraction and median insert size, then chooses the Stage B route. | minimap2, `sam_profile.py`, `coverage_profile.py` |
| **B. Organelle assembly** | **RECRUIT** route (total-DNA libraries): GetOrganelle seed-and-extend on both organelles (`-F embplant_mt,embplant_pt`). **GRAPH** route (organelle-rich libraries, or forced with `--route graph`): partial plastid depletion (BBDuk + `reformat.sh`, keeping pairs that span MTPT junctions), global SPAdes assembly, then `get_organelle_from_assembly.py` once per organelle. | GetOrganelle, BBTools, SPAdes |
| **C1. Node classification** | Classifies every graph node by depth relative to the organelle's reference depth *D* and by homology to the other organelle. Classes: mt / repeat (*n* copies) / MTPT / plastid / low-depth (NUMT candidate). Flags repeats longer than the insert size, writes Bandage colours and labels, and estimates genome size. | `graph_depth_screen.py` |
| **C2. Dead-end tracing** | Follows each dead end of the curated mitochondrial graph through GetOrganelle's unfiltered graph and reports candidate MTPT-mediated connections. These are hypotheses only and are not used to join the assembly. | `trace_dead_ends.py` |

`classify_stage_b_result.py` turns GetOrganelle's log into one plain-language status line per organelle (circular / fragmented / flagged), so that a "converged" run is not silently treated as finished.

## Installation

```bash
git clone https://github.com/LESS-UACh/duOrganelle2.0.git
cd duOrganelle2.0
conda env create -f environment.yml
conda activate duorganelle
get_organelle_config.py -a embplant_pt,embplant_mt   # GetOrganelle databases, once
chmod +x duorganelle.sh *.py
./duorganelle.sh --version
```

The script checks for every external tool before starting. All helper scripts must stay in the same folder as `duorganelle.sh`.

## Usage

```bash
./duorganelle.sh -m paired -1 R1.fq.gz -2 R2.fq.gz \
    -s_mt seed_mt.fasta -s_pt seed_pt.fasta \
    -o outdir -t 16 [options]
```

| Option | Meaning |
|---|---|
| `-s_mt`, `-s_pt` | Seed sequences (e.g. *matR*/*cox1*/*nad1*; *rbcL*/*matK*/*psbA*). Required. |
| `-r_mt`, `-r_pt` | Optional warm-start reference for either organelle. Removes the round-1 gap (see Limitations). |
| `-g_mt`, `-g_pt MIN-MAX` | NOVOPlasty genome-size range. Strongly recommended when the approximate size is known. NOVOPlasty extends towards the upper bound rather than treating it as a cap. |
| `-r N` | Maximum Stage A rounds (default 5; `-r 0` skips Stage A). |
| `-k LIST` | k-mers for Stage B (default `21,45,65,85,105,115,127`). Values that are even, > 127 or ≥ read length are dropped automatically. |
| `--route auto\|recruit\|graph` | Stage B route (default `auto`, decided from the organellar read fraction). |
| `--enrich-threshold F` | Organellar fraction at or above which `auto` chooses the graph route (default 0.5). |
| `--pt-ratio F` | Plastid:mitochondrial depth kept after depletion (default 2.5). |
| `--spades-mem G` | SPAdes memory limit in GB (default 120). |
| `-p`, `-q` | Cross-organelle screening thresholds after/before clustering (defaults 0.5 / 0.3). |

If a run stops after SPAdes, rerun the **same command with the same `-o`**: a finished SPAdes graph is reused.

### Example (*Fuchsia magellanica*, as in the paper)

Reference-free run (seed genes only):

```bash
./duorganelle.sh -m paired -1 FM18725_1.fq.gz -2 FM18725_2.fq.gz \
    -s_mt Seed_matR_MZ934755_Ovillaricae.fasta \
    -s_pt Seed_rbcL_Fmagellanica_plastome.fasta \
    -g_mt 200000-450000 --route graph --spades-mem 200 \
    -o duorganelle_B -t 16
```

### Validation

v2.0.1 was run twice on the *F. magellanica* library (DNBSEQ, 2 × 150 bp, 30.5 M reads): from the seed genes alone (run B, command above), and with an earlier plastome and a preliminary mitochondrial assembly supplied as `-r_pt` and `-r_mt` (run A). Both runs reproduced the deposited assemblies exactly:

| | Published | Run A (references) | Run B (seeds only) |
|---|---|---|---|
| Median depth pt / mt | 2,090× / 220× | 2,090× / 220× | 2,079× / 221× |
| Organellar read fraction | 9.0% | 9.0% | 8.4% |
| Plastid subsampling rate | 0.263 | 0.263 | 0.266 |
| Mitochondrial graph | 19 nodes, 22 links, 399,916 bp | identical (19/19 nodes exact) | identical (19/19 nodes exact) |
| Plastome | 155,942 bp (PZ340948) | identical | identical |

Tool versions used: NOVOPlasty 4.3.5, GetOrganelle 1.7.7.1, SPAdes 4.3.0, BLAST+ 2.17.0, minimap2 2.31, SAMtools 1.24, BBTools 40.05, Python 3.14.7.

`validation/compare_runs.py` repeats the comparison for any set of run folders against the published graph and plastome in `validation/`.

## Output

```
outdir/
├── report/run.log              # full log, including Stage B status lines
├── report/versions.txt         # DuOrganelle version, exact command and tool versions
├── novoplasty_runs/{mt,pt}/    # Stage A1 rounds
├── consensus/{mt,pt}/          # per-round screening tables and accumulated sequence
├── profile/                    # Stage A2 depth, insert size and route decision (profile.env)
├── depletion/                  # graph route: BBDuk/reformat logs and depleted reads
├── spades/                     # graph route: global assembly
├── getorganelle/{mt,pt}/       # organellar graphs (.gfa) and GetOrganelle logs
└── bandage/{mt,pt}/            # curated graph, node classes, repeat table, Bandage colours, dead-end report
```

After the run, the script recommends one final bidirectional BLASTn of the final mitochondrial assembly against the final plastome.

## Utility

`extract_neighborhood.py` extracts the neighbourhood of selected nodes from a large SPAdes GFA so it can be inspected in Bandage. It is not called by the pipeline.

```bash
python3 extract_neighborhood.py --gfa outdir/spades/assembly_graph_with_scaffolds.gfa \
    --nodes 7564413,7570939 --distance 3 --out tips_neighbourhood
```

## Limitations

- **Round-1 gap.** Without `-r_pt`, the plastid track runs unprotected in round 1. Since v2.0.1 it runs first, so the mitochondrial track is screened against it from round 1 onward. Supplying `-r_pt` or `-r_mt` removes the gap entirely.
- **Short reads cannot resolve long repeats.** Repeats longer than the insert size cannot be spanned. A fragmented mitochondrial graph may be the correct graph-level result for a recombining, multipartite mitogenome. DuOrganelle reports this but cannot decide between biology and data limitation; long reads are needed for that.
- **Runaway extension.** With wide default genome ranges, NOVOPlasty can extend the mitochondrial track into repeats or nuclear insertions. Use `-g_mt` when a size estimate is available.
- Tested on three datasets so far: *Fuchsia magellanica* and *Oenothera villaricae* (Onagraceae), and *Punica granatum* (Lythraceae).

## Authors

Millaray Aranda-Díaz, Fernanda E. Angulo, Anaía Soto-Ivanoff and José J. Nuñez — Instituto de Ciencias Marinas y Limnológicas, Facultad de Ciencias, Universidad Austral de Chile (LESS-UACh).

## Acknowledgments

Parts of the code and documentation were prepared with the assistance of Claude (Anthropic). The authors reviewed and tested all code and take full responsibility for it.

## Citation

If you use DuOrganelle, please cite the paper above and the tools it wraps: NOVOPlasty (Dierckxsens et al. 2017), GetOrganelle (Jin et al. 2020), SPAdes (Prjibelski et al. 2020), BLAST+ (Camacho et al. 2009), minimap2 (Li 2018), SAMtools (Danecek et al. 2021) and BBTools (Bushnell).

## License

MIT, see [LICENSE](LICENSE).
