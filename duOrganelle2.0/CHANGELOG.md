# Changelog

## 2.0.1 (2026-10-29) — version used in Soto-Ivanoff et al. 2026 (IJMS)
- Stage A: the plastid track now runs before the mitochondrial track in every round. In 2.0.0 the mt track ran first without a plastid filter, extended into the (~10x deeper) plastome, and its output then caused all plastid contigs of the same round to be rejected, leaving both accumulations unusable on the *Fuchsia magellanica* data. The mt track now starts round 1 with the plastid sequence as its anti-reference.
- Stage A2: without `-r_mt`, the mitochondrial depth is measured on the Stage A mt accumulation instead of the seed genes alone.
- Validated on *Fuchsia magellanica*: reference-free and reference-supplied runs give identical assemblies, identical to the deposited ones (README, "Validation").
- `validation/compare_runs.py`: compares run folders with a published graph and plastome.

## 2.0.0 (2026-10-28)
- Stage A2 library profiling (minimap2 + `sam_profile.py` + `coverage_profile.py`) chooses the Stage B route (RECRUIT or GRAPH). minimap2 replaced BBMap here because BBMap 40.05 stalled indefinitely while indexing the small reference.
- GRAPH route for organelle-rich libraries: partial plastid depletion that keeps MTPT-spanning pairs, global SPAdes assembly, `get_organelle_from_assembly.py` per organelle.
- Stage C node classification merged into `graph_depth_screen.py` (mt / repeat / MTPT / plastid / low-depth), with repeat spannability, Bandage colours and genome-size estimate.
- Stage C2 dead-end tracing (`trace_dead_ends.py`).
- `classify_stage_b_result.py`: plain-language status for each organelle after Stage B.
- New options: `-g_mt`/`-g_pt` (NOVOPlasty genome range), `--route`, `--enrich-threshold`, `--pt-ratio`, `--spades-mem`, `--version`.
- k-mer list cleaned automatically (odd, ≤ 127, < read length).
- Input files and external tools are checked before anything starts.
- A finished SPAdes graph is reused when a run is restarted with the same `-o`.
- `report/versions.txt` records the exact command and tool versions.

## 1.x
- Iterative, bidirectionally screened NOVOPlasty assembly (Stage A) followed by a joint GetOrganelle run (`-F embplant_mt,embplant_pt`).
