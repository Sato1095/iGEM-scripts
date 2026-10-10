# Final report — file manifest

Project: Boolean model of MDM2 and KRAS degradation in NSCLC (A549).
All checksums are SHA-256 of the files as generated.

---

## TIER 1 — Main figures and tables

| file | role | sha256 (first 12) |
|---|---|---|
| `network_edges.csv` | **Table 1.** The model network: 29 nodes, 99 signed directed edges, each with sign, sign basis, mechanism, resources, PMIDs, curation effort and provenance. | 1e3c8bd2bc52 |
| `network_nodes.csv` | **Table 2.** Node list with role (input / intermediate / readout) and in/out degree. | 341c3d5f9e92 |
| `rules_structural.txt` | **Table 3.** Boolean rules. Every departure from the default convention is marked OVERRIDE with its reason inline. | beb357dc68e5 |
| `attractor_report.md` | **Table 4 / Figure 1.** Exact attractor per condition by value percolation. All 29 nodes determined in all four conditions; unique attractor provably reachable from all 2^26 states. | 2f8951d56012 |
| `screen_report.md` | **Table 5 / Figure 2.** Rule-variant screen: full 32-row grid, which rule forms are determined, failure modes, predictions of passing variants, discriminating experiment. | d7c608f241b2 |

## TIER 2 — Supplementary methods and data

| file | role | sha256 (first 12) |
|---|---|---|
| `model_scope_v1.txt` | **Supp. Table S1.** Declared node scope from literature, with module comments. The stated boundary of the model. | de1bc39496de |
| `decisions_final.csv` | **Supp. Table S2.** 25 manual edge adjudications: 5 sign calls, 20 drops, each with rationale and reference basis. | 6f60eab683d9 |
| `screen_criteria.md` | **Supp. Methods.** Pre-registered pass criteria, fixed before any variant was run, including the documented revisions (XIAP removal, R4 redefinition). | 7b7778d87cdf |
| `rules_justification.md` | **Supp. Methods.** Three-tier justification of every rule: mechanical, structural override, or fitted. Contains the explicit statement of which rules were data-influenced. | 038e95b20bdd |
| `validation_nutlin_A549.csv` | **Supp. Table S3.** The MDM2-inhibition proteomics mapped onto model nodes, with log2 fold change and directional calls. | b73b1021a262 |
| `screen_results.json` | **Supp. Data.** Machine-readable full screen output: every variant, every condition, full attractor, failed criteria. | 975e5d480d13 |
| `network_changelog.csv` | **Supp. Table S4.** All 28 changes from the 117-edge network to the curated network, each with original sign, action and reason. | 1cdfd5266397 |

## TIER 3 — Code (repository / Zenodo)

| file | role | sha256 (first 12) |
|---|---|---|
| `build_network.py` | **Step 1.** Network construction from OmniPath with provenance tracking, evidence thresholds and a run manifest. Produces `final/edges.csv` (31 nodes / 117 edges). | 1b360d1c67ba |
| `curate_network.py` | **Step 2.** Applies the documented literature corrections to produce the curated network (29 nodes / 99 edges) plus `network_changelog.csv`. Verified to reproduce `network_edges.csv` exactly. | 30d4c0ad242e |
| `attractor_analysis.py` | **Step 3.** Exact attractor analysis: value percolation + exhaustive terminal-SCC search + random-start sampling for corroboration. | 357c8e6e5d51 |
| `run_screen.py` | **Step 4.** The 32-variant rule screen executing the pre-registered criteria; also regenerates `screen_report.md`. | 105640b77073 |
| `make_rules.py` | **Step 3.** Rule drafting, sign-consistency validation, and MaBoSS/BoolNet export. | c32e57a384af |

### Pipeline, end to end

Run `run_pipeline.ps1` (Windows) or `run_pipeline.sh` (POSIX) to execute all
of it in one command. The individual steps are:

```
build_network.py   --anchors KRAS MDM2 --sinks CASP3 CASP7 PARP1 CDKN1A \
                   --nodes model_scope_v1.txt --min-resources 2 \
                   --min-references 1 --min-curation-effort 5 \
                   --terminal-sinks --decisions decisions_final.csv -d final
curate_network.py --in final/edges.csv --out .
make_rules.py      --edges network_edges.csv --nodes network_nodes.csv \
                   --fixed KRAS MDM2 MDM4 --out rules_draft.txt
                   # hand-edited to rules_structural.txt; overrides documented
                   # inline and in rules_justification.md
attractor_analysis.py --rules rules_structural.txt --clamp ... --out attractor_report.md
run_screen.py      --rules rules_structural.txt \
                   --out-json screen_results.json --out-md screen_report.md
```

Only the rules step involves manual work; every other artifact regenerates
from its inputs.

**Also required in the repository, not generated here:**
- `run01/raw/` — the archived OmniPath snapshot (raw API responses with URL, timestamp and checksum per file). Without this the network cannot be regenerated, because OmniPath updates continuously.
- `final/manifest.json` — the run manifest for the 31-node build.

---

## Claims the files support

1. A provenance-tracked signed network for KRAS/MDM2 signalling and intrinsic
   apoptosis, with all manual curation documented.
2. Exact attractor analysis: one fixed point per condition, no cyclic
   attractors, outcome independent of initial state.
3. MDM2 degradation alone never produces apoptosis — 0 of 32 rule variants.
4. Combined degradation produces apoptosis in 28 of 32 variants.
5. True synergy (combination kills, neither single agent does) in 4 of 32,
   requiring both BCL-XL as a coincidence detector and p53-dependent PUMA.
6. The discriminating experiment: PUMA induction after KRAS degradation in
   p53-intact A549, or FOXO3 knockdown in the same background.

## Claims the files do NOT support

- That the model is validated. Consistency with one drug arm, using rules
  selected after inspecting model–data disagreements, is not validation.
- Anything about KRAS degradation. No dataset here perturbs KRAS, so that
  arm of the model is untested.
- Anything about the combination in cells. No dataset in this project perturbs
  both targets.

## Held back / outstanding

- Bortezomib and seven other drugs in the A549 proteomics remain unused — the
  only genuine out-of-sample test available.
- The in-silico single-node perturbation screen (29 nodes x 2 states x 4
  conditions) has not been run.
- Nutlin inhibits MDM2–p53 binding; the model clamps MDM2 to 0. These are not
  the same perturbation, and the clamped model cannot reproduce the observed
  rise in MDM2 protein.
