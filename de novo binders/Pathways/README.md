# Pathways — Boolean model of apoptotic commitment in A549

Asks one question: does degrading MDM2 and KRAS together kill a cell that
survives either alone?

A549 is TP53 wild-type, KRAS G12S, CDKN2A-deleted, which is what makes it the
right background — p53 is intact and available to be released, and the KRAS
arm is independently druggable. This is a separate strand from the structural
work; it shares the project's targets, not its methods.

---

## The pipeline, in order

Run it all with `run_pipeline.sh` (POSIX) or `run_pipeline.ps1` (Windows).
Python 3.8+, standard library only.

### 1 · Network construction — `build_network.py` → `final/`

Builds a signed, directed, provenance-tracked network from OmniPath. No
biology is hard-coded in the script: the scientific choices arrive as the
declared node scope, the evidence thresholds and the reviewed sign calls.

We ran it in **scope mode**: `model_scope_v1.txt` declares 31 nodes as the
model's stated boundary, and the script fetches every qualified edge among
them. The alternative — growing the network outward from the anchors — would
have made the boundary a side effect of search parameters rather than a
decision anyone could review.

**Evidence thresholds:** directed, ≥2 independent resources, ≥1 reference,
curation effort ≥5. Every rejected edge goes to an exclusion file with its
reason.

→ 203,001 interactions fetched, 28,146 qualified, **117 edges among 31 nodes**.

### 2 · Curation — `curate_network.py` → `network_*.csv`

Evidence thresholds measure how well an edge is attested, not whether it
points the right way. Forty edges were adjudicated by hand across steps 1 and
2 — 25 as reviewed sign calls inside step 1, 15 as corrections here — falling
into five recurring failure modes of an aggregated interactome:

| failure mode | example |
|---|---|
| reverse records | `APAF1→CYCS`: cytochrome c is released from mitochondria, not produced by APAF1 |
| signs contradicting their own annotations | `AKT1→FOXO3` recorded as activation, when the phosphorylation excludes FOXO3 from the nucleus |
| reciprocal sequestration | one binding event recorded both ways, which in Boolean form is a mutual-inhibition toggle |
| text-mining artefacts | records resting only on automated readers, which routinely swap enzyme and substrate |
| transcriptional bypasses | `TP53→CASP3`, letting p53 fire the executioner caspase without mitochondrial permeabilisation |

Also here: **PIK3CA and PIK3R1 merge into one PI3K node** (recorded as
mutually activating, which latches ON permanently in Boolean form), and
**PARP1 is removed** (with `CASP3→PARP1` gone as a readout out-edge, its only
remaining in-edges are transcriptional, so it cannot act as a cleavage
readout — apoptosis is read from CASP3/CASP7 instead).

**One edge was restored against the thresholds.** `KRAS→PI3K` is
single-resource in OmniPath and was rejected in step 1, which left the entire
PI3K/AKT arm disconnected from the perturbed node — the model could not have
addressed its own question. It is tagged `provenance=MANUAL_ADD` in the edge
list; that tag should stay visible wherever the network is shown.

→ **29 nodes, 99 edges**, with all 28 changes in `network_changelog.csv`.

### 3 · Rules — `make_rules.py` → `rules_structural.txt`

An edge list cannot say whether two activators combine with AND or OR. The
script drafts a starting point under one default convention — *any one
activator suffices, any one inhibitor suffices* — flags the nodes where that
is likely wrong, and then checks the hand-edited result back against the
network.

**This is the one step with human judgement in the middle of it.**
`rules_structural.txt` was written by hand from the draft and is an *input*
to the pipeline, not something it regenerates.

The check verifies, for every input of every rule, that the rule's actual
monotone behaviour matches the sign of the edge, and reports edges no rule
uses. 29 rules, 0 problems, **17 warnings** — all of the same kind, an edge
the rules deliberately leave unused. Each is accounted for by an `OVERRIDE`
comment in the rules file naming the edge and the reason. An eighteenth would
mean something had been dropped by accident; that is what the check is for.

Exports `model.bnd`/`.cfg` (MaBoSS) and `model.bnet` (BoolNet). Steps 4
and 5 parse the rules file directly and do not read either.

### 4 · Attractors — `attractor_analysis.py` → `attractor_report.md`

**Value percolation is exact, not sampled.** Clamp the perturbed inputs, then
repeatedly fix any node whose rule has become constant over every assignment
of its still-free inputs, and iterate.

In all four conditions every one of the 29 nodes fixes within **six rounds**.
Each condition therefore has a single provable fixed point: it is the only
attractor, and it is reachable from every one of the 2^26 = 67,108,864 states
of the free subspace. No cyclic attractors, and the outcome does not depend
on the initial state.

Updating is **asynchronous**, one node per step. Fixed points are identical
under either scheme; only cyclic attractors would differ, and there are none.
Random-start sampling was run alongside and agreed in every condition — it is
corroboration, and proves nothing on its own.

### 5 · Variant screen — `run_screen.py` → `screen_report.md`

Five rules are defensible in two forms each, so all 2^5 = 32 combinations are
scored against the pre-registered criteria, under all four conditions.

**The criteria were frozen before any variant was run.** That ordering is the
whole point, and it must never be described the other way round.
`screen_criteria.md` is the frozen document; `run_screen.py` transcribes
its C1 and C2 into code. Only the control and MDM2-degraded conditions
carry criteria; no data in hand perturbs KRAS, so `KRAS_degraded` and
`both_degraded` are predictions.

| | |
|---|---|
| satisfy all criteria | **4 of 32** |
| show true synergy | **4 of 32** |
| apoptosis under control | 0 of 32 |
| under MDM2 degradation alone | 0 of 32 |
| under KRAS degradation alone | 24 of 32 |
| under the combination | 28 of 32 |

The two fours are **overlapping but different sets**, sharing two members.
They are not the same result and should not be quoted as one.

All four passers need R1 form B and R5 form B. **R2 — the PUMA rule (BBC3) —
is left undetermined**: two passers carry each form, because FOXO3 is off in
every condition the criteria constrain, so the data cannot separate them. The
synergy argument turns on R2, and that it is unresolved is reported, not
smoothed over.

**R1 is reported, not resolved.** The screen needs form B, a coincidence
detector requiring input from both arms; the binding literature leans toward
form A, since BAD is the canonical high-affinity BCL-XL antagonist. That
tension is stated rather than settled.

---

## What the model does not claim

- **It is not validated.** Consistency with one drug arm, using rules
  selected after inspecting model–data disagreements, is a constraint, not
  evidence of correctness.
- **Nothing tests the KRAS arm.** No dataset in hand perturbs KRAS, so the
  `KRAS_degraded` and `both_degraded` outcomes are predictions rather than
  results, and neither condition contributes a pass criterion.
- **Bortezomib and seven other drugs in the panel are held back** as the only
  genuine out-of-sample test available.
- **Nutlin is not degradation.** Nutlin inhibits MDM2–p53 binding; the model
  clamps MDM2 to 0. MDM2 protein rose 2.84 log2 under nutlin precisely
  because p53 transcribes it — behaviour a clamped model cannot show.
- **BCL2L11, PMAIP1, MCL1 and MDM4 were undetected** and are unconstrained by
  any criterion.

---

## Files

| file | role |
|---|---|
| `model_scope_v1.txt` | the 31 declared nodes — the model's stated boundary |
| `decisions_final.csv` | 25 reviewed sign calls (5 signs, 20 drops), each with a rationale |
| `final/` | step 1 output, including `manifest.json` and the 117-edge list |
| `run01/raw/` | **the pinned OmniPath snapshot** |
| `network_{edges,nodes,changelog}.csv` | the model network, 29/99, and all 28 changes |
| `rules_structural.txt` | the 29 Boolean rules, overrides marked inline |
| `screen_criteria.md` | the pre-registered criteria |
| `rules_justification.md` | tier-by-tier justification of the rule set |
| `validation_nutlin_A549.csv` | the MDM2-inhibition proteomics, mapped onto model nodes |
| `attractor_report.md`, `screen_report.md`, `screen_results.json` | results |
| `model.{bnd,cfg,bnet}` | MaBoSS and BoolNet exports |
| `CHECKSUMS.md`, `FINAL_REPORT_MANIFEST.md`, `checksums.txt` | artefact ledgers |

`model_scope_v1.txt` keeps its version suffix on purpose: `final/manifest.json`
records that exact filename and its SHA-256 as an input to the committed build.
Renaming it would break that chain, and editing the manifest would falsify a
provenance record.

### The OmniPath snapshot is pinned — do not re-query

`run01/raw/` holds the raw API responses, each with its URL, fetch timestamp,
byte count and SHA-256 in a `.meta.json` beside it. `build_network.py` reuses
a cached response whenever one matches the request, so the network rebuilds
against the snapshot rather than against live OmniPath.

OmniPath updates continuously. Re-fetching would give a different network and
silently invalidate everything downstream. Do not pass `--refresh`, and do
not delete this directory.

Scope mode uses the `interactions` and `enzsub` responses. The six
`annotations_*` responses are from exploratory runs that discovered phenotype
sinks by annotation regex; they are kept as part of the archive but are not
read when the node set is declared.

---

## Known issues

The scripts, their outputs and the supporting documents are consistent with
each other. Two things are worth knowing.

- **`model.cfg` sets every node to 0**, which is what the MaBoSS exporter
  writes. If you open the model in MaBoSS, set the initial states to the
  control attractor given in `rules_justification.md` §6 rather than using
  the file as shipped.
- **Ten nodes disagree between the structural rules and the Nutlin arm**,
  and only two of them became screened rules. The other eight are listed in
  `rules_justification.md` §5 and are open — two of them (`GSK3B`, and
  `APAF1` once BCL-XL is required to survive) cannot be satisfied by any
  rule set built on this network, which is a limit of the network rather
  than of the rules.
